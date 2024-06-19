# This code is used to generate the stimuli for the connected graph experiments


using ArgParse        # Argument parsing
using Arrow           # Compressed storage
using CSV
using DelimitedFiles
using Distributions
using DataFrames
using Random

include("generative_helpers.jl")
include("defaults.jl")

#SPACING_FOLDER_NAME = Dict("even" => "evenly", "expectation" => "expectation")

arg_settings = ArgParseSettings()
@add_arg_table arg_settings begin
    # Options
    "--seeds"
        help = "Seeds to run"
        action = :store_arg
        nargs = '*'
        arg_type = Int64

    "--users-per-pair"
        help = "Number of users per pair"
        action = :store_arg
        nargs = '*'
        arg_type = Int64

     "--pair-counts"
        help = "Pair counts to run"
        action = :store_arg
        nargs = '*'
        arg_type = Int64

    "--width"
        help="Width to run"
        arg_type = Int64
        default = 100
        range_tester = (x -> x > 0)

    "--stimuli-count"
        help="Stimuli count to run"
        arg_type = Int64
        default = 200
        range_tester = (x -> x > 2)

    "--spacing"
        help="Spacing mechanism (even or expectation)"
        arg_type = String
        default = "even"
        range_tester = (x -> x in ["even", "expectation"])

    # IGNORE THIS
    "--reliability"
        help="User reliability (provide directory)"
        arg_type = String

    "--format"
        help="Format for the preferences file (arrow or csv)"
        arg_type = String
        default = "arrow"
        range_tester = (x -> x in ["arrow", "csv"])
end

parsed_args = parse_args(ARGS, arg_settings)

stimuli_count = parsed_args["stimuli-count"]
width = parsed_args["width"]
spacing = parsed_args["spacing"]

format = parsed_args["format"] 
extension = format == "csv" ? "csv" : "arrow.lz4"

seeds = parsed_args["seeds"]
if length(seeds) < 1
    seeds = DEFAULT_SEEDS
end

users_per_pairs = parsed_args["users-per-pair"]
if length(users_per_pairs) < 1
    users_per_pairs = DEFAULT_USERS_PER_PAIRS
end

pair_counts = parsed_args["pair-counts"]
if length(pair_counts) < 1
    pair_counts = DEFAULT_PAIR_COUNTS
end

reliability = parsed_args["reliability"]

# Check that the pair counts can create a connected graph, even if it is degenerate
invalid_pair_counts = false
for x in pair_counts
    if x < stimuli_count
        @error("$(x) < $(stimuli_count).  There must be at least $(stimuli_count - 1) pairs requested.")
        global invalid_pair_counts = true
    end
end
if invalid_pair_counts
    exit()
end

output_directory = "data/same/$(spacing)/$(width)" 
mkpath(output_directory)

# Get the number of available pairs
pairs_per_attribute = round(Int, ((stimuli_count * (stimuli_count - 1)) / 2))

if spacing == "even"
    w = width / 100
    stimuli_distances = Dict((i) => i * w for i in 1:stimuli_count)
end

summary_names = ["seed", "count", "less_count"]
names = ["attribute", "hit_id", "user_id", "font_A_name", "font_B_name", "user_choice"]

for users_per_pair in users_per_pairs
    user_weights = Dict((i) => 1.0 for i in 1:users_per_pair)

    #summary = Tuple{Int, Int, Int}[]
    base_folder = "$(output_directory)/$(users_per_pair)"

    summary_filename = "$(base_folder)/summary.csv"
    if isfile(summary_filename)
        # if the summary already exists, load it
        summary = DataFrame(CSV.File(summary_filename))
    else 
        # if it doesn't, start a new one
        summary = DataFrame([name => Int[] for name in summary_names])
    end

    for seed in seeds
        if spacing == "expectation"
            # Load the appropriate file
            file = "data/stimuli/$(width)/seed_$(seed)_width_$(width)_stimuli.csv"
            #data = readdlm(file, ',', Float64)
            #global stimuli_distances = Dict(round(Int, data[i, 1]) => data[i, 2] for i in 1:size(data)[1])
            global stimuli_distances = CSV.File(file; header=false) |> Dict
        end

        if reliability !== nothing
            # Generally speaking, make sure you have enough users to conver...
            file = "$(reliability)/seed_$(seed)_users.csv"
            user_weights = CSV.File(file) |> Dict
        end

        rng = MersenneTwister(seed + 2)

        neighbors = shuffle(rng, 1:stimuli_count)
        neighbor_assignments = zeros(Int, stimuli_count - 1)
        for idx in 1:(stimuli_count - 1)
            n_0 = neighbors[idx]
            n_1 = neighbors[idx + 1]
            k = triLPair(max(n_0, n_1), min(n_0, n_1))
            # Add k to the assigned values 
            neighbor_assignments[idx] = k
        end
        # These are the remaining pairs, shuffled
        available_pairs = shuffle(rng, setdiff(1:pairs_per_attribute, neighbor_assignments))
        @assert(length(neighbor_assignments) == stimuli_count - 1)
        
        folder = "$(base_folder)/$(seed)/"

        mkpath(folder)
        for count in pair_counts
            assignment_results = Dict{Int, Array{Int}}()
            selected_pairs = available_pairs[1:(count - stimuli_count + 1)]

            # Assign according to the sorted array and remove from the shuffle
            all_pairs = vcat(neighbor_assignments, selected_pairs)
            @assert(length(all_pairs) == count)

            map(x -> assignment_results[x] = all_pairs, 1:users_per_pair)

            choice_rng = MersenneTwister(seed + 15)
            attribute = 1
            ch = Tuple{Int, Int, Int, Int, Int, String}[]
            user_keys = collect(keys(assignment_results))
            for user_id in user_keys
                assign = assignment_results[user_id]
                weight = user_weights[user_id]
                d = decisions(assign, stimuli_count, choice_rng; weight = weight, stimuli_distances = stimuli_distances)
                append!(ch, map(x -> (attribute, x[1], user_id, x[2], x[3], x[4]), d))
            end

            base_name = "$(folder)/seed_$(seed)_count_$(count)_per_pair_$(users_per_pair)"
            df = DataFrame(map(idx -> getindex.(ch, idx), eachindex(first(ch))), names)

            # Count the number of mismatches in the entire set
            less_occurances = sum(map(x -> x == "less", df.user_choice))
            if size(summary[.&(summary.seed .== seed, summary.count .== count), :])[1] > 0
                # update the dataframe
                summary[.&(summary.seed .== seed, summary.count .== count), :less_count] .= less_occurances
            else
                # Create a new row in the dataframe
                push!(summary, [seed, count, less_occurances])
            end

            # Write the file
            if format == "csv"
                CSV.write("$(base_name).csv", df; header=false)
            else
                df.winner = map(x -> x ? 1 : 2, df.user_choice .== "more")
                select!(df, Not("user_choice"))
        
                f = open("$(base_name).arrow.lz4", "w")
                Arrow.write(f, df; compress = :lz4)
                close(f)
            end
        end
    end
    CSV.write(summary_filename, summary)
end
