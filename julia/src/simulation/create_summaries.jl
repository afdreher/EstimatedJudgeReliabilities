# This code is used to generate the stimuli for the connected graph experiments


using ArgParse        # Argument parsing
using CSV
using DataFrames

include("../data_loading.jl")
include("filepaths.jl")
include("defaults.jl")

arg_settings = ArgParseSettings()
@add_arg_table arg_settings begin
    # Options
    "--seeds"
        help = "Seeds to process"
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

    "--spacing"
        help="Spacing mechanism (even or expectation)"
        arg_type = String
        default = "even"
        range_tester = (x -> x in ["even", "expectation"])

    "--series"
        help="Series to run (same or convenience)"
        arg_type = String
        default = "same"
        range_tester = (x -> x in ["same", "convenience"])
end

parsed_args = parse_args(ARGS, arg_settings)

width = parsed_args["width"]
spacing = parsed_args["spacing"]
series = parsed_args["series"]

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

output_directory = joinpath("data", series, spacing, string(width))
mkpath(output_directory)

summary_names = ["seed", "count", "less_count"]

for users_per_pair in users_per_pairs
    #summary = Tuple{Int, Int, Int}[]
    base_folder = joinpath(output_directory, string(users_per_pair))

    summary_filename = joinpath(base_folder, "summary.csv")
    if isfile(summary_filename)
        # if the summary already exists, load it
        summary_df = DataFrame(CSV.File(summary_filename))
    else 
        # if it doesn't, start a new one
        summary_df = DataFrame([name => Int[] for name in summary_names])
    end

    for seed in seeds
        for count in pair_counts

            base_file = data_file(series, spacing, width, users_per_pair, seed, count, "arrow.lz4")
            if isfile(base_file)
                df = load_preferences(base_file)
            else
                df = load_preferences(data_file(series, spacing, width, users_per_pair, seed, count, "csv"))
            end 

            # Count the number of mismatches in the entire set
            less_occurances = sum(map(x -> x == "less", df.user_choice))
            if size(summary_df[.&(summary_df.seed .== seed, summary_df.count .== count), :])[1] > 0
                # update the dataframe
                summary_df[.&(summary_df.seed .== seed, summary_df.count .== count), :less_count] .= less_occurances
            else
                # Create a new row in the dataframe
                push!(summary_df, [seed, count, less_occurances])
            end
        end
    end
    CSV.write(summary_filename, summary_df)
end
