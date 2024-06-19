using ArgParse            # Argument parsing
using CSV
using Distributions
using DataFrames
using Random

include("../data_loading.jl")
include("randomized_assign_pairs.jl")
include("defaults.jl")


arg_settings = ArgParseSettings()
@add_arg_table arg_settings begin
    # Options
    "--seeds"
        help = "Seeds to run"
        action = :store_arg
        nargs = '*'
        arg_type = Int64

    "--pair-counts"
        help = "Pair counts to run"
        action = :store_arg
        nargs = '*'
        arg_type = Int64

    "--users-per-pair"
        help = "Users per pair to run"
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

    "--spacing"
        help="Spacing mechanism (even or expectation)"
        arg_type = String
        default = "even"
        range_tester = (x -> x in ["even", "expectation"])

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

@info("Extension: $(extension)")

# @info("SEEDS:    $(parsed_args["seeds"])")
# @info("PER-PAIR: $(parsed_args["users-per-pair"])")
# @info("COUNTS:   $(parsed_args["counts"])")

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
    pair_counts = DEFAULT_COUNTS
end

print("Seeds:    $(seeds)\n")
print("Per-pair: $(users_per_pairs)\n")
print("Counts:   $(pair_counts)\n")
print("Spacing:  $(spacing)\n")

# It's probably more correct to use an exponential becuase that more naturally models the Poisson 
# process that's likely a reasonable assumption of how people choose to continue picking assignments 
# on MTurk.  However, Gamma is not a bad or unreasonable choice either.  Exponential is a special 
# case of Gamma.  This fit seems to have the heavier tail that we want.
user_pair_distribution = Gamma(0.3655)          # Solved using Scipy
user_pair_distribution_scale = 834.8795         # Multiplier to expand scale

for users_per_pair in users_per_pairs
    @info("Processing users $(users_per_pair)")
    for seed in seeds
        @info("\tProcessing seed $(seed)")

        input_user_directory = "data/same/$(spacing)/$(width)/$(users_per_pair)/$(seed)"
        output_directory = "data/convenience/$(spacing)/$(width)/$(users_per_pair)/$(seed)"
        
        mkpath(output_directory)
        for count in pair_counts
            @info("\t\tProcessing count $(count)")
            user_pair_distribution_offset = min(4.0, count) # Minimum is 4 in the data, so let's set that
            user_pair_minimum = min(4.0, count)
            user_pair_maximum = count
    
            rng = MersenneTwister(seed)
            
            base_file = "$(input_user_directory)/seed_$(seed)_count_$(count)_per_pair_$(users_per_pair)"
            @info("\tBase file $(base_file)")
            if isfile("$(base_file).arrow.lz4")
                df = load_preferences("$(base_file).arrow.lz4")
                rename!(df,:hit_id => :pair_id)
            else
                df = load_preferences("$(base_file).csv")
                # df = DataFrame(CSV.File(".csv", header=false))
                # rename!(df, ["attribute", "pair_id", "user_id", "item_a", "item_b", "choice"])
            end 
            unique_pairs = combine(groupby(df, [:pair_id]), nrow => :count)

            # Verify...
            @assert(all(unique_pairs.count .== users_per_pair) && size(unique_pairs, 1) == count)

            target = 1 * count * users_per_pair

            accepted = Int[]

            total = 0
            while total < target
                # Grab a new random value
                value = min(round(Int, rand(rng, user_pair_distribution)
                                * user_pair_distribution_scale
                                + user_pair_distribution_offset), user_pair_maximum)
                total += value
                push!(accepted, value)
                if total >= target
                    break
                end
            end
            accepted_sum = sum(accepted)
            @assert(total == accepted_sum, "Initial: $(total) ≠ $(accepted_sum)")
            overshoot = total - target

            removed = 0
            user_counts = map(idx -> accepted[idx], 1:length(accepted))
            if (length(user_counts) * user_pair_minimum) > target      
                @info("Must remove entire participant")
                # If any match exactly, remove that participant
                min_idx = 1
                min_value = 1E6
                for idx in 1:length(accepted)
                    if user_counts[idx] < min_value
                        min_value = user_counts[idx]
                        min_idx = idx
                    end

                    if user_counts[idx] == overshoot
                        removed = overshoot
                        overshoot = 0
                        deleteat!(accepted, idx)
                        break
                    end
                end

                if removed == 0
                    # Otherwise, remove the minimum
                    removed = user_counts[min_idx]
                    overshoot -= removed
                    if overshoot < 0
                        for rep in 1:(-10*overshoot)
                            # Randomly select one for removal
                            idx = rand(rng, Set(1:length(accepted)))
                            if idx == min_idx || accepted[idx] == user_pair_maximum
                                continue
                            end

                            accepted[idx] += 1
                            overshoot += 1
                        
                            if overshoot == 0
                                break
                            end
                        end
                    end

                    deleteat!(accepted, min_idx)
                end
            end

            for rep in 1:(10*overshoot)
                # Randomly select one for removal
                idx = rand(rng, Set(1:length(accepted)))
                if accepted[idx] > user_pair_minimum
                    accepted[idx] -= 1
                    removed += 1
                end
            
                if removed >= overshoot
                    break
                end
            end
            filtered_sum = sum(accepted)
            @assert(target == filtered_sum, "Filtered: $(target) ≠ $(filtered_sum)")
            #@assert(target == sum(accepted))

            pairs_assignments = assignPairs(count, accepted, users_per_pair, rng)
            if (!verifyAssignments(pairs_assignments, count, accepted, users_per_pair))
                @warn("seed_$(seed)_count_$(count)_per_pair_$(users_per_pair) failed verification. Skipping.")
                continue
            end

            users_with_pair = Dict{Int, Set{Int}}()

            for v in keys(pairs_assignments)
                for l in pairs_assignments[v]
                    k = unique_pairs.pair_id[l]
                    if !haskey(users_with_pair, k)
                        users_with_pair[k] = Set(v)
                    else
                        push!(users_with_pair[k], v)
                    end
                end
            end

            new_users = Int[]
            for pid in df.pair_id
                user = rand(rng, users_with_pair[pid])
                delete!(users_with_pair[pid], user)
                push!(new_users, user)
                #print("Assigning $(user) to $(pid). $(length(users_with_pair[pid])) remain.\n")
            end

            df.user_id = new_users
            # Write the file
            base_name = "$(output_directory)/seed_$(seed)_count_$(count)_per_pair_$(users_per_pair)"
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
end
