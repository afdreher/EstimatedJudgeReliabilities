# This code is used to generate user summary data for each file.

using ArgParse        # Argument parsing
using CSV
using DataFrames

include("../data_loading.jl")
include("filepaths.jl")
include("defaults.jl")


#SPACING_FOLDER_NAME = Dict("even" => "evenly", "expectation" => "expectation")

lesscount(a) = sum(map(x -> x == "less", a))


arg_settings = ArgParseSettings()
@add_arg_table arg_settings begin
    # Options
    "--seeds"
        help = "Seed to process"
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

for users_per_pair in [8] #users_per_pairs
    for seed in seeds
        for count in pair_counts
            base_folder = joinpath(output_directory, string(users_per_pair), string(seed), )
            summary_filename =  joinpath(base_folder,"seed_$(seed)_count_$(count)_user_summary.csv")

            base_file = data_file(series, spacing, width, users_per_pair, seed, count, "arrow.lz4")
            if isfile(base_file)
                df = load_preferences(base_file)
            else
                df = load_preferences(data_file(series, spacing, width, users_per_pair, seed, count, "csv"))
            end 

            summary_df = combine(DataFrames.groupby(df, :user_id), :user_choice => lesscount => :disagree, :user_id => length => :total)
            CSV.write(summary_filename, summary_df)
        end
    end
end
