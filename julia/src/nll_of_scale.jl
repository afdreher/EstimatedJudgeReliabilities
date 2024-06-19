# This file computes the NLL of a scale optionally with user weights.  Note that
# this does not perform exclusion.  It's not a held out set test; it's just 
# goodness of fit.
#
# To process a scale use the command:
#
#    julia src/nll_of_scale.jl [path to preferences] [path to scale file] --users [path to user weights file]

using ArgParse        # Argument parsing
using Arrow           # Optional
using CSV             # CSV loading
using DataFrames      # Data frames (table storage)
using Formatting      # Python-like print formatting
using Printf          # printf / @printf

# include("./stats.jl")  # Compute Bayes factors against random or alternatives

include("negative_log.jl")
include("data_loading.jl")


arg_settings = ArgParseSettings()
@add_arg_table! arg_settings begin
    # PROBLEM STATEMENT
    "preferences"
        help = "Initial preferences (CSV or arrow.lz4)"
        required = true
        arg_type = String

    "scales"
        help = "Scales solution (CSV)"
        required = true
        arg_type = String

    "--users"
        help = "User weights (CSV)"
        arg_type = String
end

parsed_args = parse_args(ARGS, arg_settings)

# PREFERENCE_FILE = "/Users/afdreher/Google Drive/Code/odonovan_preferences.csv"
# USER_FILE = "/Users/afdreher/Google Drive/Code/odonovan_score_results.csv"
# SCALE_FILE = "/Users/afdreher/Google Drive/Code/odonovan_weighted_results_raw.csv"

PREFERENCE_FILE = parsed_args["preferences"]
SCALE_FILE = parsed_args["scales"]
USER_FILE = parsed_args["users"]

if PREFERENCE_FILE !== nothing
    # Try to read the preference file
    #preferences = CSV.read(PREFERENCE_FILE, DataFrame)
    preferences = load_preferences(PREFERENCE_FILE, true)
else
    print("ERROR: No preference file")
    exit()
end

scales = Dict{String, Dict{String, Float64}}()
if SCALE_FILE !== nothing
    # Try to read the preference file
    raw_data = CSV.read(SCALE_FILE, DataFrame)
    attributes = names(raw_data)[2:end]
    for attribute in attributes
        scales[string(attribute)] = Dict{String, Float64}()
    end

    for row in eachrow(raw_data)
        name = string(row[:name])
        for attribute in attributes
            scales[string(attribute)][name] = row[attribute]
        end
    end
else
    print("ERROR: No scale file")
    exit()
end

if USER_FILE !== nothing
    HAS_USER_WEIGHTS = true
    raw_data = CSV.read(USER_FILE, DataFrame)
    user_weights = Dict{Int, Float64}()
    for row in eachrow(raw_data)
        user_weights[row[:userNumber]] = row[:value]
    end
else
    HAS_USER_WEIGHTS = false
end

nll_of_attributes = Dict{String, Float64}()
for row in eachrow(preferences)
    attrib = string(row[:attribute])
    if !(attrib in keys(nll_of_attributes))
        nll_of_attributes[attrib] = 0.0
    end

    winner = string(row[:winner])
    loser = string(row[:loser])
    user_number = row[:user_id]

    winner_value = scales[attrib][winner]
    loser_value = scales[attrib][loser]

    user_weight = 1.0
    if HAS_USER_WEIGHTS
        user_weight = user_weights[user_number]
    end

    #println("$(winner) $(loser) => $(user_weight), $(winner_value), $(loser_value)")
    nll_of_attributes[attrib] += neg_log_P(user_weight, winner_value, loser_value)
end

for attribute in sort(attributes)
    @printf("%-20s: %0.3f\n", attribute, nll_of_attributes[attribute])
end
print("------------------------------\n")
#@printf("%-20s: %0.3f\n", "Total", sum(values(nll_of_attributes)))
@printf("%-20s: %0.13f\n", "Total", sum(values(nll_of_attributes)))