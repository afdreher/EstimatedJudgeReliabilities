# ------------------------------------------------------------------------------
#  FILE PROCESSING
# ------------------------------------------------------------------------------

# This is the directory for the input data
data_directory(series::String, 
               spacing::String, 
               width::Int, 
               per_pair::Int, 
               seed::Int, 
               pair_count::Int)::String =
    joinpath("data", 
        series, 
        spacing, 
        string(width),
        string(per_pair), 
        string(seed))

# This is the resulting user file for O'Donovan or Chen (CrowdBT)
results_directory(series::String, 
                  spacing::String, 
                  width::Int, 
                  per_pair::Int, 
                  seed::Int, 
                  pair_count::Int, 
                  method::String)::String =
    joinpath("results",             
        series, 
        spacing, 
        method, 
        string(width),
        string(per_pair), 
        string(seed))

data_file(series::String, 
          spacing::String, 
          width::Int, 
          per_pair::Int, 
          seed::Int, 
          pair_count::Int, 
          extension::String)::String =
    joinpath(
        data_directory(series, spacing, width, per_pair, seed, pair_count),
        "seed_$(seed)_count_$(pair_count)_per_pair_$(per_pair).$(extension)")
            
file_base(series::String, 
          spacing::String, 
          width::Int, 
          per_pair::Int, 
          seed::Int, 
          pair_count::Int, 
          method::String, 
          appendix::String)::String =
    joinpath(
        results_directory(series, spacing, width, per_pair, seed, pair_count, method),
        "seed_$(seed)_count_$(pair_count)_per_pair_$(per_pair)$(appendix)")

# This is the resulting user file for O'Donovan or Chen (CrowdBT)
scale_file(series::String, 
           spacing::String, 
           width::Int, 
           per_pair::Int, 
           seed::Int, 
           pair_count::Int, 
           method::String)::String =
    file_base(series, spacing, width, per_pair, seed, pair_count, method, ".csv")

# This is the resulting user file for O'Donovan or Chen (CrowdBT)
users_file(series::String, 
           spacing::String, 
           width::Int, 
           per_pair::Int, 
           seed::Int, 
           pair_count::Int, 
           method::String)::String =
    file_base(series, spacing, width, per_pair, seed, pair_count, method, "_users.csv")

user_summary_file(series::String, 
                  spacing::String, 
                  width::Int, 
                  per_pair::Int, 
                  seed::Int, 
                  pair_count::Int, 
                  method::String)::String = 
    file_base(series, spacing, width, per_pair, seed, pair_count, method, "_user_summary.csv")