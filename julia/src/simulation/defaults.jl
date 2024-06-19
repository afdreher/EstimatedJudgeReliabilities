# ------------------------------------------------------------------------------
#  DEFAULTS
# ------------------------------------------------------------------------------

DEFAULT_SEEDS = [
    103, 117, 236, 763, 375, 
    167, 192, 3209, 27, 947, 
    783, 365, 723,  12, 823, 
    713, 832, 512, 356, 707, 
    310, 916, 650, 209, 503, 
    480, 408, 714, 323, 805
]

# Default pair counts are 200-900 by 100, 1000-19500 by 500, and 19900 (all pairs)
DEFAULT_PAIR_COUNTS = vcat(collect(200:100:999), collect(1000:500:19900), 19900)
DEFAULT_USERS_PER_PAIRS = collect(8:8:32)
DEFAULT_WIDTHS = [2, 10, 50, 100]  # Divide this by 100 to match the paper

# Options are "connected" and "beta_10_2".  Connected is the one used in the
# paper, which grows a single connnected component.
DEFAULT_EXPERIMENTS = ["connected"]

# This is the sampling method
DEFAULT_SERIES = ["same", "convenience"]

# Spacings.  Options are "even" and "expectation".  "Even" is the one used in
# The paper
DEFAULT_SPACINGS = ["even"]

# All of the methods to test
ALL_METHODS = ["btl", "odonovan", "crowdbt"]


# ------------------------------------------------------------------------------
#  HELPERS
# ------------------------------------------------------------------------------

function or_default(args, key, default)
    if haskey(args, key)
        arr = parsed_args[key]
        if length(arr) < 1
            arr = default
        end
        return arr
    else
        return default
    end
end