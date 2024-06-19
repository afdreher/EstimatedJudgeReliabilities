#!/usr/bin/env julia

using ArgParse            # Argument parsing
using CSV                 # CSV loading
using DataFrames          # Data frames (table storage)
using Distributions       # Compute inverse normal
using IndexedTables
using Ipopt               # Coin-OR interior point solver
using JSON                # JSON support
using JuMP                # Julia for Mathematical Optimization
using Logging             # info, error, etc.
using MathOptInterface    # MOI for inspection of results
using NLopt               # Non-linear package for JuMP
using OrderedCollections  # OrderedDict
using Printf              # Printing
using Random              # Random number generator
using YAML

include("data_loading.jl")

##
#  This implements the core of the Crowd-BT model (minus the regularizing term) 
#  of Pairwise Ranking Aggregation in a Crowdsourced Setting by Xi Chen and 
#  Paul N. Bennett.
#
#  usage: julia crowdbt.jl <preferences.csv>
##

# --------------------------------------------------------------------------------------------------
#  ARGUMENT / SETUP
# --------------------------------------------------------------------------------------------------

arg_settings = ArgParseSettings()
@add_arg_table arg_settings begin
    "preferences"
        help = "Raw preference file (CSV format)"
        required = true
    "--experiment"
        help = "Use experiment definition file (YAML).  Options provided on the command line will be ignored."
        arg_type = String
    "--edges"
        help = "Edges to enforce (CSV format)"
        arg_type = String
    "--initial-values"
        help = "Initial font/attribute values (JSON)"
        arg_type = String
    "--output", "-o", "--o"
        help = "Output the results to CSV files"
        arg_type = String
    "--print-level"
        help = "IPOPT print level"
        arg_type = Int64
        default = 5
    "--verbose", "-v", "--v"
        help = "Run in verbose mode"
        action = :store_true
    "--dependent-users"
        help = "Use one reliability score for user"
        action = :store_true
    "--minimum-acceptable", "--ma"
        help = "Minimum number of acceptable votes to count for preferences."
        arg_type = Int64
        default = 1
    "--minimum-weight", "--mw"
        help = "Minimum edge weight to be added to edge enforcement constraints."
        arg_type = Float64
        default = 0.0
    "--edge-tolerance"
        help = "Minimum edge distance"
        arg_type = Float64
        default = 1E-6
    "--user-start-value"
        help = "A fixed starting value for all users"
        arg_type = Float64
    "--scale-min"
        help = "Minimum scale value"
        arg_type = Float64
        default = -100.0
    "--scale-max"
        help = "Maximum scale value"
        arg_type = Float64
        default = 100.0
    "--omit-attributes"
        help = "Attributes to skip.  Use either omit or only, but not both."
        action = :store_arg
        nargs = '+'
    "--only-attributes"
        help = "Only run with these attributes; skip others.  Use either omit or only, but not both."
        action = :store_arg
        nargs = '+'
    "--omit-users"
        help = "Users to skip.  Use either omit or only, but not both."
        action = :store_arg
        nargs = '+'
    "--only-users"
        help = "Only run with these users; skip others.  Use either omit or only, but not both."
        action = :store_arg
        nargs = '+'
    "--invert-users"
        help = "Invert the preferences from these users to see if it matters"
        action = :store_arg
        nargs = '+'
    "--randomize-users"
        help = "Randomize the preferences from these users to see if it matters"
        action = :store_arg
        nargs = '+'
    "--seed"
        help = "Random seed"
        arg_type = Int64
        default = 117
    "--lambda"
        help = "Use centering fencepost.  If <= 0, no dummy vote will be used"
        arg_type = Float64
        default = -1.0
    "--solver"
        help = "The solver to use"
        arg_type = String
        default = "IPOPT"
        range_tester = (x -> x in ["IPOPT", "SLSQP"])
    "--tol"
        help = "Solver tolerance"
        arg_type = Float64
    "--starting-scale-file"
        help = "Starting scale values (optional)"
        arg_type = String

end

parsed_args = parse_args(ARGS, arg_settings)

PREFERENCE_FILE = parsed_args["preferences"]
EXPERIMENT_FILE = parsed_args["experiment"]

if EXPERIMENT_FILE !== nothing
  # Load the definitions from the experiment file
  config = YAML.load_file(EXPERIMENT_FILE; dicttype=Dict{String,Any})
 
  OUTPUT = get(config, "output", nothing)

  EDGE_FILE = get(config, "edges", nothing)
  INITIAL_VALUES_FILE = get(config, "initial-values", nothing)

  PRINT_LEVEL = get(config, "print-level", 5)
  VERBOSE = get(config, "verbose", false)

  DEPENDENT_USERS = get(config, "dependent-users", false)

  MIN_ACCEPTABLE = get(config, "minimum-acceptable", 1)
  MIN_WEIGHT = get(config, "minimum-weight", 0.0)
  EDGE_TOLERANCE = get(config, "edge-tolerance", 1E-6)
  USER_MIN = get(config, "user-min", -1.0)
  USER_MAX = get(config, "user-max", 1.0)
  USER_START_VALUE = get(config, "user-start-value", nothing)
  SCALE_MIN = get(config, "scale-min", -100.0)
  SCALE_MAX = get(config, "scale-max", 100.0)

  OMITTED_ATTRIBUTES = get(config, "omit-attributes", nothing)
  ONLY_ATTRIBUTES = get(config, "only-attributes", nothing)

  OMITTED_USERS = get(config, "omit-users", nothing)
  ONLY_USERS = get(config, "only-users", nothing)
  INVERTED_USERS = get(config, "invert-users", nothing)
  RANDOMIZED_USERS = get(config, "randomize-users", nothing)

  SEED = get(config, "seed", 117)
  LAMBDA = get(config, "lambda", -1.0)
  SOLVER = get(config, "solver", "IPOPT")
  TOL = get(config, "tol", nothing)

  SCALE_FILE = get(config, "starting-scale-file", nothing)
else
  EDGE_FILE = parsed_args["edges"]
  OUTPUT = parsed_args["output"]
  INITIAL_VALUES_FILE = parsed_args["initial-values"]

  PRINT_LEVEL = parsed_args["print-level"]
  VERBOSE = parsed_args["verbose"]
  DEPENDENT_USERS = parsed_args["dependent-users"]
  MIN_ACCEPTABLE = parsed_args["minimum-acceptable"]
  MIN_WEIGHT = parsed_args["minimum-weight"]
  EDGE_TOLERANCE = parsed_args["edge-tolerance"]
  USER_START_VALUE = parsed_args["user-start-value"]
  SCALE_MIN = parsed_args["scale-min"]
  SCALE_MAX = parsed_args["scale-max"]

  OMITTED_ATTRIBUTES = parsed_args["omit-attributes"]
  ONLY_ATTRIBUTES = parsed_args["only-attributes"]

  OMITTED_USERS = parsed_args["omit-users"]
  ONLY_USERS = parsed_args["only-users"]
  INVERTED_USERS = parsed_args["invert-users"]
  RANDOMIZED_USERS = parsed_args["randomize-users"]

  SEED = parsed_args["seed"]
  LAMBDA = parsed_args["lambda"]
  SOLVER = parsed_args["solver"]
  TOL = parsed_args["tol"]

  SCALE_FILE = parsed_args["starting-scale-file"]
end


if VERBOSE
  println("Parsed arguments:")
  @printf "           preferences: %s\n" PREFERENCE_FILE
  @printf "                 edges: %s\n" EDGE_FILE
  @printf "                output: %s\n" (OUTPUT === nothing ? "console" : OUTPUT)
  @printf "\n"
  @printf "               verbose: %s\n" (VERBOSE ? "true" : "false")
  @printf "       dependent-users: %s\n" (DEPENDENT_USERS ? "true" : "false")
  @printf "    minimum-acceptable: %d\n" MIN_ACCEPTABLE
  @printf "        minimum-weight: %f\n" MIN_WEIGHT
  @printf "        edge-tolerance: %s\n" EDGE_TOLERANCE
  if USER_START_VALUE !== nothing
      @printf "      user-start-value: %f\n" USER_START_VALUE
  end
  @printf "             scale-min: %s\n" SCALE_MIN
  @printf "             scale-max: %s\n" SCALE_MAX
  @printf "                lambda: %s\n" LAMBDA
  @printf "                sovler: %s\n" SOLVER
  if TOL !== nothing
    @printf "             tolerance: %f\n" TOL
  end
  @printf "\n"
  @printf "       omit-attributes: %s\n" OMITTED_ATTRIBUTES
  @printf "       only-attributes: %s\n" ONLY_ATTRIBUTES
  @printf "\n"
  @printf "            omit-users: %s\n" OMITTED_USERS
  @printf "            only-users: %s\n" ONLY_USERS
  @printf "          invert-users: %s\n" INVERTED_USERS
  @printf "       randomize-users: %s\n" RANDOMIZED_USERS
  @printf "   starting-scale-file: %s\n" (SCALE_FILE === nothing ? "none" : SCALE_FILE)
end

intersection(a, b) = findall(in(b),a)

# --------------------------------------------------------------------------------------------------
# DATA LOADING
# --------------------------------------------------------------------------------------------------

if VERBOSE
  println("Loading data...")
end

# Create associative array structure
# Yes, this should be a class, but whatever...
attributes = []        # Convert from idx -> attribute_name
attributes2index = []  # Convert from attribute_name -> idx
stimuli = []           # Convert from idx -> font_name
stimuli2index = []     # Convert from font_name -> idx
users = []             # Convert from idx -> user_id
users2index = []       # Convert from user_id -> idx
hits = []              # Convert from idx -> hit_id
hits2index = []        # Convert from hit_id -> idx

edges = []
A = 0
N = 0
U = 0
H = 0


USE_OMITTED_ATTRIBUTES = (OMITTED_ATTRIBUTES !== nothing && length(OMITTED_ATTRIBUTES) > 0)
USE_ONLY_ATTRIBUTES = (USE_OMITTED_ATTRIBUTES == false && ONLY_ATTRIBUTES !== nothing && length(ONLY_ATTRIBUTES) > 0)

USE_OMITTED_USERS = (OMITTED_USERS !== nothing)
USE_ONLY_USERS = (USE_OMITTED_USERS == false && ONLY_USERS !== nothing)
USE_LAMBDA = LAMBDA > 0

userParticipation = []


if PREFERENCE_FILE !== nothing
  preferences = load_preferences(PREFERENCE_FILE)

  # Filter the attributes, if desired
  attributes = unique(preferences[!, :attribute])
  if USE_OMITTED_ATTRIBUTES
    if VERBOSE
      println("Using omitted attributes")
    end
    attributes = filter!(e -> !(e in OMITTED_ATTRIBUTES), attributes)
  elseif USE_ONLY_ATTRIBUTES
    if VERBOSE
      println("Using only attributes")
    end
    attributes = filter!(e -> (e in ONLY_ATTRIBUTES), attributes)
  end
  if VERBOSE
    print("Using attributes", attributes, "\n")
  end

  # Get the unique font and user names / counts
  stimuli = unique(vcat(unique(preferences[!, :font_A_name]), unique(preferences[!, :font_B_name])))
  users = unique(preferences[!, :user_id])
  if USE_OMITTED_USERS
    users = filter!(e -> !(e in OMITTED_USERS), users)
  elseif USE_ONLY_USERS
    users = filter!(e -> (e in ONLY_USERS), users)
  end

  A = length(attributes)  # Number of attributes
  N = length(stimuli)     # Number of stimuli
  U = length(users)       # Number of users

  attributes2index = Dict(collect(zip(attributes, 1:A)))
  stimuli2index = Dict(collect(zip(stimuli, 1:N)))
  users2index = Dict(collect(zip(users, 1:U)))

  userParticipation = falses(U, A)
else
  print("ERROR: No preference file")
  exit()
end


# --------------------------------------------------------------------------------------------------
#  MODEL
# --------------------------------------------------------------------------------------------------

if VERBOSE
  println("Creating optimization model...")
end

# Don't use GN_DIRECT_L, GN_CRS2_LM, LN_SBPLX (Seg fault)
# GD_STOGO, GD_STOGO_RAND, GN_ESCH all seem too slow

# IPOPT works fine
# LD_MMA, LD_SLSQP, GN_ISRES, LN_COBYLA, LD_CCSAQ are okay


if SOLVER == "SLSQP"
  # Careful! This has not been tested since the update!
  m = Model(NLopt.Optimizer)
  set_optimizer_attribute(m, "algorithm", :LD_SLSQP)
else
  m = Model(Ipopt.Optimizer)
  set_optimizer_attribute(m, "max_iter", 50000)
  set_optimizer_attribute(m, "print_level", PRINT_LEVEL)
  if TOL !== nothing
     set_optimizer_attributes(m, "tol" => TOL)
  end
end

###
#  LD_SLSQP is very fast and generally finds very good minima
#


# --------------------------------------------------------------------------------------------------
#  VARIABLES
# --------------------------------------------------------------------------------------------------

LAMBDA_OFFSET = USE_LAMBDA ? 1 : 0
@variable(m, SCALE_MIN <= v[i=1:(N + LAMBDA_OFFSET), 1:A] <= SCALE_MAX)   # Scale values

if DEPENDENT_USERS
    @variable(m, 0 <= r[1:U] <= 1, start=1.0)
else
    @variable(m, 0 <= r[1:U, 1:A] <= 1, start=1.0)
end

if INITIAL_VALUES_FILE !== nothing
    loaded_values = JSON.parsefile(INITIAL_VALUES_FILE)
    if haskey(loaded_values, "workers")
        for (key, value) in loaded_values["workers"]
            if haskey(users2index, key)
                idx = users2index[key]
                set_start_value(r[idx], float(value))
            end
        end
    end
elseif USER_START_VALUE !== nothing
    if DEPENDENT_USERS
        r_fixed = fill(USER_START_VALUE, U)
    else
        r_fixed = fill(USER_START_VALUE, U, A)
    end
    set_start_value.(r, r_fixed)
else
    r_rng = MersenneTwister(SEED)
    if DEPENDENT_USERS
        r_random = rand(r_rng, U)
    else
        r_random = rand(r_rng, U, A)
    end
    set_start_value.(r, r_random)  # Does this work these days?
end

# Set the starting stimuli positions to random
v_rng = MersenneTwister(SEED + 10)
v_start = rand(v_rng, N + LAMBDA_OFFSET, A)

if SCALE_FILE !== nothing
  starting_scales = CSV.read(SCALE_FILE, DataFrame)

  # We want the intersection of attributes available in the file and to the
  # optimizer. Everything else is 0.
  available_names = names(starting_scales)
  valid_attributes = map(x -> available_names[x], intersection(available_names, attributes))
  
  available_stimuli = starting_scales["name"]
  valid_stimuli = map(x -> available_stimuli[x], intersection(available_stimuli, stimuli))

  for stim in valid_stimuli
    n = stimuli2index[stim]
    my_row = filter(row -> row["name"] == stim, starting_scales)
    for attrib in valid_attributes
      a = attributes2index[attrib]
      v_start[n, a] = my_row[attrib][1] 
      #print("Setting ", stim, " (", n, ") and attribute ", attrib, "(", a, ") to ", my_row[attrib][1], "\n")
    end
  end
end
if USE_LAMBDA
  v_start[N + LAMBDA_OFFSET, :] .= 0
end
set_start_value.(v, v_start)  # Does this work these days?


# --------------------------------------------------------------------------------------------------
#  OBJECTIVE: Minimize the difference between the real x_jk and the projected x_jk
# --------------------------------------------------------------------------------------------------

# NLopt doesn't support vector indexing, so we have to do some kind of metaprogramming. Yes, it's
# not quite ideal, but oh well.

terms = []  # This will collect all of the terms

HAS_INVERTED_USERS = INVERTED_USERS !== nothing
HAS_RANDOMIZED_USERS = RANDOMIZED_USERS !== nothing

usr_rng = MersenneTwister(SEED + 5)
votes = Dict{Int, Float64}()
for row in eachrow(preferences)
  # This allows for filtered attributes
  attrib = row[:attribute]
  if !(attrib in attributes)
    continue
  end

  user_id = row[:user_id]
  if !(user_id in users)
    continue
  end

  a = attributes2index[attrib]
  i = stimuli2index[row[:font_A_name]]
  j = stimuli2index[row[:font_B_name]]
  u = users2index[user_id]
  choice = row[:user_choice]

  userParticipation[u, a] = true

  if haskey(votes, u)
    votes[u] += 1.0
  else
    votes[u] = 1.0
  end

  if HAS_INVERTED_USERS && user_id in INVERTED_USERS
    if choice == "more"
      choice = "less"
    else
      choice = "more"
    end
  elseif HAS_RANDOMIZED_USERS && user_id in RANDOMIZED_USERS
    choice = rand(usr_rng, ["more", "less"])
  end

    if choice == "more"
      #print("i (A) MORE THAN j (B)\n")
      forward_v = Expr(:call, :(-), v[j, a], v[i, a])   # (j - i)
      reverse_v = Expr(:call, :(-), v[i, a], v[j, a])   # (i - j)
    else
      #print("j (B) LESS THAN i (A)\n")
      forward_v = Expr(:call, :(-), v[i, a], v[j, a])   # (i - j)
      reverse_v = Expr(:call, :(-), v[j, a], v[i, a])   # (j - i)
    end

    # 1 + exp(-x)
    forward_denom = Expr(:call, :(+), 1.0, Expr(:call, :(exp), forward_v))
    reverse_denom = Expr(:call, :(+), 1.0, Expr(:call, :(exp), reverse_v))


    # u / (1 + exp(-x)) + (1-u)/(1 + exp(x))
    if DEPENDENT_USERS
        forward_p = Expr(:call, :(/), r[u], forward_denom)
        reverse_p = Expr(:call, :(/), Expr(:call, :(-), 1.0, r[u]), reverse_denom)
    else
        forward_p = Expr(:call, :(/), r[u, a], forward_denom)
        reverse_p = Expr(:call, :(/), Expr(:call, :(-), 1.0, r[u, a]), reverse_denom)
    end
    p =  Expr(:call, :(+), forward_p, reverse_p)

    # USE NLL  (Max sum log in original is min -log)
    loss = Expr(:call, :(-), Expr(:call, :(log), p))
    push!(terms, loss)  # Add the expression to an array
end

# Lambda R
if USE_LAMBDA
  # Lambda R
  for a in 1:A, n in 1:N
      # One forward / one reverse
      vs = [Expr(:call, :(-), v[n, a], v[N + 1, a]), Expr(:call, :(-), v[N + 1, a], v[n, a])]
      for v in vs
          denom = Expr(:call, :(+), 1.0, Expr(:call, :(exp), v))
          p = Expr(:call, :(/), 1.0, denom) 
          log_p = Expr(:call, :(-), Expr(:call, :(log), p)) 
          loss = Expr(:call, :(*), LAMBDA, log_p) # Mutliply by lambda 
          push!(terms, loss)
      end
  end
end

exp = Expr(:call, :(+))   # Expression is sum over all components
append!(exp.args, terms)  # Add the array to the argument list

# exp now encodes the sum over all of the valid pairs.  We want to minimize that value
set_nonlinear_objective(m, MIN_SENSE, exp)


# --------------------------------------------------------------------------------------------------
#  CONSTRAINTS
# --------------------------------------------------------------------------------------------------

# @constraint(m, v[1, 1] == 1) # Set a stake in the ground
# for u in 1:U
#     @constraint(m, r[u] == 3)
# end

# Add in the edge constraints
# for each edge
if EDGE_FILE !== nothing
  # Try to load the file
  # Process the file into tuples (attribute, winner_idx, loser_idx, weight, count)
  edges = CSV.read(EDGE_FILE)

  println("Adding edges complete...")
  for i in 1:length(edges[:winner])
    # This allows for filtered attributes
    attrib = edges[:attribute][i]
    if !(attrib in attributes)
      continue
    end

    attr_idx = attributes2index[attrib]
    win_idx = stimuli2index[edges[:winner][i]]  # Get the index of the winner
    los_idx = stimuli2index[edges[:loser][i]]   # Get the index of the loser
    score = edges[:score][i]
    count = edges[:count][i]

    if count >= MIN_ACCEPTABLE && score >= MIN_WEIGHT
        # Winner has at least as much of the attribute as the loser (ASC)
        @constraint(m, v[los_idx, attr_idx] - v[win_idx, attr_idx] + EDGE_TOLERANCE <= 0)
    end
  end
end

if USE_LAMBDA
  # Set the last entry to be the fixed point (0)
  @constraint(m, v[N+1, :] .== 0)
end


if VERBOSE
  println("Model complete...")
  println(m)
end


# --------------------------------------------------------------------------------------------------
#  SOLVE
# --------------------------------------------------------------------------------------------------

if VERBOSE
  println("Solving model...")
end

# Solve for a feasible solution
JuMP.optimize!(m)

term_status = JuMP.termination_status(m)
primal_status = JuMP.primal_status(m)
is_optimal = term_status == MOI.OPTIMAL


println("Objective value: ", JuMP.objective_value(m))

arr = JuMP.value.(v)  # (N+1) x A matrix
ff = DataFrame(;OrderedDict(:name => stimuli, (Symbol(attributes[i])=>arr[1:N, i] for i=1:A)...)...)

# Convert the raw FF into final scales [0, 100] ea.

if  DEPENDENT_USERS
    # Write out all of the r's into a separate _cov file
    u = Array{Union{Float64, Missing}}(missing, U)  # use missing rather than nothing because dataframes complains otherwise
    user_arr = JuMP.value.(r)  # U x A matrix
    # Copy over only values where participant has some decision.  It should be
    # all, but it might not be because of filtering
    for row in 1:U
        if any(userParticipation[row, :])
            u[row] = user_arr[row]
        end
    end

    uf = DataFrame(userNumber = users, value = u)
else
    u = Array{Union{Float64, Missing}}(missing, U, A)
    user_arr = JuMP.value.(r)  # U x A matrix
    # Copy over only values where participant has some decision
    for row in 1:U, col in 1:A
        if userParticipation[row, col]
            u[row, col] = user_arr[row, col]
        end
    end
    uf = DataFrame(;OrderedDict(:userNumber => users, (Symbol(attributes[i])=>u[:, i] for i=1:A)...)...)
end

if OUTPUT !== nothing
    open("$(OUTPUT[1:end-4])_status.txt", "w") do file
      write(file, "$(string(primal_status))\n")
    end

    # Write the dataframe
    CSV.write(OUTPUT, ff)

    CSV.write("$(OUTPUT[1:end-4])_users.csv", uf)
else
    # Print all of the variables, scale, etc.
    println("Final Solution:")
    println(ff)
    println("----------")
    println(uf)
    println()
end