#!/usr/bin/env julia

using ArgParse            # Argument parsing
using CSV                 # CSV loading
using DataFrames          # Data frames (table storage)
using Distributions       # Compute inverse normal
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

using LinearAlgebra
using SparseArrays

include("data_loading.jl")

##
#  This implements the BTL model of Exploratory Font Selection Using Crowdsourced Attributes by
#  Peter O'Donovan, Janis Libeks, Aseem Agarwala, and Aaron Hertzmann
#
#  usage: julia odonovan.jl <preferences.csv>
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
    "--model"
        help = "Use either 2PL or 3PL. 2PL is the standard 2 parameter sigmoid model used in the paper.  3PL is the IRT version that includes a per-pair parameter."
        arg_type = String
        default = "2PL"
        range_tester = (x -> x in ["2PL", "3PL"])
    "--print-level"
        help = "IPOPT print level"
        arg_type = Int64
        default = 5
    "--verbose", "-v", "--v"
        help = "Run in verbose mode"
        action = :store_true
    "--estimate-user-reliability"
        help = "Estimate user reliability (Eq. 2 of section 4.3).  False (default) is the standard Bradley-Terry-Luce model."
        action = :store_true
    "--estimate-attribute-difficulty"
        help = "Estimate attribute difficulty.  This adds a additional variable"
        action = :store_true
    "--independent-users"
        help = "Use independent reliability scores for each user per attribute"
        action = :store_true
    "--fixed-scale-values"
        help = "File with fixed scale values.  Only the user weights will be estimated.  Values outside the scale will result in a error."
        arg_type = String
    "--fixed-user-weights"
        help = "File with fixed user weights.  Only the scale values will be estimated.  Values outside the scale will result in a error.  Requires estimate-user-reliability to be true for it to have an effect."
        arg_type = String
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
    "--user-min"
        help = "Minimum user reliability rating"
        arg_type = Float64
        default = -1.0
    "--user-max"
        help = "Maximum user reliability rating"
        arg_type = Float64
        default = 1.0
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
    "--use-flipped-indices"
        help = "Use v_i-v_j instead of 1-P(...)"
        action = :store_true
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
    "--loss-function"
        help = "Loss function to utilize.  Options are NLL (default), L1, L2, SQRT-L2, SQ-L2"
        arg_type = String
        default = "NLL"
        range_tester = (x -> x in ["NLL", "L1", "L2", "SQRT-L2", "SQ-L2"])
    "--seed"
        help = "Random seed"
        arg_type = Int64
        default = 117
    "--lambda"
        help = "Use centering fencepost.  If <= 0, no dummy vote will be used"
        arg_type = Float64
        default = -1.0
    "--gamma"
        help = "Ridge regression coefficient.  If <= 0, no ridge regression will be used"
        arg_type = Float64
        default = -1.0
    "--scale-gamma"
        help = "Scale gamma by the number of votes a participant provides"
        action = :store_true
    "--gamma-center"
        arg_type = Float64
        default = 1.0
    "--solver"
        help = "The solver to use"
        arg_type = String
        default = "IPOPT"
        range_tester = (x -> x in ["IPOPT", "SLSQP"])
    "--tol"
        help = "Solver tolerance"
        arg_type = Float64
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

  MODEL = get(config, "model", "2PL")

  PRINT_LEVEL = get(config, "print-level", 5)
  VERBOSE = get(config, "verbose", false)

  ESTIMATE_USER_RELIABILITY = get(config, "estimate-user-reliability", false)
  ESTIMATE_ATTRIBUTE_DIFFICULTY = get(config, "estimate-attribute-difficulty", false)
  INDEPENDENT_USERS = get(config, "independent-users", false)

  FIXED_SCALE_VALUES_FILE = get(config, "fixed-scale-values", nothing)
  FIXED_USER_WEIGHTS_FILE = get(config, "fixed-user-weights", nothing)

  MIN_ACCEPTABLE = get(config, "minimum-acceptable", 1)
  MIN_WEIGHT = get(config, "minimum-weight", 0.0)
  EDGE_TOLERANCE = get(config, "edge-tolerance", 1E-6)
  USER_MIN = get(config, "user-min", -1.0)
  USER_MAX = get(config, "user-max", 1.0)
  USER_START_VALUE = get(config, "user-start-value", nothing)
  SCALE_MIN = get(config, "scale-min", -100.0)
  SCALE_MAX = get(config, "scale-max", 100.0)
  FLIP = get(config, "use-flipped-indices", false)

  OMITTED_ATTRIBUTES = get(config, "omit-attributes", nothing)
  ONLY_ATTRIBUTES = get(config, "only-attributes", nothing)

  OMITTED_USERS = get(config, "omit-users", nothing)
  ONLY_USERS = get(config, "only-users", nothing)
  INVERTED_USERS = get(config, "invert-users", nothing)
  RANDOMIZED_USERS = get(config, "randomize-users", nothing)

  LOSS_FUNCTION =  get(config, "loss-function", "NLL")

  SEED = get(config, "seed", 117)
  LAMBDA = get(config, "lambda", -1.0)
  GAMMA = get(config, "gamma", -1.0)
  SCALE_GAMMA = get(config, "scale-gamma", false)
  GAMMA_CENTER = get(config, "gamma-center", 0.0)
  SOLVER = get(config, "solver", "IPOPT")
  TOL = get(config, "tol", nothing)

else
  EDGE_FILE = parsed_args["edges"]
  OUTPUT = parsed_args["output"]
  INITIAL_VALUES_FILE = parsed_args["initial-values"]

  MODEL = parsed_args["model"]

  PRINT_LEVEL = parsed_args["print-level"]
  VERBOSE = parsed_args["verbose"]
  ESTIMATE_USER_RELIABILITY = parsed_args["estimate-user-reliability"]
  ESTIMATE_ATTRIBUTE_DIFFICULTY = parsed_args["estimate-attribute-difficulty"]
  INDEPENDENT_USERS = parsed_args["independent-users"]
  FIXED_SCALE_VALUES_FILE = parsed_args["fixed-scale-values"]
  FIXED_USER_WEIGHTS_FILE = parsed_args["fixed-user-weights"]
  MIN_ACCEPTABLE = parsed_args["minimum-acceptable"]
  MIN_WEIGHT = parsed_args["minimum-weight"]
  EDGE_TOLERANCE = parsed_args["edge-tolerance"]
  USER_MIN = parsed_args["user-min"]
  USER_MAX = parsed_args["user-max"]
  USER_START_VALUE = parsed_args["user-start-value"]
  SCALE_MIN = parsed_args["scale-min"]
  SCALE_MAX = parsed_args["scale-max"]
  FLIP = parsed_args["use-flipped-indices"]

  OMITTED_ATTRIBUTES = parsed_args["omit-attributes"]
  ONLY_ATTRIBUTES = parsed_args["only-attributes"]

  OMITTED_USERS = parsed_args["omit-users"]
  ONLY_USERS = parsed_args["only-users"]
  INVERTED_USERS = parsed_args["invert-users"]
  RANDOMIZED_USERS = parsed_args["randomize-users"]

  LOSS_FUNCTION =  parsed_args["loss-function"]

  SEED = parsed_args["seed"]
  LAMBDA = parsed_args["lambda"]
  GAMMA = parsed_args["gamma"]
  SCALE_GAMMA = parsed_args["scale-gamma"]
  GAMMA_CENTER = parsed_args["gamma-center"]
  SOLVER = parsed_args["solver"]
  TOL = parsed_args["tol"]
end


if VERBOSE
  println("Parsed arguments:")
  @printf "           preferences: %s\n" PREFERENCE_FILE
  @printf "                 edges: %s\n" EDGE_FILE
  @printf "                output: %s\n" (OUTPUT === nothing ? "console" : OUTPUT)
  @printf "\n"
  @printf "                 model: %s\n" MODEL
  @printf "\n"
  @printf "               verbose: %s\n" (VERBOSE ? "true" : "false")
  @printf "      user-reliability: %s\n" (ESTIMATE_USER_RELIABILITY ? "true" : "false")
  @printf "  attribute-difficulty: %s\n" (ESTIMATE_ATTRIBUTE_DIFFICULTY ? "true" : "false")
  @printf "     independent-users: %s\n" (INDEPENDENT_USERS ? "true" : "false")
  if FIXED_SCALE_VALUES_FILE !== nothing
      @printf "    fixed-scale-values: %f\n" FIXED_SCALE_VALUES_FILE
  end
  if FIXED_USER_WEIGHTS_FILE !== nothing
      @printf "    fixed-user-weights: %f\n" FIXED_USER_WEIGHTS_FILE
  end
  @printf "    minimum-acceptable: %d\n" MIN_ACCEPTABLE
  @printf "        minimum-weight: %f\n" MIN_WEIGHT
  @printf "        edge-tolerance: %s\n" EDGE_TOLERANCE
  @printf "              user-min: %s\n" USER_MIN
  @printf "              user-max: %s\n" USER_MAX
  if USER_START_VALUE !== nothing
      @printf "      user-start-value: %f\n" USER_START_VALUE
  end
  @printf "             scale-min: %s\n" SCALE_MIN
  @printf "             scale-max: %s\n" SCALE_MAX
  @printf "                lambda: %s\n" LAMBDA
  @printf "                gamma: %s\n" GAMMA
  @printf "          scale-gamma: %s\n" (SCALE_GAMMA ? "true" : "false")
  @printf "         gamma-center: %s\n" GAMMA_CENTER
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
  @printf "\n"
  @printf "         loss-function: %s\n" LOSS_FUNCTION
end

function compute_optimal_hessian(model)
  d = NLPEvaluator(model)
  MOI.initialize(d, [:Hess])
  hessian_sparsity = MOI.hessian_lagrangian_structure(d)
  I = [i for (i, _) in hessian_sparsity]
  J = [j for (_, j) in hessian_sparsity]
  V = zeros(length(hessian_sparsity))
  x = all_variables(model)
  x_optimal = value.(x)
  y_optimal = dual.(all_nonlinear_constraints(model))
  MOI.eval_hessian_lagrangian(d, V, x_optimal, 1.0, y_optimal)
  n = num_variables(model)
  H = SparseArrays.sparse(I, J, V, n, n)
  vmap = Dict(x[i] => i for i in 1:n)
  add_to_hessian(H, f::Any, μ) = nothing
  function add_to_hessian(H, f::QuadExpr, μ)
      for (vars, coef) in f.terms
          if vars.a != vars.b
              H[vmap[vars.a], vmap[vars.b]] += μ * coef
          else
              H[vmap[vars.a], vmap[vars.b]] += 2 * μ * coef
          end
      end
  end
  for (F, S) in list_of_constraint_types(model)
      for cref in all_constraints(model, F, S)
          add_to_hessian(H, constraint_object(cref).func, dual(cref))
      end
  end
  add_to_hessian(H, objective_function(model), 1.0)
  return Matrix(fill_off_diagonal(H))
end


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

userParticipation = []

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

if PREFERENCE_FILE !== nothing
  # Try to read the preference file. This should be updated to take the winner / loser format, but that's not for now
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

  if MODEL == "3PL"
    # Hits only matter for the 3PL model, so avoid extra work & storage
    hits = unique(preferences[!, :hit_id])
    H = length(hits)
    hits2index = Dict(collect(zip(hits, 1:H)))
  end
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
  set_optimizer_attribute(m, "max_iter", 8000)
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

#
#  NOTE THAT YOU HAVE TO HAVE SCALE VALUES LIMITED TO USE ESTIMATE_USER_RELIABILITY
#

LAMBDA_OFFSET = USE_LAMBDA ? 1 : 0
@variable(m, SCALE_MIN <= v[i=1:(N + LAMBDA_OFFSET), 1:A] <= SCALE_MAX)   # Scale values


# Scale values are set randomly [0, 1] in the original O'Donovan / Libek code
if MODEL == "3PL"
  @variable(m, 0 <= c[i=1:H] <= 1, start=0)  # Item random guessing
end
if ESTIMATE_USER_RELIABILITY  # User rating
  if INDEPENDENT_USERS
    @variable(m, USER_MIN <= r[1:U, 1:A] <= USER_MAX, start=1.0)
  else
    @variable(m, USER_MIN <= r[1:U] <= USER_MAX, start=1.0)
  end
end
if ESTIMATE_ATTRIBUTE_DIFFICULTY
  @variable(m, 0 <= d[1:A] <= 1.0, start=1.0)
end

if ESTIMATE_USER_RELIABILITY 
  if INITIAL_VALUES_FILE !== nothing
      loaded_values = JSON.parsefile(INITIAL_VALUES_FILE)
      if haskey(loaded_values, "userNumber")
          for (key, value) in loaded_values["userNumber"]
              if haskey(users2index, key)
                  idx = users2index[key]
                  set_start_value(r[idx], float(value))
              end
          end
      end
  elseif USER_START_VALUE !== nothing
      if INDEPENDENT_USERS
        r_fixed = fill(USER_START_VALUE, U, A)
      else
        r_fixed = fill(USER_START_VALUE, U)
      end
      set_start_value.(r, r_fixed)
  else
      r_rng = MersenneTwister(SEED)

      if INDEPENDENT_USERS
        r_random = rand(r_rng, U, A)
      else
        r_random = rand(r_rng, U)
      end
      set_start_value.(r, r_random)
  end
end

# Set the starting stimuli positions to random
v_rng = MersenneTwister(SEED + 10)
v_random = rand(v_rng, N + LAMBDA_OFFSET, A)
if USE_LAMBDA
  v_random[N + LAMBDA_OFFSET, :] .= 0
end
set_start_value.(v, v_random)  # Does this work these days?


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

  if FLIP
    if choice == "more"
      #print("i (A) MORE THAN j (B)\n")
      diff_v = Expr(:call, :(-), v[j, a], v[i, a])
    else
      #print("j (B) LESS THAN i (A)\n")
      diff_v = Expr(:call, :(-), v[i, a], v[j, a])
    end
  else
    # If FLIP isn't set, we'll adjust this later
    diff_v = Expr(:call, :(-), v[j, a], v[i, a])
  end

  if ESTIMATE_USER_RELIABILITY
    if INDEPENDENT_USERS
      expon = Expr(:call, :(*), r[u, a], diff_v)
    else
      expon = Expr(:call, :(*), r[u], diff_v)
    end
  else
    expon = diff_v
  end

  if ESTIMATE_ATTRIBUTE_DIFFICULTY
    expon = Expr(:call, :(*), d[a], expon)
  end

  denom = Expr(:call, :(+), 1.0, Expr(:call, :(exp), expon))
  if MODEL == "3PL"
      h = hits2index[row[:hit_id]]
      p = Expr(:call, :(+), c[h], Expr(:call, :(/), Expr(:call, :(-), 1.0, c[h]), denom))
  else
      p = Expr(:call, :(/), 1.0, denom)
  end


  if !FLIP && choice == "less"
    p = Expr(:call, :(-), 1.0, p)
  end
  
  if LOSS_FUNCTION == "NLL"
      # USE NLL
      loss = Expr(:call, :(-), Expr(:call, :(log), p))
  else
    if LOSS_FUNCTION == "L1"
      loss = Expr(:call, :(abs), Expr(:call, :(-), 1, p))
    else
        if LOSS_FUNCTION == "SQRT-L2"
          intern = Expr(:call, :(^), p, 0.5)
        elseif LOSS_FUNCTION == "SQ-L2"
          intern = Expr(:call, :(^), p, 2)
        else
          # Vanilla L2
          intern = p
         end

        loss = Expr(:call, :(^), Expr(:call, :(-), 1, intern), 2)
    end
  end
  push!(terms, loss)  # Add the expression to an array
end

if ESTIMATE_USER_RELIABILITY  && GAMMA > 0
  for u in 1:U
    gamma = GAMMA
    if SCALE_GAMMA
      gamma = GAMMA * votes[u]
    end
    loss = Expr(:call, :(*), Expr(:call, :(^), Expr(:call, :(-), GAMMA_CENTER, r[u]), 2.0), gamma)
    push!(terms, loss)  # Add the expression to an array
  end
end

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
  edges = CSV.read(EDGE_FILE, DataFrame)

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
  @constraint(m, v[N+1, :] .== 0)
end

if FIXED_SCALE_VALUES_FILE !== nothing
  # Load the file and add the constraints
  scale_values = CSV.read(FIXED_SCALE_VALUES_FILE, DataFrame)
  println("Setting scale values...")
  # Weights file should read <stimulus> <attribute> <weight> 
  @error("Not yet implemented")

  for stimulus_row in eachrow(scale_values)
    # Grab the attribute and stimrNulus
    stimulus = stimulus_row[:stimulus]
    attribute = stimulus_row[:attribute]
    value = stimulus_row[:weight]

    if (value < SCALE_MIN || value > SCALE_MAX)
      @error("$(value) is out of range for ($(stimulus), $(attribute)). SKIPPED.")
      continue
    end

    a = attributes2index[attrib]
    s = stimuli2index[stimulus]

    @constraint(m, v[s, a] == value)
  end
end

if ESTIMATE_USER_RELIABILITY && FIXED_USER_WEIGHTS_FILE !== nothing
  # Load the file and add the constraints
  user_values = CSV.read(FIXED_USER_WEIGHTS_FILE, DataFrame)
  println("Setting user values...")
  # Weights file should read <userNumber> <weight> <attribute?>

  for user_row in eachrow(user_values)
    userNumber = user_row[:userNumber]
    value = user_row[:weight]
    
    if (value < USER_MIN || value > USER_MAX)
      @error("$(value) is out of range for user $(userNumber). SKIPPED.")
      continue
    end

    u = users2index[userNumber]
    if INDEPENDENT_USERS
      attribute = user_row[:attribute]
      a = attributes2index[attribute]
      @constraint(m, r[u, a] == value)
    else
      @constraint(m, r[u] == value)
    end
  end
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

arr = JuMP.value.(v)  # N x A matrix
ff = DataFrame(;OrderedDict(:name => stimuli, (Symbol(attributes[i])=>arr[1:N, i] for i=1:A)...)...)

if ESTIMATE_USER_RELIABILITY
  if INDEPENDENT_USERS
      u = Array{Union{Float64, Missing}}(missing, U, A)
      user_arr = JuMP.value.(r)  # U x A matrix
      # Copy over only values where participant has some decision
      for row in 1:U, col in 1:A
          if userParticipation[row, col]
              u[row, col] = user_arr[row, col]
          end
      end
      uf = DataFrame(;OrderedDict(:userNumber => users, (Symbol(attributes[i])=>u[:, i] for i=1:A)...)...)
  else
    # Write out all of the r's into a separate _cov file
    u = Array{Union{Float64, Missing}}(missing, U)  # use missing rather than nothing because dataframes complains otherwise
    user_arr = JuMP.value.(r)  # U x A matrix
    for row in 1:U
      if any(userParticipation[row, :])
          u[row] = user_arr[row]
      end
    end
    uf = DataFrame(userNumber = users, value = u)
  end
end

if MODEL == "3PL"
  hf = DataFrame(hitNumber = hits, value = JuMP.value.(c))
end

if ESTIMATE_ATTRIBUTE_DIFFICULTY
  ad = DataFrame(attribute = attributes, value = JuMP.value.(d))
end

# H_star = compute_optimal_hessian(m)
# eigs = LinearAlgebra.eigvals(H_star)
# print("H* eigs are: $(eigs)")
# print("nonzero eigs?: $(all(>=(0), eigs))")

if OUTPUT !== nothing
    open("$(OUTPUT[1:end-4])_status.txt", "w") do file
      write(file, "$(string(primal_status))\n")
    end

    # Write the dataframe
    CSV.write(OUTPUT, ff)

  if ESTIMATE_USER_RELIABILITY
    CSV.write("$(OUTPUT[1:end-4])_users.csv", uf)
  end

  if MODEL == "3PL"
    CSV.write("$(OUTPUT[1:end-4])_hits.csv", hf)
  end

  if ESTIMATE_ATTRIBUTE_DIFFICULTY
    CSV.write("$(OUTPUT[1:end-4])_attribute_difficulties.csv", ad)
  end
else
    # Print all of the variables, scale, etc.
    println("Final Solution:")
    println(ff)
    println("----------")
    if ESTIMATE_USER_RELIABILITY
      println(uf)
    end
    if MODEL == "3PL"
        println("   ----   ")
        println(hf)
    end
    if ESTIMATE_ATTRIBUTE_DIFFICULTY
      println("   ----   ")
      println(ad)
    end
    println()
end

# PRINT_HESSIAN = true
# if PRINT_HESSIAN
#     values = zeros(9)
#     @printf "Users start at: %d\n" linearindex(r[1])
#     for ui in 1:U
#         values[linearindex(r[users2index[uf[:name][ui]]])] = uf[:value][ui]
#     end
#
#     for ai in 1:A
#         attrib = Symbol(attributes[ai])
#         for ni in 1:N
#             n_idx = stimuli2index[ff[:name][ni]]
#             values[linearindex(v[ni, 1])] = ff[attrib][ni]
#         end
#     end
#
#
#     d = JuMP.NLPEvaluator(m)
#     MathProgBase.initialize(d, [:Grad, :Hess])
#     objval = MathProgBase.eval_f(d, values)
#
#     gradients = zeros(9)
#     MathProgBase.eval_grad_f(d, gradients, values)
#
#     HI, HJ = MathProgBase.hesslag_structure(d)
#     cnt = length(HI)
#     hessians = zeros(cnt)
#     o = ones(cnt)
#     MathProgBase.eval_hesslag_prod(d, hessians, values, o, 1.0 , o)
#
# end
