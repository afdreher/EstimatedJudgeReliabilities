#!/usr/bin/env julia

using MathOptInterface
using JuMP            # Julia for Mathematical Optimization
using GLPK
using JSON
using Logging         # info, error, etc.
using Random

include("helpers.jl")
include("indexing.jl")


###
#  This script finds a feasible solution (if one exists) to the problem of
#  assigning pairs for stimuli
#
# NOTE: The outcomes will largely be random (to a point) because MOI does not
#       yet support setting a starting point.  Be forewarned.  Running the 
#       allocation script multiple times potentially yields different answers.
###

# Perform the allocation
function allocatePairs(stimuli_count::Int, pair_count::Int, rng)
    # Every stimulus is assigned an *approximately* equal number of pairs
    per_stimulus = (pair_count * 2) / stimuli_count
    min_per_stimulus = floor(Int, per_stimulus)
    max_per_stimulus = ceil(Int, per_stimulus)


    # --------------------------------------------------------------------------
    #  MODEL
    # --------------------------------------------------------------------------

    # Create a model
    model = Model(GLPK.Optimizer) # GLPK


    # --------------------------------------------------------------------
    #  VARIABLES
    # --------------------------------------------------------------------

    @variable(model, x[i=1:stimuli_count, j=1:stimuli_count], Bin)
    # GLPK throught MOI doesn't yet support this.  As a result, expect the outcomes
    # to be random, even when using the same seed.
    set_start_value.(x, rand(rng, Set([0.0, 1.0]), stimuli_count, stimuli_count))


    # --------------------------------------------------------------------
    #  OBJECTIVE: Spread costs to either randomize preferences, have users
    #             rate similar stimulli, or have users rate more distant
    #             stimuli
    # --------------------------------------------------------------------

    # Just add a bunch of random values for each person / pair value
    @objective(model, Min, sum(x .* rand(rng, Set([0.0, 1.0]), stimuli_count, stimuli_count)))


    # --------------------------------------------------------------------
    #  FEASIBILITY CONSTRAINTS
    # --------------------------------------------------------------------

    @constraints(model, begin
        # Constraint 1 - Every stimulus has an approximately equal number of assignments
        row_constraints[i in 1:stimuli_count], min_per_stimulus <= sum(x[i, :]) <= max_per_stimulus
        col_constraints[i in 1:stimuli_count], min_per_stimulus <= sum(x[:, i]) <= max_per_stimulus
        self_constraints[i in 1:stimuli_count], x[i, i] == 0

        # If the solver supports symmetry, you can use that instead.  GLPK does not, so enforce the symmetry here
        reflexive_constraints[i in 1:stimuli_count, j in 1:stimuli_count], x[i, j] == x[j, i]

        attribute_sum_constraint, sum(x) == 2*pair_count
    end)

    # Solve for a feasible solution
    JuMP.optimize!(model)

    term_status = JuMP.termination_status(model)
    primal_status = JuMP.primal_status(model)
    is_optimal = term_status == MOI.OPTIMAL

    if !is_optimal
        return Nothing
    end

    xs = JuMP.value.(x)

    # Verify that these are all independent...
    independent = true
    for i in 1:stimuli_count
        for j in (i+1):stimuli_count
            if sum(abs.(xs[i,:] - xs[j,:])) < 1
                #@warn("Rows $(i) and $(j) are not independent")
                independent = false
                break
            end
        end
    end

    if !independent
        @warn("Allocation is not independent")
        return nothing
    end
    # These are now the results where duplicates have been removed
    
    indices = findmat(x -> x > 0, xs);
    d = Tuple{Int64, Int64}[]
    for (u, v) in indices
        if u < v
            push!(d, (u, v))
        end
    end
    return d
end

# Perform the allocation
function allocatePairs(stimuli_count::Int,
                       pair_count::Int, 
                       attribute_count::Int, 
                       seed::Int)
    loop_number = 0
    number_solved = 0
    results = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
    for i in 1:(5 * attribute_count)
        rng = MersenneTwister(seed + loop_number)
        current = allocatePairs(stimuli_count, pair_count, rng)

        # Check the allocation...
        if current != nothing
            number_solved += 1
            # Check to make sure the current solution has not been seen
            @info("$(number_solved) appears OK")
            results[number_solved] = current
        end

        if number_solved == attribute_count
            break
        end

        loop_number += 1
    end
    return results
end