#!/usr/bin/env julia

using Logging         # info, error, etc.
using Random


###
#  This script attempts to find a feasible solution to assinging users to pairs.
#
#  The way this script works is to shuffle the pairs and then assign them to
#  users, largely randomly.  If a user already has a pair in that person's
#  rating set, then it will go to someone else.  In the event that a pair cannot
#  otherwise be given out (such as the last pair), then the system will randomly
#  select a donor to swap, provided that the donor can do so.
###

toDict(arr) = Dict((x) => arr[x] for x in eachindex(arr))

# Perform the assignments
assignPairs(pair_count::Int, user_counts::Array{Int}, users_per_pair::Int, rng) = assignPairs(pair_count, toDict(user_counts), users_per_pair, rng)

hasValue(d, key, target)::Bool = haskey(d, key) && target in d[key]

function assignPairs(pair_count::Int, user_counts::Dict{Int, Int}, users_per_pair::Int, rng)::Dict{Int, Array{Int}}
    # Shuffling appears wasteful, but it provides a *slightly* better simulation
    # of picking pairs at various times.
    shuffled_pairs = shuffle(rng, repeat(collect(1:pair_count), users_per_pair))

    d = Dict{Int, Set{Int}}()

    # This is not the absolute best way of doing this, but whaterver...
    available = copy(user_counts)
    for pair in shuffled_pairs
        users = collect(keys(available))

        assigned = 0        # The user who eventually got the pair
        selected_user = 0   # Current selection
        needs_untargeted_swap = false # Every eligible person full

        assign_candidates = users[map(k -> !hasValue(d, k, pair), users)]
        if length(assign_candidates) < 1
            #@info("No available candidates.  All users have $(pair).")
            needs_untargeted_swap = true
        end

        if !needs_untargeted_swap
            shuffle!(rng, assign_candidates)
            for selected_user in assign_candidates
                # If this user can accept the current pair
                if !haskey(d, selected_user)
                    assigned = selected_user
                    d[selected_user] = Set(pair)
                    break
                elseif !(pair in d[selected_user])
                    assigned = selected_user
                    push!(d[selected_user], pair)
                    break
                end
            end 
        end

        if assigned == 0
            #@info("Could not assign value: $(pair)")
            # Now, we should try to find a suitable swap
            
            # Try to randomly guess a swap
            candidates = collect(keys(d))
        
            if length(candidates) < 2
                @error("No available candidates. Cannot assign $(pair).")
                @error(d)
                continue
            end

            swap_candidates = candidates[map(k -> !hasValue(d, k, pair), candidates)]

            if length(swap_candidates) < 1
                @error("No available candidates. Cannot assign $(pair).")
                continue
            end

            shuffle!(rng, swap_candidates)

            if needs_untargeted_swap
                # Find someone with space available
                donors = shuffle(rng, users)
            else
                donors = [selected_user]
            end
            
            for selected_user in donors
                for selected_candidate in swap_candidates
                    if selected_user == selected_candidate
                        continue  # Can't select ourselves
                    end
                    
                    if pair in d[selected_candidate]
                        @warn("Pair should have already been validated")
                        continue # Not a valid target
                    end
                    
                    difference = collect(setdiff(d[selected_candidate], d[selected_user]))
                    # Select a value from the difference, if available
                    if length(difference) < 1
                        continue  # No difference
                    end
                    
                    candidate_length = length(d[selected_candidate])
                    target_length = length(d[selected_user])

                    # Valid!
                    swap_value = difference[rand(rng, 1:length(difference))]
                    delete!(d[selected_candidate], swap_value)
                    push!(d[selected_user], swap_value)
                    push!(d[selected_candidate], pair)
                    assigned = selected_user

                    if candidate_length != length(d[selected_candidate])
                        @error("Unsuccessful swap -- candidate length changed!")
                    end

                    if (target_length + 1) != length(d[selected_user])
                        @error("Unsuccessful swap!-- pair not added to user!")
                    end
                    
                    # Assigned != 0
                    break
                end

                # This will happen if the initial donor / selected user overlaps
                # completely with the others
                if (assigned != 0)
                    break
                end
            end

            if assigned == 0
                @error("Could not swap value $(pair).\nSelected: $(selected_user)\nCandidates: $(swap_candidates)\nd: $(d)")
            end
        end

        if assigned == 0
            @error("Could not swap value $(pair).")
        else
            available[assigned] -= 1
            if available[assigned] == 0
                delete!(available, assigned)
                deleteat!(users, findall(x->x==assigned, users))
            end
        end
    end

    return Dict((x) => sort(collect(d[x])) for x in keys(d) )
end

verifyAssignments(assignments::Dict{Int, Array{Int}}, pair_count::Int, user_counts::Array{Int}, users_per_pair::Int) = verifyAssignments(assignments, pair_count, toDict(user_counts), users_per_pair)

function verifyAssignments(assignments::Dict{Int, Array{Int}}, pair_count::Int, user_counts::Dict{Int, Int}, users_per_pair::Int)::Bool 
    # Verify assignments
    v = Dict{Int, Int}()
    is_valid_assignment = true
    for user in keys(assignments)
        items = collect(assignments[user])
        if length(items) != user_counts[user]
            @warn("User $(user) given an incorrect number of items. $(length(items)) != $(user_counts[user])")
            @warn(assignments)
            @warn(user_counts)
            is_valid_assignment = false
        end
        
        for item in items
            if !haskey(v, item)
                v[item] = 1
            else 
                v[item] += 1
            end
        end
    end

    for i in 1:pair_count
        if v[i] != users_per_pair
            @warn("Item $(i) has an incorrect number of pairs. $(v[i]) != $(users_per_pair)")
            is_valid_assignment = false
        end
    end

    return is_valid_assignment
end