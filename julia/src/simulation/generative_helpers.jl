# This file contains functions useful for generating the data

maxIndex(n::Int)::Int = (n * (n - 1)) / 2
checkBound(k::Int, n::Int)::Bool = k > 0 && k <= maxIndex(n)

function triUIndex(k::Int, n::Int)::Tuple{Int, Int}
    if !checkBound(k, n)
        @warn("Invalid index $(k).  Valid range is [1:$(maxIndex(n))].\n")
        return (0, 0)
    end
    
    i = n - 2 - floor(sqrt(-8 * (k - 1) + 4 * n * (n - 1) - 7) / 2.0 - 0.5)
    j = (k - 1) + i + 1 - n * (n - 1) / 2 + (n - i) * ((n - i) - 1) / 2
    return (i, j) .+ 1
end

# Convert from a lower triangular index (excluding the diagaonal) to a pair
#https://atrebas.github.io/post/2021-01-17-index_to_lower_triangular_subscripts/
function triLIndex(k::Int, n::Int)::Tuple{Int, Int}
    if !checkBound(k, n)
        @warn("Invalid index $(k).  Valid range is [1:$(maxIndex(n))].\n")
        return (0, 0)
    end
    kp = k - 1
    p  = floor((sqrt(1 + 8 * kp) - 1) / 2)
    i  = p + 2
    j  = kp - p * (p + 1) / 2 + 1
    return (i, j) 
end

# Convert from pair to lower triangular index (excluding the diagonal)
function triLPair(i::Int, j::Int)::Int
    if i <= j
        @warn("Invalid cell ($(i), $(j)). The row must be less than the column.\n")
        return 0
    end

    p  = i - 2
    kp = j - 1 + p * (p + 1) / 2
    k = kp + 1
    
    return k 
end

# I'm using O'Donovan et al.'s function to provide the *best* chance of this working.
odonovan(diff::Float64, weight::Float64 = 1.0)::Float64 = 1 / (1 + exp(- weight * diff))

choose(pair, p::Float64, rng) = ifelse(rand(rng) < p, pair[1], pair[2])

# # stimuli_distances should be Dict{Int, Float64}, but... it doesn't seem to work :/ 
# function decisions(assignments::Array{Int}, stimuli_count::Int, rng; weight = 1.0, stimuli_distances = nothing)
#     c = []
#     for k in assignments
#         pair = triLIndex(k, stimuli_count)
#         if stimuli_distances !== nothing
#             p1 = stimuli_distances[pair[1]]
#             p2 = stimuli_distances[pair[2]]
#         else            
#             p1 = pair[1]
#             p2 = pair[2]
#         end
#         choice = choose(["more", "less"], odonovan(float(p1 - p2), weight), rng)
#         push!(c, (k, pair[1], pair[2], choice))
#     end
#     return c
# end

unzip(a) = map(x->getfield.(a, x), fieldnames(eltype(a)))

function decisions(assignments::Array{Int}, stimuli_count::Int, rng; weight = 1.0, stimuli_distances = nothing)
    pairs = triLIndex.(assignments, stimuli_count)
    p1, p2 = unzip(pairs)
    differences = float(positions(p1, stimuli_distances) - positions(p2, stimuli_distances))
    choices = choose.(Ref(["more", "less"]), odonovan.(differences, weight), rng)
    return collect(zip(assignments, p1, p2, choices))
end

positions(pairs, stimuli_distances) = map(x -> stimuli_distances !== nothing ? stimuli_distances[x] : x, pairs)