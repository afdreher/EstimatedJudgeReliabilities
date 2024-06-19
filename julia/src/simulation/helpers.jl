using JSON

# ------------------------------------------------------------------------------
#  HELPER FUNCTIONS
# ------------------------------------------------------------------------------

function findmat(f, A::AbstractMatrix)
    m, n = size(A)
    out = Tuple{Int, Int}[]
    for i in 1:m, j in 1:n
      f(A[i, j]) && push!(out, (i, j))
    end
    out
  end
  
  function format_output(dict, format)
      if format == "tuple"
          # Convert the ouptut to tuples from linear indices
          td = Dict{Int64, Array{Tuple{Int64, Int64}}}()
          for (i, a) in enumerate(dict)
              td[i] = map(ind2pair, dict[i])
          end
          return JSON.json(td)
      else
          # Write the output using linear indices
          return  JSON.json(dict)
      end
  end