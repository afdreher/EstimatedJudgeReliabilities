##
# These functions implement the core of the Negative Log calculation.
#
# I am breaking them out separately so I can test and reuse them.
##

# The sign for u really just determines whether you want larger or smaller values
# begin preferred.  It's symmetric, so... whatever.

sigmoid(u, w, l) = 1.0 / (1.0 + exp(-u * (w - l)))  # Core of the BTL model
btl(w, l) = sigmoid(1.0, w, l)                      # Unweighted
nll(p) = -log(p)                                    # Basic NLL definition
neg_log_P(u, w, l) = nll(sigmoid(u, w, l))          # Loss is NLL
neg_log_P(w, l) = neg_log_P(1.0, w, l)              # Unweighted

function nll_f(x...)
  # x[1] => reliablity of the user (possibly fixed)
  # x[2] => winner's value
  # x[3] => loser's value
  return sum(neg_log_P(x[i], x[i + 1], x[i + 2]) for i in 1:3:length(x))
end