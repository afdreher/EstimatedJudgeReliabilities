function ind2pair(idx, n)
    k = idx - 1  # Convert to 0-index
    i = Int(n - 2 - floor(sqrt(-8*k + 4*n*(n-1)-7)/2.0 - 0.5))
    j = Int(k + i + 1 - n*(n-1)/2 + (n-i)*((n-i)-1)/2)
    return (i + 1, j + 1)
end
