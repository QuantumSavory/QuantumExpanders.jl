"""
Puncturing and shortening are standard techniques for constructing new linear codes
from existing ones.

Let C be an [n,k,d] linear code over GF(q) with parity-check matrix H, and let t ⊆ {1,…,n}
be a set of coordinate positions specified by `cols`. 

The shortened code C_t is obtained by restricting C to codewords that are zero on t and
then puncturing those coordinates, resulting in a code of length n − |t|.

# Example

```jldoctest
julia> using QuantumExpanders; using Nemo

julia> H = [1 0 0 0 1 1;
            0 1 0 1 0 1;
            0 0 1 1 1 0];

julia> H_new = shorten(H, [6])
3×5 Matrix{Int64}:
 1  1  1  0  0
 0  1  0  1  0
 1  0  0  0  1

julia> rank(matrix(GF(2), H_new))
3
```
"""
function shorten(H::AbstractMatrix, cols::AbstractVector{<:Integer})
    G = Matrix{Int}(lift.(dual_code(matrix(ZZ, H))))
    keep_rows = [all(G[r,c] == 0 for c in cols) for r in 1:size(G,1)]
    G_new = G[keep_rows, :]
    keep_cols = setdiff(1:size(G_short,2), cols)
    G_new = G_new[:, keep_cols]
    H = Matrix{Int}(lift.(dual_code(matrix(ZZ, G_new))))
    return H
end

"""
Let C be an [n, k, d] linear code over GF(q) with parity-check matrix H_B, and let t be a
set of coordinates given by `cols`. The punctured code Cₜ is obtained by deleting the
coordinates in t from every codeword of C. The resulting code is linear and has length
n − |t|.

Returns a parity-check matrix for the punctured code Cₜ. The construction proceeds by
computing a generator matrix G of C from H_B, deleting the columns indexed by t, and then
computing a parity-check matrix for the resulting punctured code.

Here is an example of puncturing the classical [6,3,3] code:

```jldoctest
julia> using QuantumExpanders; using Nemo

julia> H = [1 0 0 0 1 1;
            0 1 0 1 0 1;
            0 0 1 1 1 0];

julia> H_new = puncture(H, [6])
2×5 Matrix{Int64}:
 1  1  1  0  0
 1  1  0  1  1

julia> rank(matrix(GF(2), H_new))
2
```

Now, it is a [5,2,2] code. This distance is verified from [dist-m4ri](https://github.com/QEC-pages/dist-m4ri) program.

"""
function puncture(H::AbstractMatrix, cols::AbstractVector{<:Integer})
    G = Matrix{Int}(lift.(dual_code(matrix(ZZ, H))))
    keep = setdiff(1:size(G,2), cols)
    G_p = G[:, keep]
    H = Matrix{Int}(lift.(dual_code(matrix(ZZ, G_p))))
    return H
end
