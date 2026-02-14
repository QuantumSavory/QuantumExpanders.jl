"""

"The puncturing and shortening technique are two important approaches to constructing new linear
codes from old ones. Many families of linear codes with interesting parameters have been obtained
with the puncturing technique [liu2020shortenedlinearcodesfinite](@cite)". 

Let C be an [n, k, d] linear code over GF(q) with parity-check matrix H, and let t be a set
of coordinates given by `cols`. Let C(t) denote the subcode of C consisting of codewords
that are zero on t. Puncturing C(t) on t produces a linear code of length n − |t| called
the shortened code of C, denoted Cₜ.

Returns a parity-check matrix of the shortened code Cₜ, obtained by deleting from H the columns indexed by t.

```jldoctest
julia> using QuantumExpanders; using Nemo

julia> H = [1 0 0 0 1 1;
            0 1 0 1 0 1;
            0 0 1 1 1 0];

julia> H_new = shorten(H, [6])
3×5 Matrix{Int64}:
 1  0  0  0  1
 0  1  0  1  0
 0  0  1  1  1

julia> rank(matrix(GF(2), H_new))
3
 ```

Now, it is a [5,3,3] code. This distance is verified from dist-m4ri program.
"""
function shorten(H::AbstractMatrix, cols::AbstractVector{<:Integer})
    keep = setdiff(1:size(H,2), cols)
    return H[:, keep]
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

Now, it is a [5,2,2] code. This distance is verified from dist-m4ri program.

"""
function puncture(H::AbstractMatrix, cols::AbstractVector{<:Integer})
    G = Matrix{Int}(lift.(dual_code(matrix(ZZ, H))))
    keep = setdiff(1:size(G,2), cols)
    G_p = G[:, keep]
    H = Matrix{Int}(lift.(dual_code(matrix(ZZ, G_p))))
    return H
end
