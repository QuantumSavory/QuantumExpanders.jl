using Nemo
using Nemo: zzModRingElem, MatSpace, zero_matrix, base_ring, lift
using Oscar
using LinearAlgebra
using Random

using QuantumClifford

const Z2Matrix = zzModMatrix

"""
    uniformly_random_code_checkmatrix(ρ, Δ)

Generate a random binary parity check matrix with rate ρ and block length Δ

# Example

```jldoctest example
julia> using QuantumExpanders: uniformly_random_code_checkmatrix # hide

julia> H = uniformly_random_code_checkmatrix(0.5, 10)
[0   0   1   1   0   0   1   1   0   1]
[1   1   0   0   1   0   0   0   1   0]
[1   1   1   0   0   0   0   1   1   1]
[1   1   1   0   0   1   0   1   1   1]
[0   0   0   1   1   1   0   0   1   0]
```
"""
function uniformly_random_code_checkmatrix(ρ, Δ)
    @assert 0<ρ<1
    r = Int(floor(ρ*Δ))
    uniformly_random_code_checkmatrix(r, Δ)
end

"""
    uniformly_random_code_checkmatrix(r, Δ)

Generate a random r×Δ binary parity check matrix of full rank.
"""
function uniformly_random_code_checkmatrix(r::Integer, Δ)
    R, _ = residue_ring(ZZ, 2)
    M = matrix_space(R, r, Δ)
    local H
    while true
        H = rand(M)
        rank(H) == r && break
    end
    H
end

"""
    dual_code(H)

Compute the dual code generator matrix for a given parity check matrix.

# Example

```jldoctest example
julia> using QuantumExpanders: dual_code # hide

julia> H = uniformly_random_code_checkmatrix(0.5, 10);

julia> G = dual_code(H)
[0   1   1   1   1   1   0   0   0   0]
[0   1   1   0   1   0   1   0   0   0]
[0   0   1   1   1   0   0   1   0   0]
[1   0   1   1   0   0   0   0   1   0]
[0   0   0   0   1   0   0   0   0   1]

julia> G * transpose(H) == zero_matrix(base_ring(G), size(G, 1), size(H, 1))
true
```
"""
function dual_code(H::Z2Matrix)
    null = nullspace(H)[2]
    @assert all(iszero, H * null)
    transpose(null)
end
function dual_code(H)
    r, Δ = size(H)
    R, _ = residue_ring(ZZ, 2)
    H = matrix_space(R, r, Δ)(R.(H))
    dual_code(H)
end

"""
    kron(l::Z2Matrix, r::Z2Matrix)

Compute the Kronecker product of two binary matrices over GF(2).

Both `l` and `r` must be matrices of type `Z2Matrix`, i.e., `Nemo.zzModMatrix`
over the ring GF(2).

# Example

```jldoctest
julia> using Nemo; R, _ = residue_ring(ZZ, 2);

julia> A = matrix(R, 2, 2, [1, 0, 0, 1])

julia> B = matrix(R, 1, 2, [1, 1]);

julia> kron(A, B)
[1   1   0   0]
[0   0   1   1]

"""
function LinearAlgebra.kron(l::Z2Matrix, r::Z2Matrix)
    l1, l2 = size(l)
    r1, r2 = size(r)
    s1 = l1 * r1
    s2 = l2 * r2
    M = matrix_space(base_ring(l), s1, s2)
    M(kron(Matrix(l), Matrix(r)))
end

"""Check that two binary parity check matrices X and Z result in a good CSS code
(i.e., commutation constraints are fulfilled)"""
function good_css(X::Union{Z2Matrix,AbstractMatrix{Bool}}, Z::Union{Z2Matrix,AbstractMatrix{Bool}})
    x = convert_to_bool(X)
    z = convert_to_bool(Z)
    _good_css_check(x, z)
end

function good_css(stab)
    for row in stab
        !all(==(0), QuantumClifford.comm(row, stab)) && return false
    end
    return true
end

convert_to_bool(M::Z2Matrix) = Bool.(Int.(lift.(M)))
convert_to_bool(M::AbstractAlgebra.MatElem{Bool}) = Matrix(M)
convert_to_bool(M::AbstractAlgebra.Generic.MatSpaceElem{Bool}) = Matrix{Bool}(M)
convert_to_bool(M::AbstractMatrix{<:Bool}) = Matrix{Bool}(M)

function _good_css_check(X, Z)
    stab = css(X, Z)
    good_css(stab)
end

"""Create a CSS code from two binary parity check matrices, X and Z"""
function css(X, Z)
    Xb = Matrix{Bool}(convert_to_bool(X))
    Zb = Matrix{Bool}(convert_to_bool(Z))
    vcat(Stabilizer(zeros(Bool, size(Xb)), Xb), Stabilizer(Zb, zeros(Bool, size(Zb))))
end
