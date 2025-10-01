const Z2Matrix = zzModMatrix

"""
Generate a random binary parity check matrix for *classical component codes* used in
quantum Tanner code construction [leverrier2022quantum](@cite).

## Local Component Codes

The random classical code generation is essential for building quantum Tanner codes as
described in [leverrier2022quantum](@cite). *Theorem 18* of [leverrier2022quantum](@cite)
establishes that for sufficiently large ``\\Delta``, random codes with specific rates yield
*asymptotically good* quantum Tanner codes with high probability. Specifically, this theorem
requires choosing ``C_A`` via a random uniform ``r \\times \\Delta`` generator matrix and
``C_B`` via a random uniform ``r \\times \\Delta`` parity-check matrix, where``r = \\lfloor \\rho \\Delta \\rfloor``.

!!! note
    The parameters must satisfy ``-\\delta \\log_2 \\delta - (1-\\delta) \\log_2 (1-\\delta) < \\rho``
    (*Gilbert-Varshamov* bound). This condition ensures that random codes of rate ``\\rho`` can
    achieve relative minimum distance ``\\delta``, which is necessary for all four codes ``C_A``,
    ``C_B``, ``C_A^{\\perp}``, ``C_B^{\\perp}`` to have minimum distances ``\\geq \\delta\\Delta``
    with high probability. The parameter ``\\rho`` must satisfy ``0 < \\rho < 1/2`` to ensure
    non-trivial quantum codes.

The randomness is fundamental to obtaining codes with the robustness properties established in
*Theorem 9* of [leverrier2022quantum](@cite), which guarantees that with probability tending to 1
as ``\\Delta \\to \\infty``, the dual tensor code

```math
\\begin{aligned}
C_A \\otimes \\mathbb{F}_2^B + \\mathbb{F}_2^A \\otimes C_B
\\end{aligned}
```

is ``\\Delta^{3/2-\\varepsilon}``-robust with ``\\Delta^{\\gamma}``-resistance to puncturing.

# Example

```jldoctest example
julia> using QuantumExpanders: uniformly_random_code_checkmatrix;

julia> H = uniformly_random_code_checkmatrix(0.5, 10)
[0   0   1   1   0   0   1   1   0   1]
[1   1   0   0   1   0   0   0   1   0]
[1   1   1   0   0   0   0   1   1   1]
[1   1   1   0   0   1   0   1   1   1]
[0   0   0   1   1   1   0   0   1   0]
```

### Arguments
- `ρ::Real`: Target rate of the classical component code, determining the code's dimension relative to its block length.
- `Δ::Integer`: Block length of the classical component code, corresponding to the size of the generating sets in the underlying Cayley graph structure.
"""
function uniformly_random_code_checkmatrix(ρ, Δ)
    @assert 0<ρ<1
    r = Int(floor(ρ*Δ))
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
Compute the generator matrix of the *dual code* for classical component codes in quantum
Tanner construction [leverrier2022quantum](@cite). We work with pairs of classical codes and
their duals to construct quantum CSS codes. The dual code relationship is crucial because
quantum Tanner codes are defined through the pair ``(\\mathcal{C}_0, \\mathcal{C}_1)`` where

```math
\\begin{aligned}
\\mathcal{C}_0 = T(\\mathcal{G}_0^{\\square}, C_0^{\\perp}), \\quad C_0 = C_A \\otimes C_B
\\end{aligned}
```

and

```math
\\begin{aligned}
\\mathcal{C}_1 = T(\\mathcal{G}_1^{\\square}, C_1^{\\perp}), \\quad C_1 = C_A^{\\perp} \\otimes C_B^{\\perp}
\\end{aligned}
```

*Theorem 18* of [leverrier2022quantum](@cite) requires that both the component codes and
their duals have sufficiently large minimum distances, which is achieved with high probability
when the codes are randomly generated as specified.

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

### Arguments
- `H::zzModMatrix`: Parity check matrix of a classical linear code, typically one of the component codes ``C_A`` or ``C_B``.
"""
function dual_code(H)
    r, Δ = size(H)
    R, _ = residue_ring(ZZ, 2)
    H = matrix_space(R, r, Δ)(R.(H)) # TODO there must be a cleaner way to write this
    null = nullspace(H)[2]
    @assert all(iszero, H*null)
    return transpose(null)
end

"""
Compute the Kronecker product of binary matrices for constructing *tensor product codes*
in quantum Tanner codes of [leverrier2022quantum](@cite).

We form *tensor product codes* from the classical component codes. Specifically, we construct
``C_0 = C_A \\otimes C_B`` and ``C_1^{\\perp} = C_A^{\\perp} \\otimes C_B^{\\perp}``, which
are essential for defining the quantum code through Tanner codes on the left-right Cayley complex.

The *robustness* properties of these tensor codes, established in *Theorem 9* of
[leverrier2022quantum](@cite), are crucial for achieving the linear minimum distance bounds
in Theorem 18. This theorem proves that with high probability, the dual tensor code
``C_A \\otimes \\mathbb{F}_2^B + \\mathbb{F}_2^A \\otimes C_B`` is ``\\Delta^{3/2-\\varepsilon}``-robust
with ``\\Delta^{\\gamma}``-resistance to puncturing. Both `l` and `r` must be matrices of type
`Nemo.zzModMatrix` over the ring GF(2).

# Example

```jldoctest
julia> using Nemo; R, _ = residue_ring(ZZ, 2);

julia> A = matrix(R, 2, 2, [1, 0, 0, 1])

julia> B = matrix(R, 1, 2, [1, 1]);

julia> kron(A, B)
[1   1   0   0]
[0   0   1   1]

### Arguments

- `l::zzModMatrix`: Left matrix in the Kronecker product, typically a generator or parity check matrix of a component code.
- `r::zzModMatrix`: Right matrix in the Kronecker product, typically a generator or parity check matrix of the other component code.
"""
function LinearAlgebra.kron(l::Z2Matrix, r::Z2Matrix)
    l1, l2 = size(l)
    r1, r2 = size(r)
    s1 = l1*r1
    s2 = l2*r2
    M = matrix_space(base_ring(l), s1, s2)
    M(kron(Matrix(l), Matrix(r)))
end

"""Check that two binary parity check matrices X and Z result in a good CSS code
(i.e., commutation constraints are fulfilled)"""
function good_css(X::Union{Z2Matrix,AbstractMatrix{Bool}}, Z::Union{Z2Matrix,AbstractMatrix{Bool}})
    x = convert_to_bool(X)
    z = convert_to_bool(Z)
    stab = css(X, Z)
    for row in stab
        !all(==(0), QuantumClifford.comm(row, stab)) && return false
    end
    return true
end

convert_to_bool(M::Z2Matrix) = Bool.(Int.(lift.(M)))
convert_to_bool(M::AbstractAlgebra.MatElem{Bool}) = Matrix(M)
convert_to_bool(M::AbstractAlgebra.Generic.MatSpaceElem{Bool}) = Matrix{Bool}(M)
convert_to_bool(M::AbstractMatrix{<:Bool}) = Matrix{Bool}(M)

"""Create a CSS code from two binary parity check matrices, X and Z"""
function css(X, Z)
    Xb = Matrix{Bool}(convert_to_bool(X))
    Zb = Matrix{Bool}(convert_to_bool(Z))
    vcat(Stabilizer(zeros(Bool, size(Xb)), Xb), Stabilizer(Zb, zeros(Bool, size(Zb))))
end
