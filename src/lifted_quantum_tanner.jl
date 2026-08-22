"""
Left multiplication permutation matrix ``L_a`` on ``\\mathbb{F}_2^{|G|}``.

Given an element ``a \\in G`` and an indexing ``\\mathrm{idx} \\colon G \\to \\{1,\\dots,|G|\\}``
of the group elements, ``L_a`` is the permutation matrix that sends the basis vector ``e_g`` to ``e_{ag}``.
Together with [`_rho`](@ref), these are the building blocks of the lifted quantum Tanner construction of
[leverrier2025small](@cite): commuting left/right translations on the group algebra ``\\mathbb{F}_2[G]``.

Returned as a sparse matrix so that Kronecker products with the tensor-code parts stay sparse; only ``|G|``
nonzero entries per generator.
"""
function _lambda(a, idx)
    n = length(idx)
    rows = Int[]
    cols = Int[]
    for (g, i) in idx
        push!(rows, i)
        push!(cols, idx[a*g])
    end
    sparse(rows, cols, ones(Int, n), n, n)
end

"""
Right multiplication permutation matrix ``R_b`` on ``\\mathbb{F}_2^{|G|}``.

Given an element ``b \\in G``, ``R_b`` sends the basis vector ``e_g`` to ``e_{g b^{-1}}``. The
inverse is used so that ``R_b`` and [`_lambda`](@ref) commute: ``L_a R_b = R_b L_a`` for all ``a, b \\in G``.
    This commutation is what makes the lifted construction produce a valid CSS code.

See [leverrier2025small](@cite) equation (2) and the surrounding discussion.
"""
function _rho(b, idx)
    n = length(idx); b_inv = inv(b)
    rows = Int[]
    cols = Int[]
    for (g, i) in idx
        push!(rows, i)
        push!(cols, idx[g*b_inv])
    end
    sparse(rows, cols, ones(Int, n), n, n)
end

"""
Assemble the block-diagonal left/right translation operators ``L_A`` and ``R_B`` used by [`parity_matrices`](@ref).

Concretely, for generating multisets
``A = (a_1, \\dots, a_{n_A})`` and ``B = (b_1, \\dots, b_{n_B})``:

```math
L_A = \\bigoplus_{i,j} L_{a_i}, \\qquad
R_B = \\bigoplus_{i,j} R_{b_j},
```

each block-diagonal with ``n_A \\cdot n_B`` blocks of size ``|G| \\times |G|``. The ordering interleaves
``A`` and ``B`` so that a downstream Kronecker ``\\text{something} \\otimes I_{|G|}`` acts blockwise in the right way.

### Arguments
- `group`: finite group ``G`` (any Oscar group, iterated via `collect`)
- `A::Vector`: multiset of ``n_A`` generators (may repeat elements)
- `B::Vector`: multiset of ``n_B`` generators (may repeat elements)

### Returns
Tuple `(LA, RB)` of sparse ``|G| n_A n_B \\times |G| n_A n_B`` permutation matrices.
"""
function build_LA_RB(group, A, B)
    elems = collect(group)
    n_G   = length(elems)
    idx   = Dict(g => k for (k,g) in enumerate(elems))
    blocks_LA = [_lambda(A[i], idx) for i in 1:length(A) for j in 1:length(B)]
    blocks_RB = [_rho(B[j],   idx) for i in 1:length(A) for j in 1:length(B)]
    blockdiag(blocks_LA...), blockdiag(blocks_RB...)
end

"""
Build the ``X`` and ``Z`` parity-check matrices of the lifted quantum Tanner code from
equation (2) of [leverrier2025small](@cite).

Let ``\\mathcal{C}_i = \\ker H_i = \\mathrm{Im}(G_i^\\top)`` and ``\\mathcal{C}_i' = \\ker H_i' = \\mathrm{Im}(G_i'^\\top)``
be the two pairs of dual classical codes on ``\\mathbb{F}_2^{n_A}`` and ``\\mathbb{F}_2^{n_B}`` respectively.
The lifted quantum Tanner code has parity-check matrices

```math
H_X = \\begin{pmatrix}
H_0 \\otimes G_0' \\otimes I_{|G|} \\\\
\\bigl(H_1 \\otimes G_1' \\otimes I_{|G|}\\bigr) L_A R_B
\\end{pmatrix},
\\qquad
H_Z = \\begin{pmatrix}
\\bigl(G_0 \\otimes H_1' \\otimes I_{|G|}\\bigr) R_B \\\\
\\bigl(G_1 \\otimes H_0' \\otimes I_{|G|}\\bigr) L_A
\\end{pmatrix},
```

where ``L_A, R_B`` are as in [`build_LA_RB`](@ref). Commutation of ``L_A`` and ``R_B``,
combined with ``H_i G_i^\\top = 0``, gives ``H_X H_Z^\\top = 0`` — a valid CSS code.

The total block length is ``n = n_A \\, n_B \\, |G|``.

### Arguments
- `LA`, `RB`: block-diagonal translation operators from [`build_LA_RB`](@ref)
- `n_G::Int`: group order ``|G|``
- `H0, H1`: two parity-check matrices on the ``A``-side (``H_0, H_1 \\in \\mathbb{F}_2^{r \\times n_A}``)
- `G0, G1`: matching generator matrices satisfying ``H_i G_i^\\top = 0``
- `H0p, H1p, G0p, G1p`: analogous ``B``-side matrices with column count ``n_B``

### Returns
Tuple `(HX, HZ)` of dense ``\\{0,1\\}`` matrices modulo 2.

# See also
[`QuantumTannerViaLeftRightActions`](@ref) is the high-level wrapper that verifies
CSS orthogonality and assembles a full code from a group and generating sets.
"""
function parity_matrices(LA, RB, n_G, H0, H1, G0, G1, H0p, H1p, G0p, G1p)
    I_G = sparse(1:n_G, 1:n_G, ones(Int, n_G))
    HX  = mod.(vcat(
        Matrix(kron(kron(sparse(H0),  sparse(G0p)), I_G)),
        Matrix(kron(kron(sparse(H1),  sparse(G1p)), I_G)*LA*RB)), 2)
    HZ  = mod.(vcat(
        Matrix(kron(kron(sparse(G0),  sparse(H1p)), I_G)*RB),
        Matrix(kron(kron(sparse(G1),  sparse(H0p)), I_G)*LA)), 2)
    HX, HZ
end

"""
CSS orthogonality check helper: verifies ``H G^\\top = 0`` over ``\\mathbb{F}_2`` and that
``G`` spans the full kernel of ``H`` (i.e. that ``H`` and ``G`` are exact duals).
"""
function _check_duality(H, G, label)
    F2 = GF(2)
    Hm = matrix(F2, Matrix{Int}(H))
    Gm = matrix(F2, Matrix{Int}(G))
    iszero(Hm * transpose(Gm)) ||
        throw(ArgumentError("$label: H·Gᵀ ≠ 0 over F_2"))
    rank(Gm) == size(H,2) - rank(Hm) ||
        throw(ArgumentError("$label: rows of G do not span ker H " *
            "(rank G = $(rank(Gm)), n − rank H = $(size(H,2)-rank(Hm)))"))
end

"""
Quantum Tanner code from the *lifted* left–right action construction of [leverrier2025small](@cite).

This is an alternative to [`QuantumTannerCode`](@ref). Instead of the square-complex
description of Leverrier & Zémor [leverrier2022quantum](@cite), the code is presented
as two classical Tanner codes lifted through commuting left and right actions of a
finite group ``G`` on the group algebra ``\\mathbb{F}_2[G]``. The two presentations
give the same family of codes but differ in what is convenient to search over: this
construction removes the total non-conjugacy condition on ``(A, B)`` and lets the
Cayley graphs on the ``A``-slice and ``B``-slice be inspected separately, which is
useful for search heuristics.

Given a finite group ``G``, two multisets ``A, B \\subseteq G`` (may contain repeats),
and two pairs of dual classical codes on ``\\mathbb{F}_2^{n_A}`` and ``\\mathbb{F}_2^{n_B}``,
the code has

- ``n = n_A \\, n_B \\, |G|`` physical qubits
- parity-check matrices ``H_X, H_Z`` given by [`parity_matrices`](@ref)
- CSS orthogonality ``H_X H_Z^\\top = 0`` verified at construction

# Constructors

Three constructors are exported, in increasing generality:

```julia
QuantumTannerViaLeftRightActions(group, A, B, H_A, H_B; p1=…, p2=…)
```
The **simplest form**: pass a group, two generating multisets, and two parity-check
matrices ``H_A, H_B``. Generator matrices are computed via [`dual_code`](@ref) and the
column permutations `p1, p2` reorder the second pair ``(H_1, G_1)`` and ``(H_1', G_1')``
respectively; both default to the identity permutation.

```julia
QuantumTannerViaLeftRightActions(group, A, B, H_A, G_A, H_B, G_B; p1=…, p2=…)
```
Same as above but with generator matrices supplied explicitly, saving a `dual_code`
call. Column permutations `p1, p2` again reorder the second pair.

```julia
QuantumTannerViaLeftRightActions(group, A, B,
    H0_A, H1_A, G0_A, G1_A, H0_B, H1_B, G0_B, G1_B)
```
**Fully explicit form**: pass all four ``A``-side matrices and all four ``B``-side
matrices directly. Useful when the two parity checks on a side are not related by a
column permutation.

# Examples

# See also
- [`QuantumTannerCode`](@ref) — the square-complex construction, useful for exhaustive
  search over ``(A, B)`` satisfying the total non-conjugacy condition.
- [`parity_matrices`](@ref) and [`build_LA_RB`](@ref) — the underlying primitives.
- [`dual_code`](@ref) — used to compute generator matrices from parity checks.

# Fields
    $TYPEDFIELDS
"""
struct QuantumTannerViaLeftRightActions{GT, ET} <: AbstractCSSCode
    """The underlying finite group ``G``."""
    group::GT
    """First generating multiset ``A \\subseteq G`` (may contain repeats)."""
    A::Vector{ET}
    """Second generating multiset ``B \\subseteq G`` (may contain repeats)."""
    B::Vector{ET}
    """``A``-side parity check matrix ``H_0``."""
    H0_A::Matrix{Int}
    """``A``-side parity check matrix ``H_1`` (typically a column permutation of ``H_0``)."""
    H1_A::Matrix{Int}
    """``A``-side generator matrix ``G_0``, satisfying ``H_0 G_0^\\top = 0``."""
    G0_A::Matrix{Int}
    """``A``-side generator matrix ``G_1``, satisfying ``H_1 G_1^\\top = 0``."""
    G1_A::Matrix{Int}
    """``B``-side parity check matrix ``H_0'``."""
    H0_B::Matrix{Int}
    """``B``-side parity check matrix ``H_1'``."""
    H1_B::Matrix{Int}
    """``B``-side generator matrix ``G_0'``, satisfying ``H_0' G_0'^\\top = 0``."""
    G0_B::Matrix{Int}
    """``B``-side generator matrix ``G_1'``, satisfying ``H_1' G_1'^\\top = 0``."""
    G1_B::Matrix{Int}
    """``X``-check matrix ``H_X`` of the quantum code."""
    hx::Matrix{Int}
    """``Z``-check matrix ``H_Z`` of the quantum code."""
    hz::Matrix{Int}

    function QuantumTannerViaLeftRightActions(group, A::Vector, B::Vector,
            H0_A, H1_A, G0_A, G1_A, H0_B, H1_B, G0_B, G1_B)
        nA = size(H0_A, 2)
        nB = size(H0_B, 2)
        length(A) == nA ||
            throw(ArgumentError("|A| = $(length(A)) ≠ n_A = $nA"))
        length(B) == nB ||
            throw(ArgumentError("|B| = $(length(B)) ≠ n_B = $nB"))
        _check_duality(H0_A, G0_A, "A-side (H0,G0)")
        _check_duality(H1_A, G1_A, "A-side (H1,G1)")
        _check_duality(H0_B, G0_B, "B-side (H0',G0')")
        _check_duality(H1_B, G1_B, "B-side (H1',G1')")
        n_G    = length(collect(group))
        LA, RB = build_LA_RB(group, A, B)
        hx, hz = parity_matrices(LA, RB, n_G,
            H0_A, H1_A, G0_A, G1_A, H0_B, H1_B, G0_B, G1_B)
        iszero(mod.(hx * hz', 2)) ||
            throw(ArgumentError("internal error: HX·HZᵀ ≠ 0 (CSS orthogonality violated)"))
        new{typeof(group), eltype(A)}(group, A, B,
            Matrix{Int}(H0_A), Matrix{Int}(H1_A),
            Matrix{Int}(G0_A), Matrix{Int}(G1_A),
            Matrix{Int}(H0_B), Matrix{Int}(H1_B),
            Matrix{Int}(G0_B), Matrix{Int}(G1_B), hx, hz)
    end
end

# Convenience constructor: user supplies (H, G) per side; second pair is a column permutation of the first.
function QuantumTannerViaLeftRightActions(group, A::Vector, B::Vector,
        H_A, G_A, H_B, G_B;
        p1 = collect(1:size(H_A,2)),
        p2 = collect(1:size(H_B,2)))
    QuantumTannerViaLeftRightActions(group, A, B,
        H_A, H_A[:,p1], G_A, G_A[:,p1],
        H_B, H_B[:,p2], G_B, G_B[:,p2])
end

# Simplest convenience constructor: only parity checks per side; generator matrices derived via dual_code.
function QuantumTannerViaLeftRightActions(group, A::Vector, B::Vector,
        H_A::Matrix, H_B::Matrix;
        p1 = collect(1:size(H_A,2)),
        p2 = collect(1:size(H_B,2)))
    QuantumTannerViaLeftRightActions(group, A, B,
        H_A, dual_code(H_A), H_B, dual_code(H_B); p1 = p1, p2 = p2)
end


parity_matrix_x(c::QuantumTannerViaLeftRightActions) = c.hx
parity_matrix_z(c::QuantumTannerViaLeftRightActions) = c.hz
parity_matrix_xz(c::QuantumTannerViaLeftRightActions) = (c.hx, c.hz)
code_n(c::QuantumTannerViaLeftRightActions) = size(c.hx, 2)
