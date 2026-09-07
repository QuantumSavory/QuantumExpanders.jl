# [Quantum Tanner Codes via Left-Right Actions](@id quantum-tanner-left-right-actions)

[`QuantumTannerViaLeftRightActions`](@ref) constructs quantum Tanner (QT) codes
using the **lifted left-right action description**. Instead of building the code
through the square-complex geometry of the original LRCC construction, this
description lifts a small CSS template through the commuting left and right
regular actions of a finite group, following the lifting formulation of
[leverrier2025small](@cite).

This construction is particularly convenient for explicit code searches because:

- the two local-code lengths may be different;
- ``A`` and ``B`` are **multisets**, so repeated group elements are allowed;
- the total non-conjugacy condition required by the bipartite LRCC construction
  is not needed; and
- different column permutations of the local codes can be searched efficiently.

For the geometric LRCC description, see
[Quantum Tanner Codes](@ref quantum-tanner-codes).

## Construction

Let ``G`` be a finite group and let

```math
A=(a_1,\\ldots,a_{n_A}),
\\qquad
B=(b_1,\\ldots,b_{n_B})
```

be two multisets of elements of ``G``.

On the ``A`` side, choose two binary linear codes

```math
C_i=\\ker H_i=\\operatorname{Im}(G_i^\\mathsf{T}),
\\qquad i\\in\\{0,1\\},
```

of length ``n_A``. On the ``B`` side, choose

```math
C_i'=\\ker H_i'=\\operatorname{Im}(G_i'^\\mathsf{T}),
\\qquad i\\in\\{0,1\\},
```

of length ``n_B``.

The construction starts from a base CSS code on an
``n_A\\times n_B`` grid. Its parity-check matrices are

```math
H_X^{\\mathrm{base}}
=
\\begin{pmatrix}
H_0\\otimes G_0'\\
H_1\\otimes G_1'
\\end{pmatrix},
\\qquad
H_Z^{\\mathrm{base}}
=
\\begin{pmatrix}
G_0\\otimes H_1'\\
G_1\\otimes H_0'
\\end{pmatrix}.
```

The lift replaces every base qubit by a fiber of ``|G|`` qubits.

## Left and right regular actions

For ``a,b\\in G``, define permutation matrices ``L_a`` and ``R_b`` on the group
algebra ``\\mathbb{F}_2[G]`` by

```math
L_a e_g=e_{ag},
\\qquad
R_b e_g=e_{g b^{-1}}.
```

The inverse in the right action is chosen so that

```math
L_aR_b=R_bL_a
```

for all ``a,b\\in G``.

For the multisets ``A`` and ``B``, the implementation assembles block-diagonal
operators

```math
L_A
=
\\bigoplus_{i=1}^{n_A}
\\bigoplus_{j=1}^{n_B}
L_{a_i},
\\qquad
R_B
=
\\bigoplus_{i=1}^{n_A}
\\bigoplus_{j=1}^{n_B}
R_{b_j}.
```

Repeated elements in ``A`` or ``B`` simply produce repeated permutation blocks,
so multiplicities are allowed.

## Lifted parity-check matrices

Let ``I_{|G|}`` denote the identity matrix on the group-algebra fiber. The
lifted QT code has parity-check matrices

```math
H_X
=
\\begin{pmatrix}
H_0\\otimes G_0'\\otimes I_{|G|}\\
\\left(H_1\\otimes G_1'\\otimes I_{|G|}\\right)L_AR_B
\\end{pmatrix},
```

and

```math
H_Z
=
\\begin{pmatrix}
\\left(G_0\\otimes H_1'\\otimes I_{|G|}\\right)R_B\\
\\left(G_1\\otimes H_0'\\otimes I_{|G|}\\right)L_A
\\end{pmatrix}.
```

Because the left and right actions commute and each local pair satisfies

```math
H_iG_i^\\mathsf{T}=0,
```

the resulting matrices obey the CSS condition

```math
H_XH_Z^\\mathsf{T}=0.
```

The number of physical qubits is

```math
n=|G|\\,n_A n_B.
```

## Constructors

Three forms of [`QuantumTannerViaLeftRightActions`](@ref) are available.

### 1. Parity-check matrices only

```julia
QuantumTannerViaLeftRightActions(
    group,
    A,
    B,
    H_A,
    H_B;
    p1 = ...,
    p2 = ...,
)
```

This is the shortest constructor. Generator matrices are obtained internally
using [`dual_code`](@ref).

By default,

```julia
p1 = 1:size(H_A, 2)
p2 = 1:size(H_B, 2)
```

so the second local-code pair is identical to the first.

### 2. Parity-check and generator matrices

```julia
QuantumTannerViaLeftRightActions(
    group,
    A,
    B,
    H_A,
    G_A,
    H_B,
    G_B;
    p1 = ...,
    p2 = ...,
)
```

Use this form when the generator matrices are already available.

Internally, the second pair on each side is obtained by applying the column
permutations:

```math
H_1 = H_A[:,p_1],
\\qquad
G_1 = G_A[:,p_1],
```

and

```math
H_1' = H_B[:,p_2],
\\qquad
G_1' = G_B[:,p_2].
```

Thus `p1` acts on the ``A``-side local code and `p2` acts on the ``B``-side
local code.

### 3. Fully explicit local-code pairs

```julia
QuantumTannerViaLeftRightActions(
    group,
    A,
    B,
    H0_A,
    H1_A,
    G0_A,
    G1_A,
    H0_B,
    H1_B,
    G0_B,
    G1_B,
)
```

This is the most general form. Use it when ``H_1`` and ``H_1'`` are not obtained
from ``H_0`` and ``H_0'`` by simple column permutations.

The constructor verifies that every ``(H,G)`` pair is dual and checks the final
CSS orthogonality condition.

## Example

The following example constructs a ``[[288, 8, (≤15, ≤15)]]`` lifted QT code over ``S_3``.

```julia
julia> using QuantumExpanders, Oscar, QECCore

julia> G = codomain(isomorphism(PermGroup, small_group(6, 1)));

julia> A = [
           one(G),
           one(G),
           cperm(G, [2,3]),
           cperm(G, [2,3]),
           cperm(G, [1,2,3]),
           cperm(G, [1,2]),
           cperm(G, [1,2]),
           cperm(G, [1,3,2]),
       ];

julia> B = [
           one(G),
           one(G),
           cperm(G, [2,3]),
           cperm(G, [1,2,3]),
           cperm(G, [1,2]),
           cperm(G, [1,3,2]),
       ];

julia> H844 = [
           1 0 0 0 0 1 1 1;
           0 1 0 0 1 0 1 1;
           0 0 1 0 1 1 0 1;
           0 0 0 1 1 1 1 0
       ];

julia> G844 = [
           0 1 1 1 1 0 0 0;
           1 0 1 1 0 1 0 0;
           1 1 0 1 0 0 1 0;
           1 1 1 0 0 0 0 1
       ];

julia> H633 = [
           1 0 0 0 1 1;
           0 1 0 1 0 1;
           0 0 1 1 1 0
       ];

julia> G633 = [
           0 1 1 1 0 0;
           1 0 1 0 1 0;
           1 1 0 0 0 1
       ];

julia> c = QuantumTannerViaLeftRightActions(
           G,
           A,
           B,
           H844,
           G844,
           H633,
           G633;
           p1 = 1:8,
           p2 = [1,2,4,3,6,5],
       );

julia> code_n(c), code_k(c)
(288, 8)

julia> hx, hz = parity_matrix_xz(c);

julia> maximum(vec(sum(hx, dims=2))), maximum(vec(sum(hz, dims=2)))
(12, 12)
```

Here

```math
|G|=6,\\qquad n_A=8,\\qquad n_B=6,
```

so the blocklength follows immediately from

```math
n=6\\cdot8\\cdot6=288.
```

The nontrivial `p2` changes the second ``B``-side local-code pair and therefore
changes the resulting quantum code while keeping the same group and multisets.

## Multisets and repeated elements

Unlike the generating sets used by [`QuantumTannerCode`](@ref), the vectors
`A` and `B` in the lifted constructor are genuinely **multisets**.

For example,

```julia
A = [
    one(G),
    one(G),
    a,
    a,
    b,
    c,
]
```

is valid. The two copies of `one(G)` and the two copies of `a` correspond to
separate positions in the base code and therefore to separate blocks in
``L_A``.

No symmetry condition such as

```math
A=A^{-1}
```

is required, and no total non-conjugacy condition between ``A`` and ``B`` is
needed. CSS commutation instead follows from the commuting left and right
regular actions.

## Accessing the code

The resulting object implements the usual `QECCore` CSS-code interface:

```julia
julia> hx = parity_matrix_x(c);

julia> hz = parity_matrix_z(c);

julia> hx, hz = parity_matrix_xz(c);

julia> code_n(c)
288

julia> code_k(c)
8
```

The parity-check matrices may then be passed to the rest of the
`QECCore` / `QuantumClifford.ECC` ecosystem for distance estimation, decoding,
or other code analysis.

## Choosing between the two QT constructions

Use [`QuantumTannerViaLeftRightActions`](@ref) when:

- you want to work directly with the lifted algebraic description;
- ``A`` or ``B`` contain repeated elements;
- the ``A``- and ``B``-side local codes have different lengths;
- you want to search over column permutations `p1` and `p2`; or
- you do not want to impose the LRCC total non-conjugacy condition.

Use [`QuantumTannerCode`](@ref) when:

- you want the original square-complex geometry explicitly; or
- you already have symmetric LRCC generating sets satisfying TNC.

## API

See [`QuantumTannerViaLeftRightActions`](@ref) for the API reference.

## References

```@bibliography
Pages = ["quantum_tanner_left_right_actions.md"]
Canonical = false
```
