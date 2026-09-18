# [Reproducing the manuscript instances](@id paper-instances)

The code tables in *Quantum Tanner Codes at Moderate Blocklength*
[mian2026quantum](@cite) specify each QT code through a finite group, group
elements, local classical codes, and in the lifted description via column
permutations. This page explains how those input data correspond to `QuantumExpanders.jl`
and what the package tests deterministically.

## What specifies an instance?

For a lifted code, the complete deterministic input is:

1. the finite group ``G``;
2. the ordered multisets ``A`` and ``B``;
3. the local parity-check and generator matrices on the two sides; and
4. the column permutations `p1` and `p2`.

The order and multiplicity of the entries in `A` and `B` matter. Two equal group
elements in different positions correspond to two different coordinates of the
base code.

For an LRCC code, the complete input is the group, the two LRCC generating sets,
and the two local parity-check/generator pairs. The generating sets must satisfy
the LRCC conditions described in [Quantum Tanner Codes](@ref quantum-tanner-codes).

## Worked lifted example

The following reconstructs the manuscript's ``[[768,24,(\leq16,\leq16)]]``
instance. In particular, its corrected second column permutation is
`[1,2,3,4,5,6,8,7]`.

```julia
using QuantumExpanders, Oscar

G = codomain(isomorphism(PermGroup, small_group(12, 1)))
x = cperm(G, [5,6,7])
y = cperm(G, [1,4,3,2], [6,7])

A = [
    one(G), one(G),
    x, x,
    y, y*x^2, y*x, inv(y),
]

H = [
    1 0 0 0 0 1 1 1;
    0 1 0 0 1 0 1 1;
    0 0 1 0 1 1 0 1;
    0 0 0 1 1 1 1 0
]

code = QuantumTannerViaLeftRightActions(
    G,
    A,
    A,
    H,
    H;
    p1 = 1:8,
    p2 = [1,2,3,4,5,6,8,7],
)

code_n(code), code_k(code)
# (768, 24)
```

The structural checks used throughout the regression suite are easy to repeat:

```julia
hx, hz = parity_matrix_xz(code)

@assert iszero(mod.(hx * hz', 2))
@assert maximum(vec(sum(hx, dims=2))) == 16
@assert maximum(vec(sum(hz, dims=2))) == 16
```

Changing only `p2` to `[1,2,3,4,5,7,8,6]` produces a different code with ``k=16``. This is why permutations are part of the reproducibility data rather than presentation-only metadata.

## Regression-test coverage

The repository contains deterministic reconstructions of the manuscript
catalogues:

| Manuscript collection | Instances | Regression file |
|---|---:|---|
| Main lifted-QT table | 30 | `test/test_lifted_quantum_tanner_table1.jl` |
| Codes near or beyond the ``\sqrt n`` line | 8 | `test/test_lifted_quantum_tanner_table2.jl` |
| LRCC appendix | 80 | `test/test_quantum_tanner_lrcc_appendix_explicit.jl` |
| Morgenstern examples | 2 | `test/test_quantum_tanner_morgenstern_appendix.jl` |

For each deterministic construction, the relevant tests check some or all of:

- `code_n(code)` and `code_k(code)`;
- CSS orthogonality;
- stabilizer rank; and
- maximum ``X``- and ``Z``-check weights.

The constructor regression tests also cover Oscar finite-field matrices and the
matrices returned by [`dual_code`](@ref).

## Interpreting the reported distances

The manuscript writes the lifted-code parameters as

```math
[[n,k,(\leq d_X,\leq d_Z)]].
```

Randomized code distancs are reported as **upper bounds** on the true distances. More trials can tighten an upper bound.

!!! warning "Do not convert an estimate into an equality"
    If an estimator returns `(dx, dz) = (21, 18)`, report the code as
    ``[[n,k,(\leq21,\leq18)]]`` unless an exact distance computation or a
    matching lower bound proves equality.

## A reproducible workflow

For a new manuscript row:

1. record the constructor inputs in one deterministic data structure;
2. construct the code and assert $n$, $k$, orthogonality, rank, and weights;
3. run randomized distance estimation separately and record the estimator,
   trial count, seed when applicable, and returned pair; and
4. generate or copy the manuscript row from the same canonical data whenever
   possible.

## References

```@bibliography
Pages = ["paper_instances.md"]
Canonical = false
```
