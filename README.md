# QuantumExpanders.jl

<table>
    <tr>
        <td>Documentation</td>
        <td>
            <a href="https://quantumsavory.github.io/QuantumExpanders.jl/stable"><img src="https://img.shields.io/badge/docs-stable-blue.svg" alt="Documentation of latest stable version"></a>
            <a href="https://quantumsavory.github.io/QuantumExpanders.jl/dev"><img src="https://img.shields.io/badge/docs-dev-blue.svg" alt="Documentation of dev version"></a>
        </td>
    </tr><tr></tr>
    <tr>
        <td>Continuous integration</td>
        <td>
            <a href="https://github.com/QuantumSavory/QuantumExpanders.jl/actions?query=workflow%3ACI+branch%3Amaster"><img src="https://github.com/QuantumSavory/QuantumExpanders.jl/actions/workflows/ci.yml/badge.svg" alt="GitHub Workflow Status"></a>
            <a href="https://buildkite.com/quantumsavory/QuantumExpanders"><img src="https://badge.buildkite.com/8ef137151415f29c03544c5b7963f6bc6afc1022f29cfc072a.svg?branch=master" alt="Buildkite Workflow Status"></a>
        </td>
    </tr><tr></tr>
    <tr>
        <td>Code coverage</td>
        <td>
            <a href="https://codecov.io/gh/QuantumSavory/QuantumExpanders.jl"><img src="https://img.shields.io/codecov/c/gh/QuantumSavory/QuantumExpanders.jl?label=codecov" alt="Test coverage from codecov"></a>
        </td>
    </tr><tr></tr>
</table>

QuantumExpanders is a &nbsp;
    <a href="https://julialang.org">
        <img src="https://raw.githubusercontent.com/JuliaLang/julia-logo-graphics/master/images/julia.ico" width="16em">
        Julia Language
    </a>
    &nbsp; package for constructing quantum Tanner codes. To install QuantumExpanders,
    please <a href="https://docs.julialang.org/en/v1/manual/getting-started/">open
    Julia's interactive session (known as REPL)</a> and press the <kbd>]</kbd> key in the REPL to use the package mode, and then type:
</p>

```julia
pkg> add https://github.com/QuantumSavory/QuantumExpanders.jl.git
```

To update, just type `up` in the package mode.

- `random_quantum_Tanner_code` - Create a CSS code out of two Tanner codes 𝒞ᶻ and 𝒞ˣ, each constructed out of two related graphs 𝒢₀□, 𝒢₁□, the graphs build out of a Cayley complex, which itself was based on a group G=SL₂qⁱ and two generator sets A and B. For consistency of indexing of "qubits" ≈ "graph edges" ≈ "squares" and "local bits" ≈ "generator pairs" we have the `edge_*_idx` maps.

```jldoctest
julia> l = 1; i = 2;

julia> q = 2^l
2

julia> Δ = q+1
3

julia> SL₂, B = morgenstern_generators(l, i)
[ Info: |SL₂(𝔽(4))| = 60
(SL(2,4), MatrixGroupElem{FqFieldElem, FqMatrix}[[o+1 o+1; 1 o+1], [o+1 1; o+1 o+1], [o+1 o; o o+1]])

julia> A = alternative_morgenstern_generators(B, FirstOnly())
4-element Vector{MatrixGroupElem{FqFieldElem, FqMatrix}}:
 [0 1; 1 o+1]
 [o+1 1; 1 0]
 [o+1 o+1; o 0]
 [0 o+1; o o+1]

julia> rng = MersenneTwister(10);

julia> hx, hz = random_quantum_Tanner_code(0.4, SL₂, A, B, rng=deepcopy(rng));

julia> c = CSS(hx, hz);

julia> import JuMP; import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=120))
(360, 10, 4)
```

## Quantum Tanner codes using other Frobenius groups

Expanding upon the work of [radebold2025explicit](@cite), which was confined to dihedral groups, we have constructed new explicit quantum Tanner codes based on a broader class of Frobenius groups.

### Symmetric Group

Here is the `[[150, 48, 4]]` using symmetric group of order 3.

```jldoctest
julia> rng = MersenneTwister(43);

julia> hx, hz = random_quantum_Tanner_code(0.65, G, S, S, bipartite=false, use_same_local_code=true, rng=deepcopy(rng));
(length(group), length(A), length(B)) = (6, 5, 5)
length(group) * length(A) * length(B) = 150
[ Info: |Q| = |G||A||B| = 150
Hᴬ = [1 1 0 1 1]
Hᴮ = [1 1 0 1 0; 0 0 0 1 0; 1 1 0 1 1]
Cᴬ = [1 1 0 0 0; 0 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1]
Cᴮ = [1 1 0 0 0; 0 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1]
size(Cˣ) = (16, 25)
size(Cᶻ) = (1, 25)
r1 = rank(𝒞ˣ) = 96
r2 = rank(𝒞ᶻ) = 6

julia> c = CSS(hx, hz);

julia> import JuMP; import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=120))
(150, 48, 4)
```

### Permutation Group

Here is the `[[36, 16, 3]]` code based on permutation group of order 4 that improves upon the
[[36, 8, 3]] parameters presented in `Table 5` of [radebold2025explicit](@cite), achieving
twice the number of logical qubits while maintaining the same distance.

```jldoctest
julia> rng = MersenneTwister(52);

julia> x = cperm([1,2,3,4]);

julia> G = permutation_group(4, [x]);

julia> S = normal_cayley_subset(G)
3-element Vector{PermGroupElem}:
 (1,2,3,4)
 (1,3)(2,4)
 (1,4,3,2)

julia> hx, hz = random_quantum_Tanner_code(0.65, G, S, S, bipartite=false, use_same_local_code=true, rng=deepcopy(rng));
(length(group), length(A), length(B)) = (4, 3, 3)
length(group) * length(A) * length(B) = 36
[ Info: |Q| = |G||A||B| = 36
Hᴬ = [1 1 1]
Hᴮ = [1 0 1]
Cᴬ = [1 1 0; 1 0 1]
Cᴮ = [1 1 0; 1 0 1]
size(Cˣ) = (4, 9)
size(Cᶻ) = (1, 9)
r1 = rank(𝒞ˣ) = 16
r2 = rank(𝒞ᶻ) = 4

julia> c = CSS(hx, hz);

julia> import JuMP; import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=120))
(36, 16, 3)
```

### Cyclic Group

Here is the `[[252, 70, 6]]` using cyclic group of order 7.

```jldoctest
julia> rng = MersenneTwister(54);

julia> G = cyclic_group(7);

julia> S = normal_cayley_subset(G);

julia> hx, hz = random_quantum_Tanner_code(0.7, G, S, S, bipartite=false, use_same_local_code=true, rng=deepcopy(rng));
(length(group), length(A), length(B)) = (7, 6, 6)
length(group) * length(A) * length(B) = 252
[ Info: |Q| = |G||A||B| = 252
Hᴬ = [1 1 1 1 1 1]
Hᴮ = [1 0 0 1 1 0; 1 1 0 1 1 0; 0 1 1 0 1 0; 1 0 1 1 1 0]
Cᴬ = [1 1 0 0 0 0; 1 0 1 0 0 0; 1 0 0 1 0 0; 1 0 0 0 1 0; 1 0 0 0 0 1]
Cᴮ = [1 1 0 0 0 0; 1 0 1 0 0 0; 1 0 0 1 0 0; 1 0 0 0 1 0; 1 0 0 0 0 1]
size(Cˣ) = (25, 36)
size(Cᶻ) = (1, 36)
r1 = rank(𝒞ˣ) = 175
r2 = rank(𝒞ᶻ) = 7

julia> c = CSS(hx, hz);

julia> import JuMP; import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=120))
(252, 70, 6)
```