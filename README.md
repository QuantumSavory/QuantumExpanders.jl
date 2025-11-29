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

- `random_quantum_Tanner_code` constructs a quantum CSS code by instantiating two classical
Tanner codes, 𝒞ᶻ and 𝒞ˣ, on the graphs 𝒢₀□ and 𝒢₁□ of a left-right Cayley complex [leverrier2022quantum](https://arxiv.org/pdf/2202.13641).
This complex is generated from a group G, and two generating sets A and B of sizes Δ_A and Δ_B, which,
need not satisfy the Total No-Conjugacy condition because of *quadripartite* construction of LRCC. 
The code's qubits bijectively correspond to the squares Q of the complex, with the `edge_*_idx` output
providing the essential mappings between qubit indices, graph edges in 𝒢₀□ and 𝒢₁□, and the local
coordinate sets A×B at each vertex. The Z-parity checks of the quantum code are defined as the
generators of 𝒞ᶻ = T(𝒢₀□, (C_A⊗C_B)⊥), enforcing local constraints from the dual tensor code at
each vertex of V₀. Similarly, the X-parity checks are generators of 𝒞ˣ = T(𝒢₁□, (C_A⊥⊗C_B⊥)⊥),
enforced at vertices of V₁. When the Cayley graphs Cay(G,A) and Cay(G,B) are Ramanujan and the
component codes C_A and C_B are randomly chosen with *robust* dual tensor properties, this
construction produces an asymptotically good quantum LDPC code with parameters [[n, Θ(n), Θ(n)]]. The
QT code implementation provides a simplified variant of the Panteleev-Kalachev quantum LDPC codes
[panteleev2022asymptoticallygoodquantumlocally](https://arxiv.org/pdf/2111.03654) and is related to
the locally testable code of [dinur2022locally](https://arxiv.org/pdf/2111.04808).

Here is the `[[360, 10, 4]]` quantum Tanner code constructed from [Morgenstern Ramanujan graphs](https://www.sciencedirect.com/science/article/pii/S0095895684710549)
for even prime power q.

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

Expanding upon the work of [radebold2025explicit](https://arxiv.org/pdf/2508.05095), which was confined to
[dihedral groups](https://en.wikipedia.org/wiki/Dihedral_group), we have constructed new explicit quantum Tanner
codes based on a broader class of [Frobenius groups](https://en.wikipedia.org/wiki/Frobenius_group).

### Symmetric Group

Here is the `[[150, 48, 4]]` using [symmetric group](https://en.wikipedia.org/wiki/Symmetric_group) of order 3.

```jldoctest
julia> rng = MersenneTwister(43);

julia> G = symmetric_group(3);

julia> S = normal_cayley_subset(G);

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

Here is the `[[36, 16, 3]]` code based on [permutation group](https://en.wikipedia.org/wiki/Permutation_group)
of order 4 that improves upon the `[[36, 8, 3]]` parameters presented in `Table 5` of [radebold2025explicit](@cite),
achieving twice the number of logical qubits while maintaining the same distance.

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

Here is the `[[252, 70, 6]]` using [cyclic group](https://en.wikipedia.org/wiki/Cyclic_group) of order 7.

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

### Small Groups

Here is a `[[392, 96, 5]]` code using [quaternion group](https://en.wikipedia.org/wiki/Quaternion_group) of order 8.

```jldoctest
julia> rng = MersenneTwister(42);

julia> G = small_group(8, 4);

julia> describe(G), small_group_identification(G)
("Q8", (8, 4))

julia> S = normal_cayley_subset(G);

julia> hx, hz = random_quantum_Tanner_code(0.72, G, S, S, bipartite=false, use_same_local_code=true, rng=deepcopy(rng));
(length(group), length(A), length(B)) = (8, 7, 7)
length(group) * length(A) * length(B) = 392
[ Info: |Q| = |G||A||B| = 392
Hᴬ = [0 1 0 1 1 1 1]
Hᴮ = [1 0 0 1 0 1 0; 0 0 1 0 0 0 1; 0 0 1 0 1 0 1; 0 0 0 1 0 0 1; 1 1 1 1 0 1 0]
Cᴬ = [1 0 0 0 0 0 0; 0 0 1 0 0 0 0; 0 1 0 1 0 0 0; 0 1 0 0 1 0 0; 0 1 0 0 0 1 0; 0 1 0 0 0 0 1]
Cᴮ = [1 0 0 0 0 0 0; 0 0 1 0 0 0 0; 0 1 0 1 0 0 0; 0 1 0 0 1 0 0; 0 1 0 0 0 1 0; 0 1 0 0 0 0 1]
size(Cˣ) = (36, 49)
size(Cᶻ) = (1, 49)
r1 = rank(𝒞ˣ) = 288
r2 = rank(𝒞ᶻ) = 8

julia> c = CSS(hx, hz);

julia> import JuMP; import HiGHS;

julia> c = CSS(hx, hz);

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=120))
(392, 96, 5)
```
