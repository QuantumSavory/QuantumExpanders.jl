# QuantumExpanders.jl

```@meta
DocTestSetup = quote
    using QuantumClifford, QuantumExpanders
end
```

QuantumExpanders is a Julia package for constructing quantum Tanner codes.

# Comparison with existing work

Our method produces quantum Tanner codes with **significantly higher rates** than those reported in recent literature (Leverrier et al., *[Small quantum Tanner codes from left–right Cayley complexes](https://arxiv.org/pdf/2512.20532)*). For example:

- `[[252, 70, 6]]` → higher rate than `[[252, 2, 20]]`  
- `[[392, 96, 5]]` → higher rate than `[[396, 2, 29]]`  
- `[[576, 126, 7]]` → higher rate than `[[576, 28, 24]]`

While [recent work](https://arxiv.org/pdf/2512.20532) showcases codes with high distance, our method achieves higher rates, offering a complementary approach in the quantum Tanner code design.

## Quantum Tanner codes using other Frobenius groups

Expanding upon the work of [radebold2025explicit](https://arxiv.org/pdf/2508.05095), which was confined to
[dihedral groups](https://en.wikipedia.org/wiki/Dihedral_group), we have constructed new explicit quantum Tanner
codes based on a broader class of [Frobenius groups](https://en.wikipedia.org/wiki/Frobenius_group).

### Symmetric Group

Here is the `[[150, 48, 4]]` using [symmetric group](https://en.wikipedia.org/wiki/Symmetric_group) of order 3.

```julia
julia> rng = MersenneTwister(43);

julia> G = symmetric_group(3);

julia> S = normal_cayley_subset(G);

julia> hx, hz = random_quantum_Tanner_code(0.65, G, S, S, bipartite=false, use_same_local_code=true, rng=rng);
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

```julia
julia> rng = MersenneTwister(52);

julia> x = cperm([1,2,3,4]);

julia> G = permutation_group(4, [x]);

julia> S = normal_cayley_subset(G)
3-element Vector{PermGroupElem}:
 (1,2,3,4)
 (1,3)(2,4)
 (1,4,3,2)

julia> hx, hz = random_quantum_Tanner_code(0.65, G, S, S, bipartite=false, use_same_local_code=true, rng=rng);
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

```julia
julia> rng = MersenneTwister(54);

julia> G = cyclic_group(7);

julia> S = normal_cayley_subset(G);

julia> hx, hz = random_quantum_Tanner_code(0.7, G, S, S, bipartite=false, use_same_local_code=true, rng=rng);
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

```julia
julia> rng = MersenneTwister(42);

julia> G = small_group(8, 4);

julia> describe(G), small_group_identification(G)
("Q8", (8, 4))

julia> S = normal_cayley_subset(G);

julia> hx, hz = random_quantum_Tanner_code(0.72, G, S, S, bipartite=false, use_same_local_code=true, rng=rng);
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

Here is a `[[576, 126, 7]]` code using the [direct product](https://en.wikipedia.org/wiki/Direct_product_of_groups) of two cyclic groups (**C₃ × C₃**).

```julia
julia> rng = MersenneTwister(20);

julia> G = small_group(9, 2);

julia> describe(G), small_group_identification(G)
("C3 x C3", (9, 2))

julia> S = normal_cayley_subset(G);

julia> hx, hz = random_quantum_Tanner_code(0.8, G, S, S, bipartite=false, use_same_local_code=true, rng=rng);
(length(group), length(A), length(B)) = (9, 8, 8)
length(group) * length(A) * length(B) = 576
[ Info: |Q| = |G||A||B| = 576
Hᴬ = [1 1 1 0 1 1 1 1]
Hᴮ = [0 0 0 1 0 1 1 1; 1 1 0 1 1 1 0 0; 1 0 1 0 1 1 1 1; 1 0 1 0 0 0 1 0; 1 1 1 0 1 1 1 1; 1 0 1 1 1 0 1 0]
Cᴬ = [1 1 0 0 0 0 0 0; 1 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 1 0 0 0 1 0 0 0; 1 0 0 0 0 1 0 0; 1 0 0 0 0 0 1 0; 1 0 0 0 0 0 0 1]
Cᴮ = [1 1 0 0 0 0 0 0; 1 0 1 0 0 0 0 0; 0 0 0 1 0 0 0 0; 1 0 0 0 1 0 0 0; 1 0 0 0 0 1 0 0; 1 0 0 0 0 0 1 0; 1 0 0 0 0 0 0 1]
size(Cˣ) = (49, 64)
size(Cᶻ) = (1, 64)
r1 = rank(𝒞ˣ) = 441
r2 = rank(𝒞ᶻ) = 9

julia> c = CSS(hx, hz);

julia> import JuMP; import HiGHS;

julia> c = CSS(hx, hz);

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=120))
(576, 126, 7)
```

Here is the novel `[[360, 61, 10]]` quantum Tanner code constructed from [Morgenstern Ramanujan graphs](https://www.sciencedirect.com/science/article/pii/S0095895684710549)
for even prime power q.

```julia
julia> l = 1; i = 2;

julia> q = 2^l
2

julia> Δ = q+1
3

julia> SL₂, B = morgenstern_generators(l, i)
[ Info: |SL₂(𝔽(4))| = 60
(SL(2,4), Oscar.MatrixGroupElem{Nemo.FqFieldElem, Nemo.FqMatrix}[[o+1 o+1; 1 o+1], [o+1 1; o+1 o+1], [o+1 o; o o+1]])

julia> A = alternative_morgenstern_generators(B, FirstOnly())
4-element Vector{Oscar.MatrixGroupElem{Nemo.FqFieldElem, Nemo.FqMatrix}}:
 [0 1; 1 o+1]
 [o+1 1; 1 0]
 [o+1 o+1; o 0]
 [0 o+1; o o+1]

julia> rng = MersenneTwister(21);

julia> hx, hz = random_quantum_Tanner_code(0.74, SL₂, A, B, rng=rng);
(length(group), length(A), length(B)) = (60, 4, 3)
length(group) * length(A) * length(B) = 720
[ Info: |V₀| = |V₁| = |G| = 60
[ Info: |E_A| = Δ|G| = 240, |E_B| = Δ|G| = 180
[ Info: |Q| = Δ²|G|/2 = 360
Hᴬ = [0 1 1 0]
Hᴮ = [1 1 0; 0 1 1]
Cᴬ = [1 0 0 0; 0 1 1 0; 0 0 0 1]
Cᴮ = [1 1 1]
size(Cˣ) = (3, 12)
size(Cᶻ) = (2, 12)
r1 = rank(𝒞ˣ) = 179
r2 = rank(𝒞ᶻ) = 120

julia> c = CSS(hx, hz);

julia> import JuMP; import HiGHS;

julia> code_n(c), code_k(c), distance(c, DistanceMIPAlgorithm(solver=HiGHS, time_limit=120))
(360, 61, 10)
```
