# [Morgenstern Ramanujan Graphs](@id morgenstern-graphs)

```@meta
DocTestSetup = quote
    using QuantumExpanders
    using Graphs
    using LinearAlgebra
end
```

Morgenstern's construction [morgenstern1994existence](@cite) provides explicit
``(q+1)``-regular Ramanujan graphs for **every even prime power** ``q = 2^l``,
extending the celebrated Lubotzky–Phillips–Sarnak construction
[lubotzky1988ramanujan](@cite), which requires ``q`` to be an odd prime.

A ``k``-regular graph is a **Ramanujan graph** if every eigenvalue ``\mu`` of its
adjacency matrix other than ``\pm k`` satisfies

```math
|\mu| \leq 2\sqrt{k-1}.
```

This is asymptotically optimal by the Alon–Boppana bound, which makes Ramanujan
graphs the best possible spectral expanders. Expander graphs of this kind are the
combinatorial backbone of the quantum Tanner code construction of
[dinur2022locally](@cite) implemented in this package.

## Construction

For ``q = 2^l`` and even ``i``, the Morgenstern graph is the Cayley graph

```math
\Gamma = \mathrm{Cay}\bigl(\mathrm{SL}_2(\mathbb{F}_{q^i}),\, B\bigr),
```

where ``B`` is a set of ``q+1`` involutions in ``\mathrm{SL}_2(\mathbb{F}_{q^i})``
arising from a quaternion algebra over the function field ``\mathbb{F}_q(x)``.
The resulting graph is ``(q+1)``-regular, connected, and non-bipartite, with

```math
|\Gamma| = q^{3i} - q^{i}.
```

The generator set is produced by [`morgenstern_generators`](@ref), and the Cayley
graph is assembled with [`cayley_right`](@ref).

## Example: ``l = 1,\ i = 2``

Here we construct the Morgenstern Ramanujan graph for `l = 1, i = 2`
(so ``q = 2^l = 2``) and verify that it satisfies **all** the properties
guaranteed by Theorem 5.13 of [morgenstern1994existence](@cite), as well as the
spectral expansion bounds of Claims 6.1 and 6.2 of [dinur2022locally](@cite).

```julia
julia> using QuantumExpanders, Graphs, LinearAlgebra

julia> l = 1; i = 2;

julia> q = 2^l; r = q + 1; # q = 2, so Γ is 3-regular

julia> G, B = morgenstern_generators(l, i);
[ Info: |SL₂(𝔽(4))| = 60

julia> Γ = cayley_right(G, B);
```

### Generator set ``B``

The set ``B`` contains ``q + 1`` generators, each of determinant ``1``
and order ``2`` (so ``\Gamma`` is an undirected simple graph):

```julia
julia> length(B) == q + 1
true

julia> all(det(b) == one(base_ring(b)) for b in B)
true

julia> all(matrix(b^2) == identity_matrix(base_ring(b), 2) for b in B)
true
```

### Property I: ``(q+1)``-regularity and order ``|\Gamma| = q^{3i} - q^{i}``

```julia
julia> all(degree(Γ, v) == q + 1 for v in vertices(Γ))
true

julia> nv(Γ) == q^(3i) - q^i == 60
true

julia> is_connected(Γ)
true
```

### Property II: Non-bipartiteness

```julia
julia> is_bipartite(Γ)
false
```

### Property III: Ramanujan bound

The trivial eigenvalue is ``q + 1``, and every other eigenvalue ``\mu``
satisfies ``|\mu| \leq 2\sqrt{q}``:

```julia
julia> λs = sort(real.(eigvals(Matrix(adjacency_matrix(Γ)))), rev=true);

julia> λs[1] ≈ q + 1
true

julia> all(abs(μ) ≤ 2√q + 1e-10 for μ in λs[2:end])
true
```

### Girth bound

The girth satisfies ``g(\Gamma) \geq \tfrac{2}{3}\log_q |\Gamma|``:

```julia
julia> using LibIGraph

julia> girth_lower_bound = floor(Int, (2/3)*log(q, nv(Γ)));

julia> g_igraph = IGraph(Γ); girth_val = Ref{LibIGraph.igraph_real_t}(0.0); cycle = IGVectorInt();

julia> LibIGraph.igraph_girth(g_igraph.objref, girth_val, cycle.objref);

julia> Int(girth_val[]) >= girth_lower_bound
true
```

### Property IV: Diameter bound

The diameter satisfies ``\mathrm{diam}(\Gamma) \leq 2\log_q |\Gamma| + 2``:

```julia
julia> diameter(Γ) ≤ ceil(Int, 2*log(q, nv(Γ)) + 2)
true
```

### Property V: Chromatic number

The chromatic number satisfies ``\chi(\Gamma) \geq \frac{q+1}{2\sqrt{q}} + 1``:

```julia
julia> χ_lower_bound = (q + 1)/(2√q) + 1;

julia> greedy_color(Γ; sort_degree=false, reps=1000).num_colors >= χ_lower_bound
true

julia> length(color(Γ; algorithm=DSATUR()).colors) >= χ_lower_bound
true

julia> length(color(Γ; algorithm=Greedy()).colors) >= χ_lower_bound
true
```

### Property VI: Independence number

The independence number satisfies ``i(\Gamma) \leq \frac{2\sqrt{q}}{q+1}|\Gamma|``:

```julia
julia> ind_set = independent_set(Γ, MaximalIndependentSet());

julia> length(ind_set) ≤ ceil(Int, (2√q/(q + 1))*nv(Γ))
true

julia> all(u == v || !has_edge(Γ, u, v) for u in ind_set, v in ind_set)
true
```

### Expander properties

``\Gamma`` is an ``(n, r, 1 - \lambda^2/r^2)``-expander with
``\lambda = \max_{\mu \neq q+1} |\mu|``, where ``\lambda \leq r - d^2/8r``
and ``\lambda \leq 2\sqrt{r-1}`` (optimality):

```julia
julia> λ = maximum(abs.(λs[2:end]));

julia> d = 1 - λ^2/r^2; d > 0
true

julia> λ ≤ r - d^2/(8r) + 1e-10
true

julia> λ ≤ 2√(r - 1) + 1e-10
true
```

## Spectral expansion of the alternative generator sets

The alternative generator sets produced by
[`alternative_morgenstern_generators`](@ref), which are used in the quantum Tanner
code construction, satisfy the explicit second-eigenvalue bounds of
[dinur2022locally](@cite):

- **Claim 6.1 (ii)** — the `AllPairs` generators yield a Cayley graph of degree
  ``k_1 = q^2 + q`` with normalized second eigenvalue
  ``\lambda_2 < \frac{3q-1}{q^2+q}``.
- **Claim 6.2** — the `FirstOnly` generators yield a Cayley graph of degree
  ``k_1 = 2q`` with normalized second eigenvalue
  ``\lambda_2 < \frac{3\sqrt{2q-1}}{2q}``.

```julia
julia> A₁ = alternative_morgenstern_generators(B, AllPairs());

julia> Γ₁ = cayley_right(G, A₁);

julia> λ₂ = sort(real.(eigvals(Matrix(adjacency_matrix(Γ₁))/(q^2 + q))), rev=true)[2];

julia> λ₂ < (3q - 1)/(q^2 + q)
true

julia> A₂ = alternative_morgenstern_generators(B, FirstOnly());

julia> Γ₂ = cayley_right(G, A₂);

julia> λ₂ = sort(real.(eigvals(Matrix(adjacency_matrix(Γ₂))/(2q))), rev=true)[2];

julia> λ₂ < 3√(2q - 1)/(2q)
true
```

## References

```@bibliography
Pages = ["morgenstern.md"]
Canonical = false
```
