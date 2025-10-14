"""
Computes the [Legendre symbol](https://en.wikipedia.org/wiki/Legendre_symbol) ``\\frac{a}{p}``
for an odd prime ``p``.
    
The Legendre symbol determines whether ``a`` is a [quadratic residue](https://en.wikipedia.org/wiki/Quadratic_residue)
modulo ``p``, and specifically controls which group ``\\mathrm{PGL}`` or ``\\mathrm{PSL}`` is used in the LPS construction
of the Ramanujan graph [lubotzky1988ramanujan](@cite).

- When ``\\frac{p}{q} = -1``, the graph ``X^{p,q}`` is constructed as a Cayley graph of ``\\mathrm{PGL}(2, \\mathbb{Z}/q\\mathbb{Z})``
and is *bipartite* of order ``q(q^2-1)``.
- When ``\\frac{p}{q} = 1``, the graph is constructed as a Cayley graph of ``\\mathrm{PSL}(2, \\mathbb{Z}/q\\mathbb{Z})``, is
*non-bipartite*, and has order ``q(q^2-1)/2``.
"""
function legendre_symbol(a::Int, p::Int)
    @assert is_prime(p) "p must be prime"
    ls = powermod(a, (p-1)÷2, p)
    return ls == p - 1 ? -1 : ls
end

"""
Generates the [center](https://en.wikipedia.org/wiki/Center_(group_theory)) of the general linear group
``\\mathrm{GL}(2, F)`` over a finite field ``F``, consisting of all scalar matrices

```math
\\begin{aligned}
\\begin{bmatrix}
x & 0 \\\\
0 & x
\\end{bmatrix}
\\end{aligned}
```

where ``x`` is any nonzero element of ``F``. These matrices form the center ``Z(GL(2, F))`` and
are used in the LPS construction to form the projective general linear group ``PGL(2, F) = GL(2, F)/Z(GL(2, F))``.
"""
function scalar_matrices_GL(GL2::MatrixGroup)
    F = base_ring(GL2)
    return [GL2([x 0; 0 x]) for x in F if x != 0]
end

"""
Generates the [center](https://en.wikipedia.org/wiki/Center_(group_theory)) of the special linear group
``\\mathrm{SL}(2, F)`` over a finite field ``F``, consisting of scalar matrices

```math
\\begin{aligned}
\\begin{bmatrix}
x & 0 \\\\
0 & x
\\end{bmatrix}
\\end{aligned}
```

with ``x^2 = 1``. These matrices form the center ``Z(SL(2, F))`` and are used in the LPS construction
to form the projective special linear group ``PSL(2, F) = SL(2, F)/Z(SL(2, F))``.
"""
function scalar_matrices_SL(SL2::MatrixGroup)
    F = base_ring(SL2)
    return [SL2([x 0; 0 x]) for x in F if x^2 == one(F)]
end

"""
Finds all integer solutions to the equation ``p = a^2 + b^2 + c^2 + d^2`` for a prime ``p \\equiv 1 \\pmod{4}``.
according to the [Jacobi's theorem](https://en.wikipedia.org/wiki/Jacobi%27s_four-square_theorem)."""
function solve_four_squares(p::Int)
    @assert mod(p, 4) == 1 "p must be ≡ 1 mod 4 [lubotzky1988ramanujan](@cite)"
    solutions = Tuple{Int,Int,Int,Int}[]
    max_val = isqrt(p)
    for w in 0:max_val, x in 0:max_val, y in 0:max_val
        rem = p - (w^2 + x^2 + y^2)
        rem ≥ 0 || continue
        z = isqrt(rem)
        z^2 == rem || continue
        signs_w = w == 0 ? [0] : [1, -1]
        signs_x = x == 0 ? [0] : [1, -1]
        signs_y = y == 0 ? [0] : [1, -1]
        signs_z = z == 0 ? [0] : [1, -1]
        for sw in signs_w, sx in signs_x, sy in signs_y, sz in signs_z
            push!(solutions, (sw*w, sx*x, sy*y, sz*z))
        end
    end
    @assert length(solutions) == 8*(p+1) "Jacobi's theorem: should have exactly 8(p+1) solutions [lubotzky1988ramanujan](@cite), Section 2, Page 264"
    return unique(solutions)
end

"""
Filters the solutions from `solve_four_squares` to select exactly ``p+1` solutions
with ``a > 0`` and ``b, c, d`` even. """
function process_solutions(solutions::AbstractVector, p::Int)
    @assert mod(p, 4) == 1 "p must be ≡ 1 mod 4 [lubotzky1988ramanujan](@cite)"
    @assert length(solutions) == 8*(p+1) "Jacobi's theorem: should have exactly 8(p+1) solutions [lubotzky1988ramanujan](@cite), Section 2, Page 264"
    filtered = filter(sol -> sol[1] > 0 && all(iseven, sol[2:4]), solutions)
    @assert all(sol -> sol[1] > 0, filtered) "All solutions must have a₀ > 0. Page 262 of [lubotzky1988ramanujan](@cite)"
    @assert all(sol -> isodd(sol[1]), filtered) "a₀ must be odd. Page 262 of [lubotzky1988ramanujan](@cite)"  
    @assert all(sol -> all(iseven, sol[2:4]), filtered) "a₁,a₂,a₃ must be even. Page 262 of [lubotzky1988ramanujan](@cite)"
    @assert length(filtered) == p + 1 "Should have exactly p+1 solutions. Page 262 of [lubotzky1988ramanujan](@cite)"
    @assert all(sol -> isodd(sol[1]), filtered) "After unit action, all a₀ must be odd"
    @assert all(sol -> all(iseven, sol[2:4]), filtered) "After unit action, all a₁,a₂,a₃ must be even"
    # See Section 3, Page 265 of [lubotzky1988ramanujan](@cite)
    # The set S containing p+1 elements splits into distinct conjugate pairs {α₁,ᾱ₁,α₂,ᾱ₂,...,αₛ,ᾱₛ} where s=(p+1)/2.
    conjugates = Set([(sol[1], -sol[2], -sol[3], -sol[4]) for sol in filtered])
    @assert length(conjugates) == length(filtered) "Solutions must form distinct conjugate pairs"
    @assert iseven(length(filtered)) "Number of solutions must be even to form conjugate pairs"
    return filtered
end

"""
Constructs the generator matrices from the filtered four-square solutions. For each solution (a, b, c, d), creates a matrix:

```math
\\begin{aligned}
\\begin{bmatrix}
a + i b & c + i d \\\\
-c + i d & a - i b
\\end{bmatrix}
\\end{aligned}
```

where ``i`` satisfies ``i^2 \\equiv -1 \\pmod{q}``. These matrices have determinant
``a^2 + b^2 + c^2 + d^2 = p`` and will serve as the ``p+1`` generators for the Cayley
graph.
"""
function lps_generators(solutions::AbstractVector, F::FqField, p::Int)
    u = sqrt(F(-1))  # Find u such that u² = -1 in F.
    @assert u^2 == F(-1)
    generators = MatrixElem{typeof(F(0))}[]
    for (a, b, c, d) in solutions
        mat = matrix(F, [a + u*b   c + u*d;
                        -c + u*d   a - u*b])
        @assert det(mat) == F(p) "Generator matrix must have determinant p"
        push!(generators, mat)
    end
     @assert length(generators) == p + 1 "Should have p+1 generators"
    return generators
end

"""
Check the Ramanujan property:
- For a (p+1)-regular graph, the trivial eigenvalue is p+1.
- All other eigenvalues should have absolute value ≤ 2√p.
"""
function is_ramanujan(g::SimpleGraph, p::Int)
    A = adjacency_matrix(g)
    λ = sort(eigvals(Matrix(A)), rev=true)
    bound = 2 * sqrt(p)
    non_trivial = filter(x -> abs(x - (p+1)) > 1e-6, λ)
    return all(v -> abs(v) ≤ bound + 1e-6, non_trivial)
end

function lps_graph(::Val{-1}, p::Int, q::Int)
    F = GF(q)
    GL2 = GL(2, F)
    center = scalar_matrices_GL(GL2)
    PG, morphism = quo(GL2, center)
    solutions = solve_four_squares(p)
    solutions_processed = process_solutions(solutions, p)
    generators = lps_generators(solutions_processed, F, p)
    gl_gens = [GL2(mat) for mat in generators]
    PG_generators = [morphism(gen) for gen in gl_gens]
    return cayley_right(PG, PG_generators)
end

function lps_graph(::Val{1}, p::Int, q::Int)
    F = GF(q)
    SL2 = SL(2, F)
    center = scalar_matrices_SL(SL2)
    PG, morphism = quo(SL2, center)
    solutions = solve_four_squares(p)
    solutions_processed = process_solutions(solutions, p)
    generators = lps_generators(solutions_processed, F, p)
    s = sqrt(F(p))
    s_inv = inv(s)
    generators_scaled = [s_inv * mat for mat in generators]
    sl_gens = [SL2(mat) for mat in generators_scaled]
    PG_generators = [morphism(gen) for gen in sl_gens]
    return cayley_right(PG, PG_generators)
end

"""
Constructs the Lubotzky–Phillips–Sarnak Ramanujan graph ``X^(p,q)`` as described in [lubotzky1988ramanujan](@cite).

Returns the ``(p+1)``-regular LPS Ramanujan graph ``X^{p,q}``.

# LPS Ramanujan graph ``X^(p,q)``

Let ``p`` and ``q`` be distinct primes congruent to ``1 \\pmod{4}``. The LPS Ramanujan graphs ``X^{p,q}`` are
``p+1``-regular Cayley graphs of the projective linear group ``\\mathrm{PSL}(2,\\mathbb{Z}/q\\mathbb{Z})`` when
the Legendre symbol satisfies ``\\left( \\tfrac{p}{q} \\right) = 1``, and of ``\\mathrm{PGL}(2,\\mathbb{Z}/q\\mathbb{Z})``
when ``\\left( \\tfrac{p}{q} \\right) = -1``.

The construction of these graphs relies on representing the prime p as a sum of four squares. Specifically, the p+1 generators
for the Cayley graph are derived from the distinct integer solutions to

```math
\\begin{aligned}
p = a_0^2 + a_1^2 + a_2^2 + a_3^2
\\end{aligned}
```

where ``a_0 > 0`` is odd and ``a_1, a_2, a_3`` are even. That there are exactly p+1 such representations is a consequence of
[Jacobi’s theorem](https://en.wikipedia.org/wiki/Jacobi%27s_four-square_theorem) on the number of representations ``r_4(n)``:

```math
\\begin{aligned}
r_4(n) = 8 \\sum_{\\substack{d \\mid n \\\\ 4 \\nmid d}} d
\\end{aligned}
```

To generalize this construction, [lubotzky1988ramanujan](@cite) used representations by certain quaternary quadratic forms
([dickson1927quaternary](@cite), [sarnak1990some](@cite)). Define the form

```math
\\begin{aligned}
\\mathcal{Q}_q(x_1, x_2, x_3, x_4) = x_1^2 + 4q^2x_2^2 + 4q^2x_3^2 + 4q^2x_4^2,
\\end{aligned}
```

and let ``r_Q(n)`` denote the number of integer solutions ``v \\in \\mathbb{Z}^4`` to ``\\mathcal{Q}_q(v) = n``. Unlike the
explicit formula for ``r_4(n)``, no simple closed form exists for ``r_Q(n)``. However, the Ramanujan conjecture [ramanujan1916certain](@cite)
proved in this context by Eichler [eichler1954quaternare](@cite) and Igusa [igusa1956fibre](@cite) provides a asymptotic approximation which
establishes LPS graphs ``X^{p,q}`` as Ramanujan graphs. For n = p^k with ``vk \\ge 0``v, we have

```math
\\begin{aligned}
r_Q(p^k) = C(p^k) + O_\\varepsilon\\!\\left( p^{k(1/2 + \\varepsilon)} \\right) \\quad \\text{as } k \\to \\infty, \\quad \\forall \\varepsilon > 0,
\\end{aligned}
```

where the main term C(p^k) is given by

```math
\\begin{aligned}
C(p^k) =
\\begin{cases}
c_1 \\displaystyle\\sum_{d \\mid p^k} d & \text{if } \\left( \\tfrac{p}{q} \\right) = 1 \\\\
c_2 \\displaystyle\\sum_{d \\mid p^k} d & \text{if } \\left( \\tfrac{p}{q} \\right) = -1 \\text{ and } k \text{ is even} \\\\
0 & \\text{if } \\left( \\tfrac{p}{q} \\right) = -1 \\text{ and } k \\text{ is odd}
\\end{cases}
\\end{aligned}
```

The constants ``c_1`` and ``c_2`` are determined in Section 4 of [lubotzky1988ramanujan](@cite).

!!! note
    As stated by [lubotzky1988ramanujan](@cite): "As was proved in Section 3 the girth ``g(X^{p,q})`` of our graphs
    ``\\to \\infty`` as q (and hence n) ``\\to \\infty``. Thus the spectrum of the graphs ``X^{p,q}`` lies in 
    ``[-2\\sqrt{p}, 2\\sqrt{p}]`` (besides ``\\pm(p+1)``) and it is distributed in this interval according to the
    density ``d\\mu_{p+1}`` as ``q \\to \\infty``".
    
Here, ``d\\mu_{p+1}`` denotes the limiting spectral distribution given in Proposition 4.3 of [lubotzky1988ramanujan](@cite) where
``k = p + 1`` is the degree of the graph.

```math
\\begin{aligned}
d\\mu_{p+1}(t) = 
\\begin{cases} 
\\dfrac{\\sqrt{p - t^2/4}}{\\pi(p+1)(1-(t/(p+1))^2)} \\, dt & \\text{if } |t| \\leq 2\\sqrt{p} \\\\
0 & \\text{otherwise}
\\end{cases}
\\end{aligned}
```

# Cayley graphs of Free Groups via Geometric Group Theory

As detailed in [lubotzky1988ramanujan](@cite), for a prime ``p \\equiv 1 \\pmod{4}``,
there exists a set S of p + 1 integral quaternions of norm p, unique up to units
and satisfying ``\\alpha \\equiv 1 \\pmod{2}``. [lubotzky1988ramanujan](@cite) establishes
that every quaternion ``\\alpha \\in ``H(\\mathbb{Z})`` with ``N(\\alpha) = p^k``
can be expressed uniquely in the form ``\\alpha = \\varepsilon p^r R_m(\\alpha_1, \\ldots, \\bar{\alpha}_s)``
where ``\\varepsilon`` is a unit, ``2r + m = k``, and ``R_m`` is a *reduced word*
in the elements of ``S`` and their conjugates, where "reduced" means no generator
is adjacent to its inverse (see Definition 2.3.4 of [loh2017geometric](@cite)). This
unique factorization property shows that the multiplicative group ``A(2)`` formed by
these quaternion classes (modulo the identification ``\\pm p^{v_1}\\alpha \\sim p^{v_2}\\beta``
is *freely generated* by the classes ``[\\alpha_1], [\\alpha_2], \\ldots, [\\alpha_s]``
(see Proposition 2.3.5 of [loh2017geometric](@cite)), with every element represented by
exactly one reduced word (see Corollary 2.3.6 of [loh2017geometric](@cite)). Consequently,
the Cayley graph ``\\text{Cay}(A(2), S)`` is an infinite (p + 1)-regular tree (see Theorem 2.3.1
of [loh2017geometric](@cite)). The finite Ramanujan graphs ``X^{p,q}`` are then constructed as
explicit quotients of this tree by taking the congruence subgroup ``A(2q)``, defined as the
kernel of reduction modulo ``2q``. The quotient ``A(2)/A(2q)`` is shown to be isomorphic to either
``\\PGL(2, \\mathbb{Z}/q\\mathbb{Z})`` or ``\\PSL(2, \\mathbb{Z}/q\\mathbb{Z})``, and the
resulting Cayley graph is the LPS graph ``X^{p,q}``.

### Arguments
- `p`: A prime number congruent to ``1 \\pmod{4}``.
- `q`: A prime number, distinct from `p`, also congruent to ``1 \\pmod{4}``.
"""
function LPS(p::Int, q::Int)
    @assert is_prime(p) && is_prime(q) "p and q must be primes."
    @assert p % 4 == 1 && q % 4 == 1 "p and q must be ≡ 1 mod 4."
    symbol = legendre_symbol(p, q)
    symbol == -1 && return lps_graph(Val(-1), p, q)
    symbol == 1 && return lps_graph(Val(1), p, q)
    throw(ArgumentError("Unexpected Legendre symbol value."))
end


"""
Constructs the edge-vertex incidence graph of a given graph G.
Let G = (V, E) be a graph with vertex set V and edge set E. The
edge-vertex incidence graph of G is a bipartite graph with vertex
set E ∪ V and edge set {(e, v) ∈ E × V : v is an endpoint of e}.

# Example

```jldoctest
julia> G = LPS(5, 13);

julia> B = edge_vertex_incidence_graph(G);

julia> is_unbalanced_bipartite(B)
true
```
"""
function edge_vertex_incidence_graph(G)
    n = nv(G)
    m = ne(G)
    B = SimpleGraph(n + m)
    edge_to_index = Dict{Graphs.Edge, Int}()
    for (i, e) in enumerate(edges(G))
        edge_to_index[e] = n + i
    end
    for (e, idx) in edge_to_index
        add_edge!(B, src(e), idx)
        add_edge!(B, dst(e), idx)
    end
    return B
end

"""Verify that the edge-vertex incidence graph B is an unbalanced bipartite graph."""
function is_unbalanced_bipartite(B)
    # The first n vertices in B correspond to the vertices of G
    n = 0
    while n + 1 <= nv(B) && degree(B, n + 1) > 2
        n += 1
    end
    # If no vertices satisfy degree > 2, return false
    if n == 0
        return false
    end
    # The degree of any vertex in the first n vertices of B is d
    d = degree(B, 1)
    # The remaining vertices in B correspond to the edges of G
    m = nv(B) - n
    for v in 1:n
        if degree(B, v) != d
            return false
        end
    end
    for v in (n + 1):(n + m)
        if degree(B, v) != 2
            return false
        end
    end
    if n * d != m * 2
        return false
    end
    return true
end

"""
Test the Alon–Chung Lemma [alon1988explicit](@cite) on a regular graph.

Returns whether the bound holds, edges in induced subgraph, and theoretical upper bound.

The Alon–Chung Lemma states that for a d-regular graph G with n vertices and 
second-largest eigenvalue ``\\lambda``, any vertex subset ``X \\subseteq V(G)``
with cardinality ``|X| = \\gamma n`` satisfies:

```math
\\begin{aligned}
|E(G[X])| \\leq \\frac{dn}{2} \\left( \\gamma^2 + \\frac{\\lambda}{d} \\gamma (1 - \\gamma) \\right)
\\end{aligned}
```

where ``E(G[X])`` denotes the edge set of the [induced subgraph](https://en.wikipedia.org/wiki/Induced_subgraph).

# Example

```jldoctest
julia> G = LPS(5, 13);

julia> B = edge_vertex_incidence_graph(G);

julia> alon_chung_lemma(B, 0.1)
(true, 143, 1520.6609615130378)
```

See [expandercode556667](@cite) for more details how Alon–Chung Lemma is used to characterize the
edge-vertex incidence graphs of the LPS expander graphs in the construction of classical expander codes.

### Arguments
- `g`: A d-regular graph
- `γ`: Fraction of vertices in the subset (0 < γ < 1)
"""
function alon_chung_lemma(g::AbstractGraph, γ::Real)
    n = nv(g)
    d = degree(g, 1)
    A = adjacency_matrix(g)
    Λ = eigvals(Matrix(A))
    λ = sort(Λ, rev=true)[2]
    # Select a random subset X of size γ*n.
    s = Int(round(γ*n))
    X = randperm(n)[1:s]
    subg = induced_subgraph(g, X)[1]
    # Count the number of edges in the subgraph induced by X.
    edges = ne(subg)
    bound = (d*n/2)*(γ^2+(λ/d)*γ*(1-γ))
    return edges ≤ bound, edges, bound
end
