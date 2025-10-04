"""
Finds an irreducible polynomial of the form ``f(x) = x^2 + x + \\varepsilon`` over
the finite field ``\\mathbb{F}_q``.

For quadratic polynomials, irreducibility is equivalent to having no roots in the base field:

```jldoctest
julia> using Oscar

julia> 𝔽₂, _ = finite_field(2, 1, "a");

julia> R, x = polynomial_ring(𝔽₂, "x");

julia> f = x^2 + x + 𝔽₂(1);

julia> is_irreducible(f) && isempty(roots(f))
true

julia> g = x^2 + x + 𝔽₂(0);

julia> is_irreducible(g) || !isempty(roots(g))
true
```

In characteristic 2, the standard quadratic form ``x^2 + \\varepsilon`` is insufficient because
its derivative is zero, making it inseparable if it has a root. The form ``x^2 + x + \\varepsilon``
is used instead, as its derivative is 1, ensuring separability. Its irreducibility over ``\\mathbb{F}_q``
is equivalent to it having no roots in ``\\mathbb{F}_q``.

```jldoctest
julia> using Oscar

julia> 𝔽₂, _ = finite_field(2, 1, "a");

julia> R, x = polynomial_ring(𝔽₂, "x");

julia> f = x^2 + 𝔽₂(1);

julia> derivative(f) == 0 && !is_separable(f)
true

julia> 𝔽₄, _ = finite_field(2, 2, "b");

julia> R₄, y = polynomial_ring(𝔽₄, "y");

julia> factor(y^2 + 𝔽₄(1))
1 * (y + 1)^2
```

The form ``x^2 + x + \\varepsilon`` works correctly in characteristic 2, making it
essential for applications like the construction of optimal expander graphs:

```jldoctest
julia> using Oscar

julia> 𝔽₂, _ = finite_field(2, 1, "a");

julia> R, x = polynomial_ring(𝔽₂, "x");

julia> f = x^2 + x + 𝔽₂(1);

julia> is_irreducible(f) && is_separable(f) && derivative(f) == 1
true
```

# Morgenstern's construction of Ramanujan graphs

The construction requires a [quaternion algebra](https://en.wikipedia.org/wiki/Quaternion_algebra)
over ``\\mathbb{F}_q(x)`` of the form [morgenstern1994existence](@cite):

```math
\\begin{aligned}
\\mathcal{A} = k\\mathbf{1} + k\\mathbf{i} + k\\mathbf{j} + k\\mathbf{ij}, \\quad
\\mathbf{i}^2 = \\mathbf{i} + \\varepsilon, \\quad \\mathbf{j}^2 = x, \\quad \\mathbf{ij} = \\mathbf{ji} + \\mathbf{j}
\\end{aligned}
```

where ``k = \\mathbb{F}_q(x)`` and ``f(x) = x^2 + x + \\varepsilon`` is irreducible over ``\\mathbb{F}_q``.

The norm in this algebra is given by:

```math
\\begin{aligned}
N(a + b\\mathbf{i} + c\\mathbf{j} + d\\mathbf{ij}) = a^2 + b^2\\varepsilon + ab + (c^2 + d^2\\varepsilon + cd)x
\\end{aligned}
```

The "basic norm" elements that generate the Ramanujan graphs are exactly those of the form [morgenstern1994existence](@cite):

```math
\\begin{aligned}
\\xi = 1 + \\gamma\\mathbf{j} + \\delta\\mathbf{ij}, \\quad \\text{where } \\gamma, \\delta \\in \\mathbb{F}_q
\\end{aligned}
```

satisfy

```math
\\begin{aligned}
\\gamma^2 + \\gamma\\delta + \\delta^2\\varepsilon = 1
\\end{aligned}
```

This equation has exactly ``q+1`` solutions in ``\\mathbb{F}_q``, providing the ``q+1`` generators needed 
for ``(q+1)``-regular Ramanujan graphs.

Under the isomorphism ``\\theta: \\mathcal{A} \\to M_2(k)``, these generators map to matrices [morgenstern1994existence](@cite):

```math
\\begin{pmatrix}
1 & \\gamma + \\delta i \\\\
(\\gamma + \\delta i + \\delta)x & 1
\\end{pmatrix}
```

where ``i`` is a root of ``x^2 + x + \\varepsilon = 0``.

In characteristic 2, the polynomial ``x^2 + \\varepsilon`` is inseparable because its derivative is
zero. As a result, it cannot be used to define a separable quadratic field extension, which is necessary
for constructing a division quaternion algebra in Morgenstern’s Ramanujan graph construction. In contrast, 
``x^2 + x + \\varepsilon``, being both irreducible and separable, yields a proper non-split quaternion algebra.
Without separability, the resulting graphs do not attain the Ramanujan eigenvalue bound ``2\\sqrt{q}``,
violating the optimality of the construction.

# Arguments
- `R`: Polynomial ring ``\\mathbb{F}_q[x]`` used to construct the irreducible polynomial for Morgenstern's quaternion algebra.

"""
function morgenstern_f(R::FqPolyRing)
    x = gen(R)
    count = 0
    while true
        count += 1
        ε = rand(base_ring(R))
        p = x^2+x+ε
        if is_irreducible(p)
            @debug "Found irreducible polynomial after $count attempts: $p"
            return p
        end
    end
end

"""
Find all ``q + 1`` solutions ``(\\gamma, \\delta) \\in \\mathbb{F}_q^2`` to the
``\\gamma^2 + \\gamma\\delta + \\delta^2\\varepsilon = 1``.

## Quaternion Algebra

A [quaternion algebra](https://en.wikipedia.org/wiki/Quaternion_algebra) over ``k = \\mathbb{F}_q(x)`` is
a [skewfield](https://en.wikipedia.org/wiki/Division_ring) ``\\mathcal{A}`` with [center](https://en.wikipedia.org/wiki/Center_(group_theory))
``k`` that has degree four as a vector space over ``k``. In Morgenstern's explicit construction of Ramanujan graphs
for even characteristic ``q`` [morgenstern1994existence](@cite), a specific quaternion algebra is defined as

```math
\\begin{aligned}
\\mathscr{A} = k\\mathbf{1} + k\\mathbf{i} + k\\mathbf{j} + k\\mathbf{ij}
\\end{aligned}
```

with relations

```math
\\begin{aligned}
\\mathbf{i}^2 = \\mathbf{i} + \\varepsilon, \\mathbf{j}^2 = x, \\mathbf{ij} = \\mathbf{ji} + \\mathbf{j}
\\end{aligned}
```

The parameter ``\\varepsilon \\in \\mathbb{F}_q`` is chosen so that the polynomial ``f(x) = x^2 + x + \\varepsilon``
is irreducible over ``\\mathbb{F}_q``. This ensures ``\\mathcal{A}`` is a skewfield. The algebra is [ramified](https://en.wikipedia.org/wiki/Ramification_(mathematics))
at the finite place ``x`` and at ``1/x``.

The connection to graph theory arises from studying elements in the integral set

```math
\\begin{aligned}
\\mathscr{S} = \\mathbb{F}_q[x]\\mathbf{1} + \\mathbb{F}_q[x]\\mathbf{i} + \\mathbb{F}_q[x]\\mathbf{j} + \\mathbb{F}_q[x]\\mathbf{ij}
\\end{aligned}
```

The "basic norm x+1" elements are defined as

```math
\\begin{aligned}
\\xi = 1 + \\gamma\\mathbf{j} + \\delta\\mathbf{ij}, \\quad \\text{with} \\gamma, \\delta \\in \\mathbb{F}_q,
\\end{aligned}
```

satisfying the norm equation ``N(\\xi) = \\gamma^2 + \\gamma\\delta + \\delta^2\varepsilon = 1``. This equation
has exactly ``q+1`` solutions in ``\\mathbb{F}_q``, parameterizing generators ``\\xi_1, \\dots, \\xi_{q+1}``. These
generators define a [free product](https://en.wikipedia.org/wiki/Free_product) group 

```math
\\begin{aligned}
A(x) = \\langle \\xi_1 \\rangle * \\langle \\xi_2 \\rangle * \\cdots * \\langle \\xi_{q+1} \\rangle
\\end{aligned}
```

which acts simply transitively on the ``q+1``-regular tree ``T_{x+1} = G'_{x+1}/G'_{O_{x+1}}``. Taking the
quotient by a [congruence subgroup](https://en.wikipedia.org/wiki/Congruence_subgroup) ``A(g)``, where ``g(x)``
is irreducible of even degree, yields a finite ``(q+1)``-regular graph ``\\Gamma_g = A(g) \\backslash T_{x+1}``.
This graph is the Cayley graph of ``PSL_2(\\mathbb{F}_{q^d})`` with respect to the images of the ``q+1`` generators.

# Arguments
- `R`: Polynomial ring ``\\mathbb{F}_q[x]`` where `q` is a power of 2.
"""
function morgenstern_solutions(R::FqPolyRing)
    F = base_ring(R)
    unit = gen(F)
    q = order(F)
    f = morgenstern_f(R) # random sampler
    ε = coeff(f, 0)
    # (1, 0) is always a solution: 1² + 1·0 + 0²·ε = 1
    sols = [(one(F),zero(F))]
    # Theorem 5.13 of [morgenstern1994existence](@cite): find all solutions (γ, δ) ∈ 𝔽q² to γ² + γδ + δ²ε = 1
    for δ in F
        iszero(δ) && continue
        for γ in F
            if γ^2+γ*δ+δ^2*ε == one(F)
                push!(sols, (γ, δ))
                # in GF2, if (γ, δ) is solution, so is (γ+δ, δ)
                other_γ = γ+δ
                push!(sols, (other_γ, δ))
                break
            end
        end
    end
    @assert length(unique(sols)) == q+1
    return ε, unique(sols)
end

"""
The theorem *5.13* of [morgenstern1994existence](@cite) provides a method to construct families of
(q+1)-regular Ramanujan graphs for **even prime** powers q. This explicit construction produces Cayley
graphs of the *projective special linear group* ``\\mathrm{PSL}_2(\\mathbb{F}_{q^d})`` with respect to a
specific set of ``q+1`` generators. These generators are ``2 \\times 2`` matrices of the form:

```math
\\begin{aligned}
\\begin{pmatrix}
1 & \\gamma_k + \\delta_k\\mathbb{i} \\\\
(\\gamma_k + \\delta_k\\mathbb{i} + \\delta_k)x & 1
\\end{pmatrix}, \\quad k=1,\\ldots,q+1
\\end{aligned}
```

where q = 2^l is an *even* prime power, and d is an even integer extension degree. The field ``\\mathbb{F}_{q^d}`` is
constructed as ``\\mathbb{F}_q[x]/g(x)\\mathbb{F}_q[x]`` where g(x) is an irreducible polynomial of degree d. Within
this field, ``\\mathbb{i}`` denotes a root of the irreducible polynomial ``x^2 + x + \\varepsilon = 0``. The pairs
``(\\gamma_k, \\delta_k)`` are the q+1 solutions in `\\mathbb{F}_q^2`` to the ``\\gamma_k^2 + \\gamma_k\\delta_k + \\delta_k^2\\varepsilon = 1``.
And x is the polynomial variable that represents an element of ``\\mathbb{F}_{q^d}`` in the construction.

The same theorem states that the resulting Cayley graph ``\\Gamma_g`` has the following properties: it is
a (q+1)-regular Ramanujan graph of order ``|\\Gamma_g| = q^{3d} - q^d`` and is non-bipartite. The graph has
girth at least ``\\frac{2}{3}\\log_q|\\Gamma_g|`` and diameter at most ``2\\log_q|\\Gamma_g| + 2``. Furthermore,
as per Theorem *5.11*, all eigenvalues ``\\mu`` of the adjacency matrix satisfy ``|\\mu| \\leq 2\\sqrt{q}`` for ``\\mu \\neq ``\\pm(q+1)``.

!!! note
    In the construction of Morgenstern Ramanujan graphs for even prime powers ``q = 2^l``, we utilize
    the fact that in characteristic 2, the [projective special linear group](https://en.wikipedia.org/wiki/Projective_linear_group)
    ``\\mathrm{PSL}_2(\\mathbb{F}_q)`` is isomorphic to the [special linear group](https://en.wikipedia.org/wiki/Special_linear_group)
    ``\\mathrm{SL}_2(\\mathbb{F}_q)``. This isomorphism holds because the [center](https://en.wikipedia.org/wiki/Center_(group_theory))
    of ``\\mathrm{SL}_2(\\mathbb{F}_q)`` is trivial in characteristic 2 when q is even as confirmed by identity ``Z \\cap SL(2, \\mathbb{F}) = I``. 

# Arguments
- `l`: A positive integer specifying that q = 2^l, where q is the size of the base field ``\\mathbb{F}_q``.
- `i`: An *even* positive integer specifying the extension degree for the field ``\\mathbb{F}_{q^i}``. 
"""
function morgenstern_generators(l::Int, i::Int)
    @assert iseven(i) "Extension degree i must be even for Morgenstern construction (Theorem 5.13) of [morgenstern1994existence](@cite)"
    p = 2
    q = p^l
    qⁱ = q^i
    # finite fields 𝔽_q and extension 𝔽_qⁱ as in Section 5
    𝔽q, _ = finite_field(p, l)
    𝔽qⁱ, _ = finite_field(p, l * i)
    # Embedding morphism from 𝔽_q to 𝔽_qⁱ for field extension
    morph = embed(𝔽q, 𝔽qⁱ)
    # Given irreducible polynomial x² + x + ε, find solutions to γ² + γδ + δ²ε = 1
    R𝔽q, x = polynomial_ring(𝔽q, :x)
    ε, Bsols = morgenstern_solutions(R𝔽q)
    # Theorem 5.13: There are exactly q+1 solutions (γ,δ) ∈ 𝔽_q² to γ² + γδ + δ²ε = 1
    @assert length(Bsols) == q + 1 "Expected $((q+1)) solutions as per Theorem 5.13"
    R𝔽qⁱ, y = polynomial_ring(𝔽qⁱ, :y)
    𝕚s = roots(y^2+y+morph(ε))
    @assert length(𝕚s) == 2 "Irreducible quadratic x² + x + ε must have exactly 2 roots in extension field 𝔽_qⁱ"
    𝕚 = rand(𝕚s) # selecting one of the two roots at random
    # Note: In characteristic 2, PSL₂(𝔽_qⁱ) = SL₂(𝔽_qⁱ).
    SL₂qⁱ = SL(2, 𝔽qⁱ)
    @info "|SL₂(𝔽($(qⁱ)))| = $(order(SL₂qⁱ))"
    order(SL₂qⁱ) > 10_000 && @warn "We are working with a very big group, this will take a long time."
    order(SL₂qⁱ) > 300_000 && throw(ArgumentError("The group is too big, we refuse to even try to proceed."))
    B = eltype(SL₂qⁱ)[]
    xᵥₐₗ = gen(𝔽qⁱ)
    for (γ, δ) in Bsols
        # Embed γ, δ from 𝔽_q to 𝔽_qⁱ
        γ = morph(γ)
        δ = morph(δ)
        # see Eq. 21 of [morgenstern1994existence](@cite)
        mat = matrix(𝔽qⁱ, [1              γ+δ*𝕚; 
                           (γ+δ*𝕚+δ)*xᵥₐₗ     1])
        # normalize the matrix to have determinant 1 for SL₂
        detₘₐₜ = det(mat)
        @assert !iszero(detₘₐₜ) "Generator matrix must be invertible"
        normalizedₘₐₜ = mat*inv(sqrt(detₘₐₜ))
        @assert det(normalizedₘₐₜ) == one(𝔽qⁱ) "Normalized matrix must have determinant 1 for SL₂"
        g = SL₂qⁱ(normalizedₘₐₜ)
        @assert g^2 == one(SL₂qⁱ) "Each generator must have order 2 (Theorem 5.13)"
        push!(B, g)
    end
    return SL₂qⁱ, B
end

abstract type MorgensternAlgorithm end

struct AllPairs <: MorgensternAlgorithm end
struct FirstOnly <: MorgensternAlgorithm end

"""
Morgenstern showed that for every prime q, there exist infinitely many groups
``\\G_i = \\mathrm{PGL}_2(q^i)`` or ``G_i = \\mathrm{PSL}_2(q^i)``, each with a
symmetric generating set ``B_i`` of size ``q+1``, such that the Cayley graphs
``\\mathrm{Cay}(G_i, B_i)`` are Ramanujan graphs ([morgenstern1994existence](@cite),
[dinur2022locally](@cite)). That is, the second largest eigenvalue satisfies
``\\lambda(\\mathrm{Cay}(G_i, B_i)) \\leq 2\\sqrt{q}/(q+1)``.

The construction uses an explicit *arithmetic [lattice](https://en.wikipedia.org/wiki/Lattice_(discrete_subgroup))*
``\\Gamma`` in ``\\mathrm{PSL}_2(\\mathbb{F}_q)`` that is isomorphic to the [free product](https://en.wikipedia.org/wiki/Free_product)
of ``q+1`` copies of the cyclic group of order 2, where ``B = {b_0, b_1, \\ldots, b_q}``
consists of elements of order ``2`` [dinur2022locally](@cite). The Cayley graphs are
obtained as images ``B_i = \\phi(B)`` under [epimorphisms](https://en.wikipedia.org/wiki/Epimorphism)
``\\phi: \\Gamma \\to G_i``.

Dinur provided an *alternative* construction of symmetric generating sets ``A_i`` for
``G_i = \\mathrm{PGL}_2(q^i)`` such that the pairs ``(A_i, B_i)`` satisfy the total
non-congruency condition [dinur2022locally](@cite) for Morgenstern construction of
Ramanujan graphs using even prime q.

## Generator Constructions

### All Pairs Construction

The subgroup ``\\Lambda`` is generated by the symmetric set

```math
\\begin{aligned}
A = {b_t b_s \\mid b_t, b_s \\in B, t \\neq s}
\\end{aligned}
```

which has size ``k_1 = q^2 + q``. For a finite group ``G`` with symmetric generating set
``B`` where each element has order 2, the all pairs construction generates the alternative generating set:

```math
\\begin{aligned}
A = {b_i b_j \\mid b_i, b_j \\in B, i \\neq j}
\\end{aligned}
```

This construction yields Cayley graphs ``\\text{Cay}(G_i, A_i)`` with spectral expansion satisfying:

```math
\\begin{aligned}
\\lambda(\\text{Cay}(G_i, A_i)) < \\frac{3q-1}{q^2+q} < \\frac{3\\sqrt{k_1-1}}{k_1}
\\end{aligned}
```
"""
function alternative_morgenstern_generators(B::AbstractVector, ::AllPairs)
    @assert !isempty(B)
    @assert all(x -> x isa GroupElem, B) "All elements in B must be group elements"
    A = eltype(B)[]
    N = length(B)
    sizehint!(A, N*(N-1))
    @assert N >= 2 "Generating set B must have at least 2 elements"
    for i in 1:N
        for j in 1:N
            if i!=j
                a = B[i]*B[j]
                push!(A,a)
            end
        end
    end
    @assert length(A) == N*(N-1) "Output set should have size |B|*(|B|-1) = $(N*(N-1))"
    @assert all(a -> inv(a) in A, A) "Set A must be symmetric"
    @assert !(one(B[1]) in A) "Identity element cannot be in generating set A"
    @assert (N-1 <= 4) || (length(unique(A)) == length(A)) "Claim 6.1(iii): All generators distinct for q > 4 [dinur2022locally](@cite)"
    @assert all(a -> a*a != one(B[1]), A) "Claim 6.1(iii): Each element must have order > 2 [dinur2022locally](@cite)"
    @assert all(a -> !(a in B), A) "Claim 6.1(iii): Generated elements cannot be in original set B [dinur2022locally](@cite)"
    return A
end

"""
## Generator Constructions

### First Element Construction

A more efficient construction recognizes that ``\\Lambda`` is actually a
[free group](https://en.wikipedia.org/wiki/Free_group) on the q generators
``\\b_0 b_j : j = 1, \\ldots, q``. Since ``(b_0 b_j)^{-1} = b_j b_0``, we obtain
the symmetric generating set:

```math
\\begin{aligned}
A' = {b_0 b_j, b_j b_0 \\mid j = 1, \\ldots, q}
\\end{aligned}
```

This construction produces Cayley graphs with improved spectral properties:

```math
\\begin{aligned}
\\lambda(\\text{Cay}(G_i, A'_i)) < \\frac{3\\sqrt{2q-1}}{2q}
\\end{aligned}
```

The default constructor uses the more efficient `FirstOnly` algorithm, providing better
spectral expansion with smaller generating sets of size ``k_1 = 2q`` compared to ``k_1 = q^2 + q``
for the all pairs construction.
"""
function alternative_morgenstern_generators(B::AbstractVector, ::FirstOnly)
    @assert !isempty(B) "Generating set B must not be empty"
    @assert all(x -> x isa GroupElem, B) "All elements in B must be group elements"
    A = eltype(B)[]
    N = length(B)
    sizehint!(A, 2*(N-1))
    @assert N >= 2 "Generating set B must have at least 2 elements"
    for i in 2:N
        push!(A,B[1]*B[i])
        push!(A,B[i]*B[1])
    end
    @assert length(A) == 2*(N-1) "Output set should have size 2*(|B|-1) = $(2*(N-1))"
    @assert all(a -> inv(a) in A, A) "Set A must be symmetric"
    @assert !(one(B[1]) in A) "Identity element cannot be in generating set A"
    @assert (N - 1 <= 4) || (length(unique(A)) == length(A)) "All generators distinct for q > 4"
    @assert all(a -> a*a != one(B[1]), A) "Claim 6.1(iii): Each element must have order > 2"
    @assert all(a -> !(a in B), A) "Generated elements cannot be in original set B"
    return A
end

alternative_morgenstern_generators(B) = alternative_morgenstern_generators(B, FirstOnly())
