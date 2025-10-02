"""
Randomly sample for ε∈𝔽q until you find one such that x^2+x+ε∈𝔽q[x] is irreducible.

Takes as input 𝔽q[x], the polynomial ring we want to work with.

Returns x^2+x+ε.

See [morgenstern1994existence](@cite).
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
Internally call the `morgenstern_f` sampler to find an irreducible x^2+x+ε∈𝔽q[x].
Then find all q+1 solutions (γ,δ) of γ²+γδ+δ²ε=1.

Takes as input 𝔽q[x], the polynomial ring we want to work with.

Returns (ε, sols), where sols is the list of solutions.

See [morgenstern1994existence](@cite).
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
Create alternative Morgenstern generators using all pairwise products (i≠j).

Introduced in Sec. 6.1 of [dinur2022locally](@cite).
Building upon [morgenstern1994existence](@cite).
"""
function alternative_morgenstern_generators(B::AbstractVector, ::AllPairs)
    A = eltype(B)[]
    N = length(B)
    for i in 1:N
        for j in 1:N
            if i!=j
                a = B[i]*B[j]
                push!(A,a)
            end
        end
    end
    return A
end

"""
Create alternative Morgenstern generators using products with first element only.

Introduced as the "better" alternative in Sec. 6.1 of [dinur2022locally](@cite).
Building upon [morgenstern1994existence](@cite).
"""
function alternative_morgenstern_generators(B::AbstractVector, ::FirstOnly)
    A = eltype(B)[]
    N = length(B)
    for i in 2:N
        push!(A,B[1]*B[i])
        push!(A,B[i]*B[1])
    end
    return A
end

alternative_morgenstern_generators(B) = alternative_morgenstern_generators(B, FirstOnly())
