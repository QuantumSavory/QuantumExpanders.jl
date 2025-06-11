using Nemo
using Oscar
using LinearAlgebra
using Random

"""
Randomly sample for ε∈𝔽q until you find one such that x^2+x+ε∈𝔽q[x] is irreducible.

Takes as input 𝔽q[x], the polynomial ring we want to work with.

Returns x^2+x+ε.

See [morgenstern1994existence](@cite).
"""
function morgenstern_f(R)
    x = gen(R)
    count = 0
    while true
        count += 1
        p = x^2+x+rand(R,0:0)
        if is_irreducible(p)
            @debug "`morgenstern_f` ran $(count) attempt(s)"
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
function morgenstern_solutions(R)
    F = base_ring(R)
    unit = gen(F)
    q = size(F)
    f = morgenstern_f(R) # random sampler
    ε = coefficients(f)[0]
    sols = [(one(F),zero(F))]
    for s in F
        fs = f(s)
        sfs = sqrt(fs)
        α = s*inv(sfs)
        @assert s*inv(sfs) == inv(sfs)*s
        β = inv(sfs)
        @assert α^2 + α*β + ε*β^2 == 1
        push!(sols, (α,β))
    end
    @assert length(unique(sols))==q+1
    return ε, sols
end

"""
Give all Morgenstern generators over PSL₂qⁱ, where i is even, q=2ˡ, and p is prime.

Returns (SL₂qⁱ, B), where B is the list of generators. As PSL₂qⁱ=SL₂qⁱ, we keep in SL₂qⁱ.

See [morgenstern1994existence](@cite).
"""
function morgenstern_generators(l,i)
    @assert iseven(i)
    p = 2
    q = p^l
    qⁱ = q^i
    @info "q = 2^$(l) = $(q)"
    @info "qⁱ = $(q)^$(i) = $(qⁱ)"
    𝔽q , unit = FiniteField(p,l)
    𝔽qⁱ, punit = FiniteField(p,l*i)
    morph = embed(𝔽q,𝔽qⁱ)
    R𝔽q, x = PolynomialRing(𝔽q, "x")
    R𝔽qⁱ, y = PolynomialRing(𝔽qⁱ, "y")
    ε, Bsols = morgenstern_solutions(R𝔽q)
    @assert length(Bsols) == q+1
    @info "|B| = q+1 = $(length(Bsols))"
    @info "ε = $(ε)"
    𝕚s = roots(y^2+y+morph(ε))
    𝕚 = rand(𝕚s) # selecting one of the two roots at random
    @info "𝕚 = $(𝕚)"
    # PSL₂qⁱ and PGL₂qⁱ are the same, so we are not going to try to work with the larger GL₂qⁱ
    #GL₂qⁱ = general_linear_group(2,𝔽qⁱ)
    #@info "|GL₂(𝔽(qⁱ))| = $(length(GL₂qⁱ))"
    SL₂qⁱ = special_linear_group(2,𝔽qⁱ)
    @info "|SL₂(𝔽(qⁱ))| = $(length(SL₂qⁱ))"
    if length(SL₂qⁱ)>10_000
        @warn "We are working with a very big group, this will take a long time."
    end
    if length(SL₂qⁱ)>300_000
        error("The group is too big, we refuse to even try to proceed.")
    end
    # The Center is a single element when p=2, so PSL and SL are the same,
    # therefore the computations below are not necessary. VERIFY
    #CSL₂qⁱ, Cₘₒᵣₚₕ = center(SL₂qⁱ) # seems to take time that scales with the size of SL₂qⁱ even though it is either 1 or 2 element group.
    #@info "|Center of SL₂(𝔽(qⁱ))| = $(length(CSL₂qⁱ))"
    #PSL₂qⁱ, Pₘₒᵣₚₕ = quo(SL₂qⁱ,CSL₂qⁱ)
    #@info "|PSL₂(𝔽(qⁱ))| = $(length(PSL₂qⁱ))"
    #@assert length(GL₂qⁱ) == length(SL₂qⁱ) == length(PSL₂qⁱ)

    slunit = one(SL₂qⁱ)
    B = typeof(slunit)[]
    for sol in Bsols
        γ,δ = morph.(sol)
        #γ+δ*𝕚 ∈ 𝔽qⁱ
        #(γ+δ*𝕚+δ)*morph(unit) ∈ 𝔽qⁱ
        _mat = 𝔽qⁱ[1 γ+δ*𝕚; (γ+δ*𝕚+δ)*punit 1]
        _matp = _mat * inv(sqrt(det(_mat))) # XXX This seems implicit in the papers, VERIFY
        @assert _mat * inv(sqrt(det(_mat))) == inv(sqrt(det(_mat))) * _mat
        b = SL₂qⁱ(_matp)
        @assert b^2==slunit
        push!(B,b)
    end
    SL₂qⁱ, B
end

abstract type MorgensternAlgorithm end

struct AllPairs <: MorgensternAlgorithm end
struct FirstOnly <: MorgensternAlgorithm end

"""
    alternative_morgenstern_generators(B::AbstractVector, ::AllPairs)

Create alternative Morgenstern generators using all pairwise products (i≠j).

Introduced in Sec. 6.1 of [dinur2022locally](@cite).
Building upon [morgenstern1994existence](@cite).
"""
function alternative_morgenstern_generators(B::AbstractVector, ::AllPairs)
    A = eltype(B)[]
    N = length(B)
    for i in 1:N
        for j in 1:N
            if i ≠ j
                a = B[i] * B[j]
                push!(A, a)
            end
        end
    end
    return A
end

"""
    alternative_morgenstern_generators(B::AbstractVector, ::FirstOnly)

Create alternative Morgenstern generators using products with first element only.

Introduced as the "better" alternative in Sec. 6.1 of [dinur2022locally](@cite).
Building upon [morgenstern1994existence](@cite).
"""
function alternative_morgenstern_generators(B::AbstractVector, ::FirstOnly)
    A = eltype(B)[]
    N = length(B)
    for i in 2:N
        push!(A, B[1] * B[i])
        push!(A, B[i] * B[1])
    end
    return A
end

# Convenience methods
alternative_morgenstern_generators(B) = alternative_morgenstern_generators(B, FirstOnly())
