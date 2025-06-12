using Nemo
using Oscar
using Oscar: embed
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
    q = order(F)
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
    morgenstern_generators(l::Int, i::Int)

Compute generators for the Morgenstern construction of expander graphs over PSL₂(𝔽_{qⁱ}), where
`q = 2ˡ` and `i` is even. Returns:
- `generators`: List of `q+1` matrices in PSL₂(𝔽_{qⁱ}) (note PSL₂=SL₂ in characteristic 2)
- `group_order`: The size of PSL₂(𝔽_{qⁱ}) = qⁱ(qⁱ²-1)

The `morgenstern_generators` is optimized for Cayley graph construction via BFS traversal
without full group enumeration. For graph construction, use `cayley_left(generators, group_order)`
or `cayley_right(generators, group_order)`.

See [morgenstern1994existence](@cite).
"""
function morgenstern_generators(l, i)
    @assert iseven(i)
    p = 2
    q = p^l
    qⁱ = q^i
    𝔽q = GF(q)
    𝔽qⁱ = GF(qⁱ)
    generators = Vector{Matrix{Nemo.FqFieldElem}}(undef, q+1)
    # Find irreducible polynomial x² + x + ε and solutions
    R, x = polynomial_ring(𝔽q, "x")
    ε, Bsols = morgenstern_solutions(R)
    # Find root i of x² + x + ε in 𝔽qⁱ
    Ry, y = polynomial_ring(𝔽qⁱ, "y")
    f = y^2 + y + 𝔽qⁱ(ε)
    𝕚 = roots(f)[1]
    # Generate the q+1 generators
    for (idx, (γ, δ)) in enumerate(Bsols)
        γⁱ = 𝔽qⁱ(γ)
        δⁱ = 𝔽qⁱ(δ)
        a = γⁱ + δⁱ*𝕚
        b = (γⁱ + δⁱ*𝕚 + δⁱ)*gen(𝔽qⁱ)
        mat = [one(𝔽qⁱ) a; b one(𝔽qⁱ)]
        det_mat = one(𝔽qⁱ) - a*b
        sqrt_det = sqrt(det_mat)
        mat_normalized = mat ./ sqrt_det
        generators[idx] = mat_normalized
    end
    group_order = qⁱ * (qⁱ^2 - 1)
    return generators, group_order
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
