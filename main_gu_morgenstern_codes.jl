# This file contains a "scratch pad" of various tests and demos related to generating
# Tanner codes based as described in Gu et al, but using the Morgenstern graphs.
#
# We do not have |A|=|B| so we are not exactly following the original prescription.
# In particular, we have taken a lot of liberties with the generation of the lower tensor codes.
#
# Not particularly careful with what is getting imported where, would clean up when
# turning it into a package.

using Nemo
using Oscar
using LinearAlgebra
using Random
using Graphs

using CairoMakie # Plotting libraries that are extremely slow to load due to
                 # lattency limitations in the Julia compiler. Otherwise great and fast.

include("morgenstern.jl")
include("cayley_graphs.jl")
include("tensor_codes.jl")

##

Random.seed!(42)
l=1
i=2
@time SL₂qⁱ, B = morgenstern_generators(l,i)
@time A = alternative_morgenstern_generators(B)
@show length(SL₂qⁱ), length(A), length(B)
@show length(SL₂qⁱ)*length(A)*length(B)
@assert is_nonconjugate(SL₂qⁱ, A, B)
@assert is_symmetric_gen(A)
@assert is_symmetric_gen(B)

##

𝒢₀□, 𝒢₁□, edge₀_q_idx, edge₁_q_idx, edge₀_ab_idx, edge₁_ab_idx = cayley_complex_square_graphs(SL₂qⁱ,A,B)
# Be careful with notation here as we are interchangeably using
# parity check matrices and generator matrices
# while also using codes and their dual codes.
# This can lead to confusion as
# the parity check matrix of a code is the generator matrix of its dual code.
Random.seed!(42)
Hᴬ = uniformly_random_code_checkmatrix(2,length(A))
Hᴮ = uniformly_random_code_checkmatrix(2,length(B))
Cᴬ = dual_code(Hᴬ)
Cᴮ = dual_code(Hᴮ)
C₀ = kron(Cᴬ,Cᴮ) # consider it as a generator matrix
@show size(C₀)
C₀⁺ = dual_code(C₀)
C₁ = kron(Hᴬ,Hᴮ) # consider it as a generator matrix
@show size(C₁)
C₁⁺ = dual_code(C₁)
@assert good_css(Hᴬ,Cᴬ)
@assert good_css(Hᴮ,Cᴮ)
@assert good_css(C₀,C₁)

𝒞ᶻ = tanner_code(𝒢₀□,edge₀_q_idx,edge₀_ab_idx,C₀)
𝒞ˣ = tanner_code(𝒢₁□,edge₁_q_idx,edge₁_ab_idx,C₁)
r1 = rank(𝒞ᶻ)
r2 = rank(𝒞ˣ)
@assert good_css(dual_code(𝒞ˣ),dual_code(𝒞ᶻ))
@assert good_css(𝒞ˣ,𝒞ᶻ)

##

using QuantumClifford
using QuantumCliffordPlots

stab = css(𝒞ˣ,𝒞ᶻ)
@assert good_css(stab)
for row in stab
    @assert all(==(0), QuantumClifford.comm(stab[1],stab))
end

QuantumClifford.stab_looks_good(stab) # internal function used for sanity checks
stabilizerplot(stab)
