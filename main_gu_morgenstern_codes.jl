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

𝒢₀□, 𝒢₁□, edge₀_index, edge₁_index = cayley_complex_square_graphs(SL₂qⁱ,A,B)
Hᴬ = uniformly_random_code(0.9,length(A))
Hᴮ = uniformly_random_code(0.9,length(B))
C₀ = kron(Hᴬ,Hᴮ)
C₀⁺ = dual_code(C₀)
r,Δ² = size(C₀⁺)
@show r/Δ²
C₁ = kron(dual_code(Hᴬ),dual_code(Hᴮ))
C₁⁺ = dual_code(C₁)

𝒞ᶻ = tanner_code(𝒢₀□,edge₀_index,C₀⁺)
𝒞ˣ = tanner_code(𝒢₁□,edge₁_index,C₁⁺)
