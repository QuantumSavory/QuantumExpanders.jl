# This file contains a "scratch pad" of various tests and demos.
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

l=2
i=2
@time SL₂qⁱ, B = morgenstern_generators(l,i)
@time graphB = cayley_right(SL₂qⁱ, B)
@time evalsB = adjacency_spectrum(graphB) # slow, dense, there should be a better way to do it, especially if we care about only two eigvals
@time A = alternative_morgenstern_generators(B)

@time @assert is_self_nonconjugate(SL₂qⁱ, B)

# Checks from [morgenstern1994existence](@cite). TODO
q = 2^l
N = size(graphB,1)
@assert evals[end-1] <= 2 * sqrt(q) # is Ramanujan
# @assert is q+1 regular
N == q^(3i)-q^i
# @assert is not bipartite
# @assert girth >= 2/3 * log(q,N)
@assert diameter(graphB) <= 2*log(q,N)+2
# @assert chromatic number >= (q+1) / (2*sqrt(q)) + 1
# @assert independence number <= 2*N*sqrt(q) / (q+1)

##

rgraph = random_regular_graph(length(SL₂qⁱ),length(B))
revals = adjacency_spectrum(rgraph)
@assert diameter(rgraph) <= 2*log(q,N)+2

##
# To plot the graph... not very useful

using GraphMakie    # To plot graphs
using NetworkLayout # To do spectral layout

graphplot(graphB, layout=Spectral(dim=2)) # slow and not really useful

##

p = 5
l = 3
q = p^l
𝔽q , unit = FiniteField(p,l)
@time SL₂q = special_linear_group(2,𝔽q);
length(SL₂q)
@time CSL₂q, Cₘₒᵣₚₕ = Oscar.center(SL₂q); # GETTING SLOWER AND SLOWER
length(CSL₂q)
@time PSL₂q, Pₘₒᵣₚₕ = quo(SL₂q,CSL₂q);
length(PSL₂q)
#@time collect(PSL₂q);
