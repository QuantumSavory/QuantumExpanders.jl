# This file contains a "scratch pad" of various tests and demos related to generating
# expander graphs based on the Morgenstern construction.
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
# Making a Cayley graph for even q using Morgenstern's construction
##

l=2
i=2
@time SL₂qⁱ, B = morgenstern_generators(l,i)
@time graphB = cayley_right(SL₂qⁱ, B)
@time evalsB = adjacency_spectrum(graphB) # slow, dense, there should be a better way to do it, especially if we care about only two eigvals

# Checks from [morgenstern1994existence](@cite). TODO
q = 2^l
d = q+1
N = size(graphB,1)
@assert evalsB[end-1] <= 2 * sqrt(q) # is Ramanujan
@assert evalsB[end] ≈ d # q+1 regular
@assert N == q^(3i)-q^i
# @assert is not bipartite
# @assert girth >= 2/3 * log(q,N)
@assert diameter(graphB) <= 2*log(q,N)+2
# @assert chromatic number >= (q+1) / (2*sqrt(q)) + 1
# @assert independence number <= 2*N*sqrt(q) / (q+1)

# Make the second set of generators for the Cayley complex
@time Al = alternative_long_morgenstern_generators(B)
@assert is_nonconjugate(SL₂qⁱ, Al, B)
@assert length(unique(Al)) == length(Al)
@assert all(Al.^2 .!= (one(Al[1]),))

# Make the second set of generators for the Cayley complex
@time A = alternative_morgenstern_generators(B)
@assert is_nonconjugate(SL₂qⁱ, A, B)
@assert length(unique(A)) == length(A)
@assert all(A.^2 .!= (one(A[1]),))

##
# Compare spectral properties of Morgenstern vs Random graphs
##

l=2 # 1 or 2
i=2 # 2 or 4
q = 2^l
d = q+1
dₐₗ = q^2+q
dₐ = 2q
N = q^(3i)-q^i

second_eval(graph) = adjacency_spectrum(graph)[end-1]
morg_graph(l,i) = cayley_right(morgenstern_generators(l,i)...)
function morg_alt_graph(l,i)
    SL₂qⁱ, B = morgenstern_generators(l,i)
    A = alternative_morgenstern_generators(B)
    cayley_left(SL₂qⁱ,A)
end
function morg_alt_long_graph(l,i)
    SL₂qⁱ, B = morgenstern_generators(l,i)
    A = alternative_long_morgenstern_generators(B)
    cayley_left(SL₂qⁱ,A)
end

##

samples = 3
morg_spec = [second_eval(morg_graph(l,i)) for _ in 1:samples]
morg_alt_spec = [second_eval(morg_alt_graph(l,i)) for _ in 1:samples]
morg_alt_long_spec = [second_eval(morg_alt_long_graph(l,i)) for _ in 1:samples]
rand_spec = [second_eval(random_regular_graph(N,d)) for _ in 1:samples]
rand_alt_spec = [second_eval(random_regular_graph(N,dₐ)) for _ in 1:samples]
rand_alt_long_spec = [second_eval(random_regular_graph(N,dₐₗ)) for _ in 1:samples]

##

fig = Figure(resolution=(1000,300))
ax1 = Axis(fig[1,1],
    title="Morgenstern vs Random graphs\n N=$(N) deg=$(d)", ylabel="λ₂",
    xticksvisible = false,
    xticks=([1,2],["Morgenstern", "Random"]))
ax2 = Axis(fig[1,2],
    title="Alternative\nMorgenstern vs Random graphs\n N=$(N) deg=$(dₐ)", ylabel="λ₂",
    xticksvisible = false,
    xticks=([1,2],["Morgenstern", "Random"]))
ax3 = Axis(fig[1,3],
    title="Alternative long\nMorgenstern vs Random graphs\n N=$(N) deg=$(dₐₗ)", ylabel="λ₂",
    xticksvisible = false,
    xticks=([1,2],["Morgenstern", "Random"]))
ramabound = 2 * sqrt(d-1)
ramaboundₐ = 2 * sqrt(dₐ-1)
ramaboundₐₗ = 2 * sqrt(dₐₗ-1)
scatter!(ax1, one.(morg_spec),morg_spec, label="Morgenstern")
scatter!(ax1, one.(rand_spec)*2,rand_spec, label="Random")
hdeg = hlines!(ax1, [d], label="degree", color=:black)
hram = hlines!(ax1, [ramabound], label="Ramanujan gap", color=:black, linestyle=:dash)
ylims!(ax1, ramabound/2, d*1.05)
xlims!(ax1, 0.5, 2.5)
scatter!(ax2, one.(morg_alt_spec),morg_alt_spec, label="Morgenstern")
scatter!(ax2, one.(rand_alt_spec)*2,rand_alt_spec, label="Random")
hdeg = hlines!(ax2, [dₐ], label="degree", color=:black)
hram = hlines!(ax2, [ramaboundₐ], label="Ramanujan gap", color=:black, linestyle=:dash)
ylims!(ax2, ramaboundₐ/2, dₐ*1.05)
xlims!(ax2, 0.5, 2.5)
scatter!(ax3, one.(morg_alt_long_spec),morg_alt_long_spec, label="Morgenstern")
scatter!(ax3, one.(rand_alt_long_spec)*2,rand_alt_long_spec, label="Random")
hdeg = hlines!(ax3, [dₐₗ], label="degree", color=:black)
hram = hlines!(ax3, [ramaboundₐₗ], label="Ramanujan gap", color=:black, linestyle=:dash)
ylims!(ax3, ramaboundₐₗ/2, dₐₗ*1.05)
xlims!(ax3, 0.5, 2.5)
Legend(fig[1,4], [hdeg,hram], ["degree", "Ramanujan gap"])
display(fig)

##
# To plot the graph... not very useful
##

using GraphMakie    # To plot graphs
using NetworkLayout # To do spectral layout

graphplot(graphB, layout=Spectral(dim=2)) # slow and not really useful
