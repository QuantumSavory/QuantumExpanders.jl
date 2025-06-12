"""
A module for constructing quantum LDPC codes
"""

module QuantumExpanders

using Nemo
using Oscar
using LinearAlgebra
using Random
using Graphs
using Oscar:coefficients
using Graphs:add_edge!, nv, ne, neighbors, Graphs, edges, Edge, src, dst, degree, adjacency_matrix

include("cayley_graphs.jl")
include("morgenstern.jl")
include("ramanujan.jl")
include("tensor_codes.jl")

export gen_code, gen_good_code, tanner_code, tanner_code_quadripartite,
cayley_complex_square_graphs, cayley_complex_square_graphs_quadripartite,
morgenstern_generators, alternative_morgenstern_generators, scalar_matrices_GL,
scalar_matrices_SL, solve_four_squares, process_solutions, create_generators,
construct_cayley_graph, ramanujan_graph, is_ramanujan, legendre_symbol,
edge_vertex_incidence_graph, is_unbalanced_bipartite, alon_chung_lemma,
expander_code_parity_matrix

##

Random.seed!(42)
function gen_code(ra,rb,bipartite=true,generators="A,B",use_same_local_code=false)
    l=1
    i=2
    @time SL₂qⁱ, B = morgenstern_generators(l,i)
    @time A = alternative_morgenstern_generators(B)
    @show length(SL₂qⁱ), length(A), length(B)
    @show length(SL₂qⁱ)*length(A)*length(B)
    if bipartite
        @assert is_nonconjugate(SL₂qⁱ, A, B)
        @assert is_symmetric_gen(A)
        @assert is_symmetric_gen(B)
    end

    ##
    if !bipartite && generators=="B"
        A = B 
    end
    
    if !bipartite && generators=="A"
        B = A
    end

    𝒢₀□, 𝒢₁□, edge₀_q_idx, edge₁_q_idx, edge₀_ab_idx, edge₁_ab_idx = bipartite ? cayley_complex_square_graphs(SL₂qⁱ,A,B) : cayley_complex_square_graphs_quadripartite(SL₂qⁱ,A, B) 
    # Be careful with notation here as we are interchangeably using
    # parity check matrices and generator matrices
    # while also using codes and their dual codes.
    # This can lead to confusion as
    # the parity check matrix of a code is the generator matrix of its dual code.
    Hᴬ = uniformly_random_code_checkmatrix(ra,length(A))
    Hᴮ = uniformly_random_code_checkmatrix(rb,length(B)) # TODO - make use of use_same_local_code
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

    𝒞ᶻ = bipartite ? tanner_code(𝒢₀□,edge₀_q_idx,edge₀_ab_idx,C₀) : tanner_code_quadripartite(𝒢₀□,edge₀_q_idx,edge₀_ab_idx,C₀)
    𝒞ˣ = bipartite ? tanner_code(𝒢₁□,edge₁_q_idx,edge₁_ab_idx,C₁) : tanner_code_quadripartite(𝒢₁□,edge₁_q_idx,edge₁_ab_idx,C₁)  
    @show r1 = rank(𝒞ᶻ)
    @show r2 = rank(𝒞ˣ)
    #@assert good_css(dual_code(𝒞ˣ),dual_code(𝒞ᶻ))
    @assert good_css(𝒞ˣ,𝒞ᶻ)
    𝒞ˣ,𝒞ᶻ
end

##

function gen_good_code(ra,rb,minweighta=1,minweightb=1,bipartite=true,generators="A,B",use_same_local_code=false,number_iterations=100)
    for i in 1:number_iterations
        𝒞ˣ,𝒞ᶻ=gen_code(ra,rb,bipartite,generators,use_same_local_code)
        zm = minimum(unique(sum(𝒞ᶻ, dims=1)))
        xm = minimum(unique(sum(𝒞ˣ, dims=1)))
        if zm>=minweighta && xm>=minweightb
            return 𝒞ˣ,𝒞ᶻ
        end
    end
end

end #module