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
using Graphs: add_edge!, nv, ne, neighbors, Graphs, edges, Edge, src, dst, degree, adjacency_matrix, add_vertex!, has_edge

include("cayley_graphs.jl")
include("quantum_tanner_code.jl")
include("tensor_codes.jl")
include("morgenstern.jl")
include("ramanujan.jl")

export
    gen_code, gen_good_code, tanner_code, tanner_code_quadripartite,
    # Cayley graphs
    cayley_right, cayley_left, is_nonconjugate, is_ramanujan, is_symmetric_gen,
    cayley_complex_square_graphs, cayley_complex_square_graphs_quadripartite,
    # morgenstern Ramanujan
    morgenstern_generators, alternative_morgenstern_generators, AllPairs, FirstOnly,
    # Lubotzky-Phillips-Sarnak Ramanujan
    ramanujan_graph, scalar_matrices_GL, scalar_matrices_SL, solve_four_squares,
    process_solutions, create_generators, construct_cayley_graph, legendre_symbol,
    edge_vertex_incidence_graph, is_unbalanced_bipartite, alon_chung_lemma

end #module
