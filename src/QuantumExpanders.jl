"""
A module for constructing quantum LDPC codes
"""

module QuantumExpanders

using Nemo
using Oscar
using LinearAlgebra
using Random
using Graphs
using Graphs: add_edge!, nv, ne, neighbors, Graphs, edges, Edge, src, dst, degree, adjacency_matrix, add_vertex!, has_edge,
vertices
using Oscar: coefficients, zzModRingElem, MatSpace, zero_matrix, base_ring, lift, matrix_space, zzModMatrix,
residue_ring, ZZ, nullspace, transpose, base_ring, kron, Matrix, embed, GroupElem, MatrixGroup, FqField
using Graphs: add_edge!, nv, ne, neighbors, Graphs
using Multigraphs
using ProgressMeter

include("cayley_graphs.jl")
include("quantum_tanner_code.jl")
include("tensor_codes.jl")
include("morgenstern.jl")
include("lubotzky_phillips_sarnak_ramanujan.jl")

export
    gen_code, gen_good_code, tanner_code, tanner_code_quadripartite,
    # Cayley graphs
    cayley_right, cayley_left, is_nonconjugate, is_ramanujan, is_symmetric_gen,
    cayley_complex_square_graphs, cayley_complex_square_graphs_quadripartite,
    # Morgenstern Ramanujan
    morgenstern_generators, alternative_morgenstern_generators, AllPairs, FirstOnly,
    # Lubotzky-Phillips-Sarnak Ramanujan
    ramanujan_graph, scalar_matrices_GL, scalar_matrices_SL, solve_four_squares,
    process_solutions, lps_generators, legendre_symbol, LPS, lps_graph,
    edge_vertex_incidence_graph, is_unbalanced_bipartite, alon_chung_lemma

end #module
