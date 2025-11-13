"""
A module for constructing quantum LDPC codes
"""

module QuantumExpanders

using Nemo
using Oscar
using QECCore
import QECCore: code_n, code_k, parity_matrix, parity_matrix_z, parity_matrix_x, distance, AbstractCSSCode
using QuantumClifford
using QuantumClifford: Stabilizer, comm
using LinearAlgebra
using Random
using Graphs
using Graphs: add_edge!, nv, ne, neighbors, Graphs, edges, Edge, src, dst, degree, adjacency_matrix, add_vertex!, has_edge,
vertices, induced_subgraph, AbstractGraph, is_bipartite, bipartite_map, has_edge
using Oscar: coefficients, zzModRingElem, MatSpace, zero_matrix, base_ring, lift, matrix_space, zzModMatrix,
residue_ring, ZZ, nullspace, transpose, base_ring, kron, Matrix, embed, GroupElem, MatrixGroup, FqField, transpose, matrix,
rank, GF, FPGroup, FPGroupElem, Group, GroupElem
using Graphs: add_edge!, nv, ne, neighbors, Graphs
using Multigraphs
using ProgressMeter

include("cayley_graphs.jl")
include("tensor_codes.jl")
include("morgenstern.jl")
include("quantum_tanner_codes.jl")
include("quantum_tanner_code_multigraphs.jl")
include("lubotzky_phillips_sarnak_ramanujan.jl")

export
    gen_code, gen_good_code, tanner_code, tanner_code_quadripartite,
    # Cayley graphs
    cayley_right, cayley_left, is_nonconjugate, is_ramanujan, is_symmetric_gen,
    cayley_complex_square_graphs, cayley_complex_square_graphs_quadripartite,
    # Morgenstern Ramanujan
    morgenstern_solutions, morgenstern_generators, alternative_morgenstern_generators,
    AllPairs, FirstOnly,
    # Lubotzky-Phillips-Sarnak Ramanujan
    scalar_matrices_GL, scalar_matrices_SL, solve_four_squares, process_solutions,
    lps_generators, legendre_symbol, lps_graph, LPS, is_ramanujan,
    # Tensor codes
    uniformly_random_code_checkmatrix, dual_code,
    # Quantum Tanner codes
    enumerate_squares, random_code_pair, convert_squares_to_incidence_matrix, QuantumTannerCode,
    parity_matrix, parity_matrix_x, parity_matrix_z, parity_matrix_xz, code_n, code_k,
    # tensor codes
    uniformly_random_code_checkmatrix, dual_code, good_css

end #module
