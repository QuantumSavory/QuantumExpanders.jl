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
using Graphs:add_edge!,nv,ne,neighbors,Graphs

include("morgenstern.jl")
include("cayley_graphs.jl")
include("quantum_tanner_code.jl")
include("tensor_codes.jl")

export gen_code, gen_good_code, tanner_code, tanner_code_quadripartite, cayley_complex_square_graphs, cayley_complex_square_graphs_quadripartite, morgenstern_generators, alternative_morgenstern_generators, AllPairs, FirstOnly

end #module