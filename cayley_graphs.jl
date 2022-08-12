using Nemo
using Oscar
using LinearAlgebra
using Random
using Graphs

"""Construct the Cayleyʳⁱᵍʰᵗ graph for a given group and set of generators."""
function cayley_right(group,generators)
    idx_to_mat = collect(group); # TODO see if there is a better (lazy?) way to enumerate
    mat_to_idx = Dict(mat=>i for (i,mat) in pairs(idx_to_mat))

    N = length(group)
    graph = SimpleGraph(N)
    for (i,g) in pairs(idx_to_mat)
        for b in generators
            j = mat_to_idx[g*b]
            add_edge!(graph,i,j)
        end
    end
    graph
end

"""Construct the Cayleyˡᵉᶠᵗ graph for a given group and set of generators."""
function cayley_left(group,generators)
    idx_to_mat = collect(group); # TODO see if there is a better (lazy?) way to enumerate
    mat_to_idx = Dict(mat=>i for (i,mat) in pairs(idx_to_mat))

    N = length(group)
    graph = SimpleGraph(N)
    for (i,g) in pairs(idx_to_mat)
        for b in generators
            j = mat_to_idx[b*g]
            add_edge!(graph,i,j)
        end
    end
    graph
end

"""Check the TNC condition of [dinur2022locally](@cite)."""
function is_nonconjugate(group,genA,genB)
    genset = Set(genB)
    for g in group
        for b in genA
            if inv(g)*b*g ∈ genset
                return false
            end
        end
    end
    true
end
