using Nemo
using Oscar
using LinearAlgebra
using Random
using Graphs
using Multigraphs
using ProgressMeter

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

"""Construct the Cayley complex quare graphs 𝒢₀□ and 𝒢₁□ as presented in [gu2022efficient](@cite)."""
function cayley_complex_square_graphs(G,A,B,GraphType=Multigraph)
    idx_to_mat = collect(G); # TODO see if there is a better (lazy?) way to enumerate
    mat_to_idx = Dict(mat=>i for (i,mat) in pairs(idx_to_mat))

    N = length(G)
    𝒢₀□ = GraphType(N)
    𝒢₁□ = GraphType(N)
    edge₀_index = Dict{Tuple{Int,Int,Int},Int}()
    edge₁_index = Dict{Tuple{Int,Int,Int},Int}()
    count = 1
    doneset = Set{Tuple{eltype(A),eltype(B)}}()
    @showprogress for (_,g) in pairs(idx_to_mat)
        for a in A
            inva = inv(a)
            ag = a*g
            i = mat_to_idx[ag]
            for b in B
                invb = inv(b)
                (inva, invb) ∈ doneset && continue # TODO there should be a better way to avoid double counting
                push!(doneset, (a,b))
                agb = a*g*b
                @assert agb != g
                j = mat_to_idx[agb]
                add_edge!(𝒢₀□,i,j)
                edge₀_index[(minmax(i,j)...,Multigraphs.mul(𝒢₀□,i,j))] = count
                gb = g*b
                @assert ag != gb
                j = mat_to_idx[gb]
                add_edge!(𝒢₁□,i,j)
                edge₁_index[(minmax(i,j)...,Multigraphs.mul(𝒢₁□,i,j))] = count
                count += 1
            end
        end
    end
    𝒢₀□, 𝒢₁□, edge₀_index, edge₁_index
end

"""Construct the Tanner code for a given multigraph, edge numbering and local code.

The edge numbering is a map from (vertex, vertex, multiplicity) to index.
Most convenient when used with [`cayley_complex_square_graphs`](@ref).

As depicted in [dinur2022locally](@cite), [leverrier2022quantum](@cite), and [gu2022efficient](@cite)."""
function tanner_code(mgraph,edge_index,local_code)
    V = nv(mgraph)
    E = ne(mgraph, count_mul=true)
    r, Δ = size(local_code)
    code = zeros(Bool, r*V, E)
    for v in Graphs.vertices(mgraph)
        neigh = Graphs.neighbors(mgraph,v)
        col = 1
        for v2 in neigh
            multiplicity = Multigraphs.mul(mgraph,v,v2)
            for m in 1:multiplicity
                e = edge_index[(minmax(v,v2)...,m)]
                for row in 1:r
                    @assert col <= Δ
                    code[(v-1)*r+row,e] = local_code[row,col].data # TODO □.data is bad way to write this
                end
                col += 1
            end
        end
        @assert col == Δ+1
    end
    code
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

"""Check the generating set is symmetric."""
is_symmetric_gen(gens) = Set(inv.(gens)) == Set(gens)
