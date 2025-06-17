using Nemo
using Oscar
using LinearAlgebra
using Random
using Graphs
using Multigraphs
using ProgressMeter

"""Construct the Cayleyʳⁱᵍʰᵗ graph for a given group and set of generators."""
function cayley_right(generators::Vector{Matrix{FqFieldElem}}, group_order::Int)
    identity_mat = [one(generators[1][1,1]) zero(generators[1][1,1]);
                   zero(generators[1][1,1]) one(generators[1][1,1])]
    graph = SimpleGraph(1)
    mat_to_idx = Dict{typeof(identity_mat), Int}()
    mat_to_idx[identity_mat] = 1
    queue = [identity_mat]
    while !isempty(queue)
        current_mat = popfirst!(queue)
        current_idx = mat_to_idx[current_mat]
        for gen in generators
            new_mat = current_mat * gen
            if !haskey(mat_to_idx, new_mat)
                new_idx = length(mat_to_idx) + 1
                mat_to_idx[new_mat] = new_idx
                add_vertex!(graph)
                push!(queue, new_mat)
            end
            target_idx = mat_to_idx[new_mat]
            add_edge!(graph, current_idx, target_idx)
        end
    end
    @assert nv(graph) == group_order "Graph size doesn't match expected group order"
    return graph
end

"""Construct the Cayleyˡᵉᶠᵗ graph for a given group and set of generators."""
function cayley_left(generators::Vector{Matrix{FqFieldElem}}, group_order::Int)
    identity_mat = [one(generators[1][1,1]) zero(generators[1][1,1]);
                   zero(generators[1][1,1]) one(generators[1][1,1])]
    graph = SimpleDiGraph(1)
    mat_to_idx = Dict{typeof(identity_mat), Int}()
    mat_to_idx[identity_mat] = 1
    queue = [identity_mat]
    while !isempty(queue)
        current_mat = popfirst!(queue)
        current_idx = mat_to_idx[current_mat]
        for gen in generators
            new_mat = gen * current_mat
            if !haskey(mat_to_idx, new_mat)
                new_idx = length(mat_to_idx) + 1
                mat_to_idx[new_mat] = new_idx
                add_vertex!(graph)
                push!(queue, new_mat)
            end
            target_idx = mat_to_idx[new_mat]
            add_edge!(graph, current_idx, target_idx)
        end
    end
    @assert nv(graph) == group_order "Graph size doesn't match expected group order"
    return graph
end

"""Construct the Cayley complex square graphs 𝒢₀□ and 𝒢₁□ using group generators."""
function cayley_complex_square_graphs(
    A::Vector{Matrix{FqFieldElem}},
    B::Vector{Matrix{FqFieldElem}},
    group_order::Int,
    GraphType=DiMultigraph
)
    identity_mat = [one(A[1][1,1]) zero(A[1][1,1]);
                  zero(A[1][1,1]) one(A[1][1,1])]
    # First construct the entire group using BFS with generators A
    mat_to_idx = Dict{typeof(identity_mat), Int}()
    mat_to_idx[identity_mat] = 1
    queue = [identity_mat]
    idx_to_mat = [identity_mat]
    while !isempty(queue)
        current_mat = popfirst!(queue)
        for a in A
            new_mat = a * current_mat
            if !haskey(mat_to_idx, new_mat)
                new_idx = length(mat_to_idx) + 1
                mat_to_idx[new_mat] = new_idx
                push!(idx_to_mat, new_mat)
                push!(queue, new_mat)
            end
        end
    end
    @assert length(mat_to_idx) == group_order "Constructed group doesn't match expected order"
    # Total no-conjugacy verification
    for g in idx_to_mat, a in A, b in B
        @assert !exact_equal(a*g, g*b) "TNC violation: a*g = g*b for a=$a, b=$b, g=$g"
    end
    # Now build the square graphs
    N = group_order
    𝒢₀□ = GraphType(N)
    𝒢₁□ = GraphType(N)
    edge₀_q_idx = Dict{Tuple{Int,Int,Int},Int}()
    edge₁_q_idx = Dict{Tuple{Int,Int,Int},Int}()
    edge₀_ab_idx = Dict{Tuple{Int,Int,Int},Int}()
    edge₁_ab_idx = Dict{Tuple{Int,Int,Int},Int}()
    q_count = 0
    donedict = Dict{Tuple{Int,Int,Int,Int},Int}()
    for (iᵍ, g) in enumerate(idx_to_mat)
        ab_count = 0
        for a in A
            ag = a * g
            iᵃᵍ = mat_to_idx[ag]
            for b in B
                ab_count += 1
                agb = a * g * b
                iᵃᵍᵇ = mat_to_idx[agb]
                gb = g * b
                iᵍᵇ = mat_to_idx[gb]
                # Check for double counting
                q = (minmax(iᵍ, iᵃᵍᵇ)..., minmax(iᵍᵇ, iᵃᵍ)...)
                if !haskey(donedict, q)
                    q_count += 1
                    donedict[q] = q_count
                end
                # Add edges to 𝒢₀□
                e₀ = (iᵍ, iᵃᵍᵇ)
                add_edge!(𝒢₀□, e₀...)
                edge₀_q_idx[(e₀..., Multigraphs.mul(𝒢₀□, e₀...))] = donedict[q]
                edge₀_ab_idx[(e₀..., Multigraphs.mul(𝒢₀□, e₀...))] = ab_count
                # Add edges to 𝒢₁□
                e₁ = (iᵍᵇ, iᵃᵍ)
                add_edge!(𝒢₁□, e₁...)
                edge₁_q_idx[(e₁..., Multigraphs.mul(𝒢₁□, e₁...))] = donedict[q]
                edge₁_ab_idx[(e₁..., Multigraphs.mul(𝒢₁□, e₁...))] = ab_count
            end
        end
    end
    @info "|Q| = |G||A||B|/2 = $(q_count)"
    @assert q_count == group_order * length(A) * length(B) ÷ 2
    @assert unique(values(indegree(𝒢₀□))) == [length(A) * length(B)]
    @assert unique(values(indegree(𝒢₁□))) == [length(A) * length(B)]
    @assert unique(values(outdegree(𝒢₀□))) == [length(A) * length(B)]
    @assert unique(values(outdegree(𝒢₁□))) == [length(A) * length(B)]
    return 𝒢₀□, 𝒢₁□, edge₀_q_idx, edge₁_q_idx, edge₀_ab_idx, edge₁_ab_idx
end

"""Construct the Cayley complex square graphs 𝒢₀□ and 𝒢₁□ using the quadripartite construction as presented in [leverrier2022quantum](@cite).

The quadripartite construction removes the TNC and symmetric generator set conditions.

It is more convenient to count the edges as directional (i.e. double counting them),
as that makes it much easier to track how edge indices correspond to indices in A×B.
"""
function cayley_complex_square_graphs_quadripartite(
    A::Vector{Matrix{FqFieldElem}},
    B::Vector{Matrix{FqFieldElem}},
    group_order::Int,
    GraphType=DiMultigraph
)
    identity_mat = [one(A[1][1,1]) zero(A[1][1,1]);
                  zero(A[1][1,1]) one(A[1][1,1])]
    # Lazy group construction using BFS with generators A
    mat_to_idx = Dict{typeof(identity_mat), Int}()
    mat_to_idx[identity_mat] = 1
    queue = [identity_mat]
    idx_to_mat = [identity_mat]
    while !isempty(queue)
        current_mat = popfirst!(queue)
        for a in A
            new_mat = a * current_mat
            if !haskey(mat_to_idx, new_mat)
                new_idx = length(mat_to_idx) + 1
                mat_to_idx[new_mat] = new_idx
                push!(idx_to_mat, new_mat)
                push!(queue, new_mat)
            end
        end
    end
    @assert length(mat_to_idx) == group_order "Constructed group doesn't match expected order"
    # Build quadripartite graphs
    N = group_order
    𝒢₀□ = GraphType(2*N)  # V₀₀ ∪ V₁₁
    𝒢₁□ = GraphType(2*N)  # V₀₁ ∪ V₁₀
    edge₀_q_idx = Dict{Tuple{Int,Int,Int},Int}()
    edge₁_q_idx = Dict{Tuple{Int,Int,Int},Int}()
    edge₀_ab_idx = Dict{Tuple{Int,Int,Int},Int}()
    edge₁_ab_idx = Dict{Tuple{Int,Int,Int},Int}()
    q_count = 0
    @showprogress for (iᵍ, g) in enumerate(idx_to_mat)
        ab_count = 0
        for a in A
            ag = a * g
            iᵃᵍ = mat_to_idx[ag] + N  # Shift to V₁₀
            for b in B
                ab_count += 1
                agb = a * g * b
                iᵃᵍᵇ = mat_to_idx[agb] + N  # Shift to V₁₁
                gb = g * b
                iᵍᵇ = mat_to_idx[gb]
                # Each q is unique in quadripartite construction
                q_count += 1
                # Add edges to 𝒢₀□ (V₀₀ → V₁₁)
                e₀ = (iᵍ, iᵃᵍᵇ)
                add_edge!(𝒢₀□, e₀...)
                edge₀_q_idx[(e₀..., Multigraphs.mul(𝒢₀□, e₀...))] = q_count
                edge₀_ab_idx[(e₀..., Multigraphs.mul(𝒢₀□, e₀...))] = ab_count
                # Add edges to 𝒢₁□ (V₀₁ → V₁₀)
                e₁ = (iᵍᵇ, iᵃᵍ)
                add_edge!(𝒢₁□, e₁...)
                edge₁_q_idx[(e₁..., Multigraphs.mul(𝒢₁□, e₁...))] = q_count
                edge₁_ab_idx[(e₁..., Multigraphs.mul(𝒢₁□, e₁...))] = ab_count
            end
        end
    end
    @info "|Q| = |G||A||B| = $(q_count)"
    @assert q_count == N * length(A) * length(B)
    @assert sort!(unique(values(indegree(𝒢₀□)))) == [0, length(A)*length(B)]
    @assert sort!(unique(values(indegree(𝒢₁□)))) == [0, length(A)*length(B)]
    @assert sort!(unique(values(outdegree(𝒢₀□)))) == [0, length(A)*length(B)]
    @assert sort!(unique(values(outdegree(𝒢₁□)))) == [0, length(A)*length(B)]
    return 𝒢₀□, 𝒢₁□, edge₀_q_idx, edge₁_q_idx, edge₀_ab_idx, edge₁_ab_idx
end

"""Construct the Tanner code for a given multigraph, edge numbering and local code.

The edge numbering is a map from (vertex, vertex, multiplicity) to index.
Most convenient when used with [`cayley_complex_square_graphs`](@ref).

As depicted in [dinur2022locally](@cite), [leverrier2022quantum](@cite), and [gu2022efficient](@cite)."""
function tanner_code(mgraph,edge_q_index,edge_ab_index,local_code)
    V = nv(mgraph)
    E = ne(mgraph, count_mul=true)÷2 # edges are double counted
    r, Δ = size(local_code)
    code = zeros(Bool, r*V, E)
    for v in Graphs.vertices(mgraph)
        neigh = neighbors(mgraph,v)
        q_indices = rem.([edge_q_index[(v,v₂,m)] for v₂ in neigh for m in 1:Multigraphs.mul(mgraph,v,v₂)] .-1, E).+1
        ab_indices = rem.([edge_ab_index[(v,v₂,m)] for v₂ in neigh for m in 1:Multigraphs.mul(mgraph,v,v₂)] .-1, E).+1
        indices = q_indices[sortperm(ab_indices)] # crucial to ensure consistent local view
        @assert length(q_indices) == Δ
        @assert length(Set(ab_indices)) == Δ
        for row in 1:r
            code[(v-1)*r+row,indices] .= [e.data for e in local_code[row,:]][1,:] # TODO there must be a neater way to write this
        end
    end
    code
end

"""Construct the Tanner code for a given multigraph, edge numbering and local code.

The edge numbering is a map from (vertex, vertex, multiplicity) to index.
Most convenient when used with [`cayley_complex_square_graphs_quadripartite`](@ref).

As depicted in [dinur2022locally](@cite), [leverrier2022quantum](@cite), and [gu2022efficient](@cite)."""
function tanner_code_quadripartite(mgraph,edge_q_index,edge_ab_index,local_code)
    V = nv(mgraph)
    E = ne(mgraph, count_mul=true) # edges are not double counted here
    r, Δ = size(local_code)
    code = zeros(Bool, r*V÷2, E)
    for v in sort!(Graphs.vertices(mgraph))[1:V÷2] # only first half of vertices have outgoing edges
        neigh = neighbors(mgraph,v)
        q_indices = rem.([edge_q_index[(v,v₂,m)] for v₂ in neigh for m in 1:Multigraphs.mul(mgraph,v,v₂)] .-1, E).+1
        ab_indices = rem.([edge_ab_index[(v,v₂,m)] for v₂ in neigh for m in 1:Multigraphs.mul(mgraph,v,v₂)] .-1, E).+1
        indices = q_indices[sortperm(ab_indices)] # crucial to ensure consistent local view
        @assert length(q_indices) == Δ
        @assert length(Set(ab_indices)) == Δ
        for row in 1:r 
            code[(v-1)*r+row,indices] .= [e.data for e in local_code[row,:]][1,:] # TODO there must be a neater way to write this
        end
    end
    code
end

"""
    is_nonconjugate(group_generators, group_order, genA, genB)

Check the Totally Non-Conjugate (TNC) condition from [`dinur2022locally`](@cite) using lazy group enumeration.

The TNC condition requires that no element of `genA` is conjugate to any element of `genB` under the group action.

# Example

```jldoctest
julia> using QuantumExpanders; using Oscar

julia> using QuantumExpanders: FirstOnly, AllPairs

julia> generators, order = morgenstern_generators(1, 4)
(Matrix{FqFieldElem}[[o^3 + o^2 + 1 o^3 + o^2 + 1; o^3 + 1 o^3 + o^2 + 1], [o^3 + o^2 + 1 o^3; o^3 + o o^3 + o^2 + 1], [o^3 + o^2 + 1 o^2 + 1; o + 1 o^3 + o^2 + 1]], 4080)

julia> B = generators[1:end-1]
2-element Vector{Matrix{FqFieldElem}}:
 [o^3 + o^2 + 1 o^3 + o^2 + 1; o^3 + 1 o^3 + o^2 + 1]
 [o^3 + o^2 + 1 o^3; o^3 + o o^3 + o^2 + 1]

julia> A = alternative_morgenstern_generators(B, FirstOnly())
2-element Vector{Matrix{FqFieldElem}}:
 [o^2 + 1 o^3 + o^2; o^2 o^3 + o]
 [o^3 + o o^3 + o^2; o^2 o^2 + 1]

julia> Al = alternative_morgenstern_generators(B, AllPairs())
2-element Vector{Matrix{FqFieldElem}}:
 [o^2 + 1 o^3 + o^2; o^2 o^3 + o]
 [o^3 + o o^3 + o^2; o^2 o^2 + 1]

julia> is_nonconjugate(generators, order, A, B)
true

julia> is_nonconjugate(generators, order, Al, B)
true

julia> !is_nonconjugate(generators, order, A, A)
true
```
function is_nonconjugate(
    group_generators::Vector{Matrix{FqFieldElem}},
    group_order::Int,
    genA::Vector{Matrix{FqFieldElem}},
    genB::Vector{Matrix{FqFieldElem}}
)
    F = parent(group_generators[1][1,1])
    MS = matrix_space(F, 2, 2)
    identity_mat = MS(1)
    group_gens_nemo = [MS(g) for g in group_generators]
    genA_nemo = [MS(a) for a in genA]
    genB_nemo = [MS(b) for b in genB]
    genset = Set(genB_nemo)
    discovered = Set([identity_mat])
    queue = [identity_mat]
    while !isempty(queue)
        g = popfirst!(queue)
        for a in genA_nemo
            inv_g = inv(g)
            conjugated = inv_g * a * g
            if conjugated ∈ genset
                return false
            end
        end
        for gen in group_gens_nemo
            new_elem = g * gen
            if new_elem ∉ discovered
                push!(discovered, new_elem)
                push!(queue, new_elem)
                
                if length(discovered) > group_order
                    error("Discovered more elements than group order")
                end
            end
        end
    end
    @assert length(discovered) == group_order "Group enumeration incomplete"
    return true
end

"""Check the generating set is symmetric."""
is_symmetric_gen(gens) = Set(inv.(gens)) == Set(gens)

"""Exact equality check for matrices over finite fields"""
function exact_equal(a::Nemo.FqMatrix, b::Nemo.FqMatrix)
    n = size(a, 1)
    return all(a[i,j] == b[i,j] for i in 1:n, j in 1:n)
end
function exact_equal(a::Matrix{Nemo.FqFieldElem}, b::Matrix{Nemo.FqFieldElem})
    n = size(a, 1)
    return all(a[i,j] == b[i,j] for i in 1:n, j in 1:n)
end

"""Check if a generating set is symmetric (closed under inversion)"""
function is_symmetric_gen(gens::Vector{Matrix{Nemo.FqFieldElem}})
    # Convert all matrices to FqMatrix type first
    R = parent(gens[1][1,1])
    n = size(gens[1], 1)
    M = matrix_space(R, n, n)
    gens_fq = [M(g) for g in gens]  # Convert to FqMatrix
    
    inv_gens = [Nemo.inv(g) for g in gens_fq]
    
    for g_inv in inv_gens
        if !any(exact_equal(g_inv, g) for g in gens_fq)
            return false
        end
    end
    return true
end

"""Create symmetric generating sets of equal size"""
function equalize_generators(A, B)
    Δ_A = length(A)
    Δ_B = length(B)
    if Δ_A == Δ_B
        return A, B
    elseif Δ_A < Δ_B
         return [A; A[1:Δ_B-Δ_A]], B
    else
         return A, [B; B[1:Δ_A-Δ_B]]
    end
end

"""Check if a graph is Ramanujan by verifying the eigenvalue condition."""
function is_ramanujan(g, q)
    A = adjacency_matrix(g)
    λ = eigvals(Matrix(A))
    λ_sorted = sort(λ, rev=true)
    r = q + 1
    λ_second = λ_sorted[2]
    return λ_second <= 2 * sqrt(q)
end
