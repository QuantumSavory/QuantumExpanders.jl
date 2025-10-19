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

"""Construct the Cayley complex square graphs 𝒢₀□ and 𝒢₁□ as presented in [gu2022efficient](@cite).

It is more convenient to count the edges as directional (i.e. double counting them),
as that makes it much easier to track how edge indices correspond to indices in A×B.
"""
function cayley_complex_square_graphs(G,A,B,GraphType=DiMultigraph)
    @assert is_symmetric_gen(A) "Definition 3.1: Set A must be symmetric generating set [dinur2022locally](@cite)"
    @assert is_symmetric_gen(B) "Definition 3.1: Set B must be symmetric generating set [dinur2022locally](@cite)"
    # Identity element of G is neither in A nor in B
    @assert !(one(G) in A) "Definition 3.1: Identity must not be in A [dinur2022locally](@cite)"
    @assert !(one(G) in B) "Definition 3.1: Identity must not be in B [dinur2022locally](@cite)"
    # Total No-conjugacy Condition
    @assert is_nonconjugate(G, A, B) "Definition 3.6: ∀ a ∈ A, b ∈ B, g ∈ G, g⁻¹ag ≠ b [dinur2022locally](@cite)"
    # Mappings between group element as a matrix and as an integer enumerator
    idx_to_mat = collect(G); # TODO see if there is a better (lazy?) way to enumerate
    mat_to_idx = Dict(mat=>i for (i,mat) in pairs(idx_to_mat))

    # |Q| = |G||A||B|/2 indexed by the `count` variable below.
    # |V₀| = |V₁| = |G|

    # It is convenient if the V₀ and V₁ indexing is consistent,
    # i.e. the index for (v,0)∈V₀ and for (v,1)∈V₁ should be the same.
    # The indexing function is the `mat_to_idx` map.

    # The indexing of the edges has to be consistent with
    # the indexing of Q, i.e., the indexing of |G||A||B|/2.
    # In other words, each edge should know the value of the `q_count` variable
    # for which it was generated. That is stored in the `edgeᵢ_index` maps.

    # Even more subtly, the indexing of each neighborhood of a vertex v,
    # needs to be consistent with the indexing of A×B.
    # This is why we provide two indices:
    # - an A×B index useful for ordering
    # - a larger Q index useful for assigning qubits

    N = length(G)
    𝒢₀□ = GraphType(N) # vertices V₀=G×{0}, edges Q, |A||B|-regular multigraph
    𝒢₁□ = GraphType(N) # vertices V₁=G×{1}, edges Q, |A||B|-regular multigraph
    edge₀_q_idx = Dict{Tuple{Int,Int,Int},Int}() # maps an edge (with multiplicity) to Q qubit/square index
    edge₁_q_idx = Dict{Tuple{Int,Int,Int},Int}() # maps an edge (with multiplicity) to Q qubit/square index
    edge₀_ab_idx = Dict{Tuple{Int,Int,Int},Int}() # maps an edge (with multiplicity) to AB index
    edge₁_ab_idx = Dict{Tuple{Int,Int,Int},Int}() # maps an edge (with multiplicity) to AB index
    q_count = 0
    donedict = Dict{Tuple{Int,Int,Int,Int},Int}() # used to avoid double counting
    @showprogress for (iᵍ,g) in pairs(idx_to_mat)
        iᵍ = mat_to_idx[g]
        ab_count = 0
        for (jᵃ,a) in pairs(A)
            ag = a*g
            iᵃᵍ = mat_to_idx[ag]
            for (jᵇ,b) in pairs(B)
                ab_count += 1
                agb = a*g*b
                @assert agb != g
                iᵃᵍᵇ = mat_to_idx[agb]
                gb = g*b
                @assert ag != gb
                iᵍᵇ = mat_to_idx[gb]
                # Check for double counting
                # There are squares that share one of the two diagonals, but are otherwise not the same square
                q = (minmax(iᵍ,iᵃᵍᵇ)...,minmax(iᵍᵇ,iᵃᵍ)...)
                #q = NTuple{4,Int}(sort([iᵍ,iᵃᵍᵇ,iᵍᵇ,iᵃᵍ]))
                if !haskey(donedict,q)# TODO there should be a better way to avoid double counting
                    q_count+=1
                    donedict[q] = q_count
                end
                e₀ = iᵍ,iᵃᵍᵇ # the order is important
                Multigraphs.add_edge!(𝒢₀□,e₀...)
                edge₀_q_idx[(e₀...,Multigraphs.mul(𝒢₀□,e₀...))] = donedict[q]
                edge₀_ab_idx[(e₀...,Multigraphs.mul(𝒢₀□,e₀...))] = ab_count
                e₁ = iᵍᵇ,iᵃᵍ # the order is important
                Multigraphs.add_edge!(𝒢₁□,e₁...)
                edge₁_q_idx[(e₁...,Multigraphs.mul(𝒢₁□,e₁...))] = donedict[q]
                edge₁_ab_idx[(e₁...,Multigraphs.mul(𝒢₁□,e₁...))] = ab_count
            end
        end
    end
    @assert N == length(G) "Vertex sets V₀ and V₁ each have size |G| [gu2022efficient](@cite)"
    @info "|V₀| = |V₁| = |G| = $N"
    total_A_edges = N*length(A)
    total_B_edges = N*length(B)
    @info "|E_A| = Δ|G| = $(total_A_edges), |E_B| = Δ|G| = $(total_B_edges)"
    @assert total_A_edges == N*length(A) "|E_A| = Δ|G| [gu2022efficient](@cite)"
    @assert total_B_edges == N*length(B) "|E_B| = Δ|G| [gu2022efficient](@cite)"
    @info "|Q| = Δ²|G|/2 = $(q_count)"
    @assert q_count == N*length(A)*length(B)÷2 "|Q| = Δ²|G|/2 [gu2022efficient](@cite)"
    @assert unique(values(Multigraphs.indegree(𝒢₀□))) == [length(A)*length(B)] "𝒢₀□ is Δ²-regular multigraph [gu2022efficient](@cite)"
    @assert unique(values(Multigraphs.indegree(𝒢₁□))) == [length(A)*length(B)] "𝒢₁□ is Δ²-regular multigraph [gu2022efficient](@cite)"
    @assert unique(values(Multigraphs.outdegree(𝒢₀□))) == [length(A)*length(B)] "𝒢₀□ is Δ²-regular multigraph [gu2022efficient](@cite)" 
    @assert unique(values(Multigraphs.outdegree(𝒢₁□))) == [length(A)*length(B)] "𝒢₁□ is Δ²-regular multigraph [gu2022efficient](@cite)"
    # "By TNC, each square is guaranteed to have 4 distinct vertices [gu2022efficient](@cite)"
    for g in G, a in A, b in B
        vertices = [g, a*g, g*b, a*g*b]
        @assert length(Set(vertices)) == 4 "By TNC, each square has 4 distinct vertices [gu2022efficient](@cite)"
    end
    # "There are Δ² squares incident to a given vertex [gu2022efficient](@cite)"
    for g_idx in 1:N
        g = idx_to_mat[g_idx]
        incident_squares = Set()
        for a in A, b in B
            square_vertices = (mat_to_idx[g], mat_to_idx[a*g], mat_to_idx[g*b], mat_to_idx[a*g*b])
            push!(incident_squares, square_vertices)
        end
        @assert length(incident_squares) == length(A) * length(B) "Each vertex has Δ² incident squares [gu2022efficient](@cite)"
    end
    return 𝒢₀□, 𝒢₁□, edge₀_q_idx, edge₁_q_idx, edge₀_ab_idx, edge₁_ab_idx
end

"""Construct the Cayley complex square graphs 𝒢₀□ and 𝒢₁□ using the quadripartite construction as presented in [leverrier2022quantum](@cite).

The quadripartite construction removes the TNC and symmetric generator set conditions.

It is more convenient to count the edges as directional (i.e. double counting them),
as that makes it much easier to track how edge indices correspond to indices in A×B.
"""
function cayley_complex_square_graphs_quadripartite(G,A,B,GraphType=DiMultigraph)
    # Mappings between group element as a matrix and as an integer enumerator
    idx_to_mat = collect(G); # TODO see if there is a better (lazy?) way to enumerate
    mat_to_idx = Dict(mat=>i for (i,mat) in pairs(idx_to_mat))

    # |Q| = |G||A||B| indexed by the `count` variable below.
    # |V₀₀| = |V₀₁| = |V₁₀| = |V₁₁| = |G|

    # It is convenient if the V₀₀, V₀₁, V₁₀, and V₁₁ indexing are consistent,
    # i.e. the index for (v,00)∈V₀₀, (v,01)∈V₀₁, (v,10)∈V₁₀, and (v,11)∈V₁₁ should be the same.
    # The indexing function is the `mat_to_idx` map.

    # The indexing of the edges has to be consistent with
    # the indexing of Q, i.e., the indexing of |G||A||B|.
    # In other words, each edge should know the value of the `q_count` variable
    # for which it was generated. That is stored in the `edgeᵢ_index` maps.

    # Even more subtly, the indexing of each neighborhood of a vertex v,
    # needs to be consistent with the indexing of A×B.
    # This is why we provide two indices:
    # - an A×B index useful for ordering
    # - a larger Q index useful for assigning qubits

    N = length(G)
    𝒢₀□ = GraphType(2*N) # vertices V₀₀=G×{00} ∪ V₁₁=G×{11}, edges Q, |A||B|-regular multigraph
    𝒢₁□ = GraphType(2*N) # vertices V₀₁=G×{01} ∪ V₁₀=G×{10}, edges Q, |A||B|-regular multigraph
    edge₀_q_idx = Dict{Tuple{Int,Int,Int},Int}() # maps an edge (with multiplicity) to Q qubit/square index
    edge₁_q_idx = Dict{Tuple{Int,Int,Int},Int}() # maps an edge (with multiplicity) to Q qubit/square index
    edge₀_ab_idx = Dict{Tuple{Int,Int,Int},Int}() # maps an edge (with multiplicity) to AB index
    edge₁_ab_idx = Dict{Tuple{Int,Int,Int},Int}() # maps an edge (with multiplicity) to AB index
    q_count = 0
    @showprogress for (iᵍ,g) in pairs(idx_to_mat)
        iᵍ = mat_to_idx[g]
        ab_count = 0
        for (jᵃ,a) in pairs(A)
            ag = a*g
            iᵃᵍ = mat_to_idx[ag] + N # we add N so that iᵃᵍ is shifted from V₀₁ to V₁₀ 
            for (jᵇ,b) in pairs(B)
                ab_count += 1
                agb = a*g*b
                iᵃᵍᵇ = mat_to_idx[agb] + N # we add N so that iᵃᵍᵇ is shifted from V₀₀ to V₁₁ 
                gb = g*b
                iᵍᵇ = mat_to_idx[gb]
                q = (iᵍ,iᵃᵍᵇ,iᵍᵇ,iᵃᵍ) # note each q is unique due to the quadripartite construction
                q_count+=1
                e₀ = iᵍ,iᵃᵍᵇ # the order is important
                Multigraphs.add_edge!(𝒢₀□,e₀...)
                edge₀_q_idx[(e₀...,Multigraphs.mul(𝒢₀□,e₀...))] = q_count
                edge₀_ab_idx[(e₀...,Multigraphs.mul(𝒢₀□,e₀...))] = ab_count
                e₁ = iᵍᵇ,iᵃᵍ # the order is important
                Multigraphs.add_edge!(𝒢₁□,e₁...)
                edge₁_q_idx[(e₁...,Multigraphs.mul(𝒢₁□,e₁...))] = q_count
                edge₁_ab_idx[(e₁...,Multigraphs.mul(𝒢₁□,e₁...))] = ab_count
            end
        end
    end
    @info "|Q| = |G||A||B| = $(q_count)"
    @assert q_count==N*length(A)*length(B)
    @assert sort!(unique(values(Multigraphs.indegree(𝒢₀□)))) == [0, length(A)*length(B)]
    @assert sort!(unique(values(Multigraphs.indegree(𝒢₁□)))) == [0, length(A)*length(B)]
    @assert sort!(unique(values(Multigraphs.outdegree(𝒢₀□)))) == [0, length(A)*length(B)]
    @assert sort!(unique(values(Multigraphs.outdegree(𝒢₁□)))) == [0, length(A)*length(B)]

    𝒢₀□, 𝒢₁□, edge₀_q_idx, edge₁_q_idx, edge₀_ab_idx, edge₁_ab_idx
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
Check whether the *symmetric generating* sets A and B satisfy the *Total Non-Conjugacy* (TNC) 
condition for the given group G, as defined in [dinur2022locally](@cite).

Given a finite group with symmetric generating sets A and B (i.e., ``A = A^-1`` and ``B = B^-1``),
the TNC condition requires:

```math
\\begin{aligned}
\\forall a \\in A, \\forall b \\in B, \\forall g \\in G, \\quad g^{-1} a g \\neq b
\\end{aligned}
```

### Arguments
- `group`: The finite group G
- `genA`: Symmetric generating set A for G (``A = A^-1``)
- `genB`: Symmetric generating set B for G (``B = B^-1``)
"""
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
is_symmetric_gen(gens) = Set(Nemo.inv.(gens)) == Set(gens)
