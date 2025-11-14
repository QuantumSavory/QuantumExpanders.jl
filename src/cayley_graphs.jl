"""Construct the Cayleyʳⁱᵍʰᵗ graph for a given group and set of generators."""
function cayley_right(group::Group, generators::Vector{<:GroupElem})
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
function cayley_left(group::Group, generators::Vector{<:GroupElem})
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

Returns `(𝒢₀□, 𝒢₁□, edge₀_q_idx, edge₁_q_idx, edge₀_ab_idx, edge₁_ab_idx)` where:
- `𝒢₀□`: Square graph on ``V_0 = G × {0}`` with edges between vertices connected via squares
- `𝒢₁□`: Square graph on ``V_1 = G × {1}`` with edges between vertices connected via squares  
- ``edge₀_q_idx`: Dictionary mapping edges `(src,dst,multiplicity)` in 𝒢₀□ to their corresponding square index in Q
- `edge₁_q_idx`: Dictionary mapping edges `(src,dst,multiplicity)` in 𝒢₁□  to their corresponding square index in Q
- `edge₀_ab_idx`: Dictionary mapping edges `(src,dst,multiplicity)` in 𝒢₀□ to their corresponding position in A×B grid for 𝒢₀□
- `edge₁_ab_idx`: Dictionary mapping edges `(src,dst,multiplicity)` in 𝒢₁□ to their corresponding position in A×B grid for 𝒢₁□

# Bipartite Left-Right Cayley Complex

The bipartite left-right Cayley complex a 2-dimensional complex from a finite group G and symmetric generating
sets ``A = A^-1``, ``B = B^-1``.

The vertex set V is bipartite and partitioned as V = V₀ ∪ V₁, with ``V_0 = G × {0}`` and ``V_1 = G × {1}`` representing
two copies of the group G. The edge sets consist of A-edges ``E_A`` and B-edges ``E_B``, where ``E_A`` contains pairs
``\\{(g,0), (ag,1)\\}`` for all g ∈ G and a ∈ A, and ``E_B`` contains pairs ``\\{(g,0), (gb,1)\\}`` for all g ∈ G
and b ∈ B. The graph ``G_A = (V, E_A)`` is the double cover of the left Cayley graph Cay(G, A), while ``G_B = (V, E_B)``
is the double cover of the right Cayley graph Cay(G, B) [leverrier2022quantum](@cite).

The set Q of squares is defined as the collection of 4-subsets of vertices of the form

```math
\\begin{aligned}
{(g,0), (ag,1), (gb,1), (agb,0)}
\\end{aligned}
```

for all g ∈ G, a ∈ A, and b ∈ B. Each square contains two vertices from V₀ and two from V₁, forming the
two-dimensional cells of the complex.

The Total No-Conjugacy (TNC) condition ag ≠ gb for all a ∈ A, b ∈ B, g ∈ G ensures that every square consists
of four distinct vertices and that the local view Q(v) of squares incident to any vertex v naturally identifies
with the product set A × B [leverrier2022quantum](@cite).

By restricting to vertices in V₀, the set of squares Q defines a graph 𝒢₀□ = (V₀, Q) where edges connect
pairs (g,0) and (agb,0) that appear as opposite corners of squares. Similarly, restricting to V₁ defines
the graph 𝒢₁□ = (V₁, Q) where edges connect pairs (ag,1) and (gb,1). Both 𝒢₀□ and 𝒢₁□ are Δ²-regular multigraphs
on |G| vertices, with the total number of squares given by |Q| = Δ²|G|/2.
"""
function cayley_complex_square_graphs(G::Group, A::Vector{<:GroupElem}, B::Vector{<:GroupElem}, GraphType=DiMultigraph)
    @assert is_symmetric_gen(A) "Definition 3.1: Set A must be symmetric generating set [dinur2022locally](@cite)"
    @assert is_symmetric_gen(B) "Definition 3.1: Set B must be symmetric generating set [dinur2022locally](@cite)"
    # Identity element of G is neither in A nor in B
    @assert !(one(G) in A) "Definition 3.1: Identity must not be in A [dinur2022locally](@cite)"
    @assert !(one(G) in B) "Definition 3.1: Identity must not be in B [dinur2022locally](@cite)"
    # Total No-conjugacy Condition
    @assert is_nonconjugate(G, A, B) "Definition 3.6: ∀ a ∈ A, b ∈ B, g ∈ G, g⁻¹ag ≠ b [dinur2022locally](@cite)"
    # "By TNC, each square is guaranteed to have 4 distinct vertices [gu2022efficient](@cite)."
    for g in G, a in A, b in B
        vertices = [g, a*g, g*b, a*g*b]
        @assert length(Set(vertices)) == 4 "By TNC, each square has 4 distinct vertices [gu2022efficient](@cite)"
    end
    # Mappings between group element as a matrix and as an integer enumerator
    idx_to_mat = collect(G); # TODO see if there is a better (lazy?) way to enumerate
    mat_to_idx = Dict(mat=>i for (i,mat) in pairs(idx_to_mat))
    # "There are Δ² squares incident to a given vertex, and the set of faces incident
    # to a given vertex can be naturally identified with the set A × B [gu2022efficient](@cite)."
    N = length(G)
    for g_idx in 1:N
        g = idx_to_mat[g_idx]
        incident_squares = Set()
        for a in A, b in B
            square_vertices = (mat_to_idx[g], mat_to_idx[a*g], mat_to_idx[g*b], mat_to_idx[a*g*b])
            push!(incident_squares, square_vertices)
        end
        @assert length(incident_squares) == length(A)*length(B) "Each vertex has Δ² incident squares [gu2022efficient](@cite)"
    end
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
    for g in G, a in A, b in B
        # "q ∈ Q is present as an edge (v,v') in 𝒢̂□₀ if and only if v and v' appear as opposite V₀-corners of the square q" [gu2022efficient](@cite)
        v₀₁ = mat_to_idx[g]
        v₀₂ = mat_to_idx[a*g*b]
        @assert has_edge(𝒢₀□, v₀₁, v₀₂) || has_edge(𝒢₀□, v₀₂, v₀₁) "Edge in 𝒢̂□₀ connects opposite V₀-corners of square [gu2022efficient](@cite)"
        # "Each face q ∈ Q can be identified with its diagonal connecting its corners in V₁" [gu2022efficient](@cite)
        v₁₁ = mat_to_idx[a*g]
        v₁₂ = mat_to_idx[g*b]
        @assert has_edge(𝒢₁□, v₁₁, v₁₂) || has_edge(𝒢₁□, v₁₂, v₁₁) "Edge in 𝒢̂□₁ connects opposite V₁-corners of square [gu2022efficient](@cite)"
    end
    return 𝒢₀□, 𝒢₁□, edge₀_q_idx, edge₁_q_idx, edge₀_ab_idx, edge₁_ab_idx
end

"""Construct the Cayley complex square graphs 𝒢₀□ and 𝒢₁□ using the quadripartite construction as presented in [leverrier2022quantum](@cite).

Returns `(𝒢₀□, 𝒢₁□, edge₀_q_idx, edge₁_q_idx, edge₀_ab_idx, edge₁_ab_idx)` where:
- `𝒢₀□`: Square graph on `V₀₀ ∪ V₁₁ with edges
- `𝒢₁□`: Square graph on `V₀₁ ∪ V₁₀ with edges
- `edge₀_q_idx`: Dictionary mapping `(src,dst,multiplicity)` to their corresponding square index in 𝒢₀□
- `edge₁_q_idx`: Dictionary mapping `(src,dst,multiplicity)` to their corresponding square index in 𝒢₁□
- `edge₀_ab_idx`: Dictionary mapping `(src,dst,multiplicity)` to their corresponding position in A×B grid for 𝒢₀□
- `edge₁_ab_idx`: Dictionary mapping `(src,dst,multiplicity)` to their corresponding position in A×B grid for 𝒢₁□

!!! note
    The quadripartite construction eliminates the need for the Total No-Conjugacy and symmetric generating set
    conditions required in the bipartite version, while maintaining the essential properties needed for the
    quantum code construction [leverrier2022quantum](@cite).

It is more convenient to count the edges as directional (i.e. double counting them),
as that makes it much easier to track how edge indices correspond to indices in A×B.

# Quadripartite Left-Right Cayley Complex

The quadripartite left-right Cayley complex is built from a finite group G with two generating
sets A and B. The vertex set V is partitioned into four disjoint copies of the group:

```math
\\begin{aligned}
V = V_{00} \\cup V_{01} \\cup V_{10} \\cup V_{11}
\\end{aligned}
```

where ``V_{ij} = G \\times \\{i,j\\} for i,j \\in \\{0,1\\}``

The squares Q of the complex are defined as the set of 4-tuples

```math
\\begin{aligned}
{(g,00), (ag,01), (gb,10), (agb,11) : g \\in G, a \\in A, b \\in B\\}
\\end{aligned}
```

with total cardinality ``|Q| = |G||A||B|``.

Two square graphs are derived from this complex structure. The graph

```math
\\begin{aligned}
\\mathcal{G}_0^\\square = (V_{00} \\cup V_{11}, Q)
\\end{aligned}
```

connects vertices (g,00) to (agb,11) for each square, while

```math
\\begin{aligned}
\\mathcal{G}_1^\\square = (V_{01} \\cup V_{10}, Q)
\\end{aligned}
```

connects vertices (gb,10) to (ag,01). Both graphs are (``|A| \\times |B|``)-regular directed multigraphs,
with each vertex having exactly ``\\Delta_A \\Delta_B`` incident edges, where ``\\Delta_A = |A|`` and ``\\Delta_B = |B|``.
The local view Q(v) at any vertex v identifies with ``A \\times B``, where the squares incident to v are
in bijection with pairs ``(a,b) \\in A \\times B`` [leverrier2022quantum](@cite).
"""
function cayley_complex_square_graphs_quadripartite(G::Group, A::Vector{<:GroupElem}, B::Vector{<:GroupElem}, GraphType=DiMultigraph)
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

function _to_bool_matrix(matrix)
    rows, cols = size(matrix)
    bool_matrix = zeros(Bool, rows, cols)
    for i in 1:rows, j in 1:cols
        bool_matrix[i, j] = Bool(Int(lift(matrix[i, j])))
    end
    return bool_matrix
end

"""Construct the Tanner code for a given multigraph, edge numbering and local code.

The edge numbering is a map from (vertex, vertex, multiplicity) to index.
Most convenient when used with [`cayley_complex_square_graphs`](@ref).

Returns a binary matrix of size `(r*V) × E` representing the parity check matrix of the
Tanner code, where `r` is the number of checks in the local code, `V` is the number of vertices, and `E` is the number of edges.

# Tanner Code

Tanner code is a classical code where *bits* are placed on *edges* of a graph and
*constraints* are imposed by local codes at each vertex.

The Tanner code construction is defined as:

```math
\\begin{aligned}
T(\\mathcal{G}, C_0) = \\{ x \\in \\mathbb{F}_2^E : \\operatorname{res}_{E(v)}(x) \\in C_0 \\forall v \\in V(\\mathcal{G})) \\}
\\end{aligned}
```

the set of vectors ``\\mathbb{F_2}^E`` such that the restriction of `x` to edges incident to each vertex
`v` belongs to the local code C₀. The qubits are placed on squares of the left-right Cayley complex, with
Z-stabilizers defined by T(𝒢₀□, C_A ⊗ C_B) and X-stabilizers defined by T(𝒢₁□, C_A⊥ ⊗ C_B⊥). The commuting
condition `H_X H_Zᵀ = 0` essential for CSS codes follows naturally from the incidence structure of the complex.

As depicted in [dinur2022locally](@cite), [leverrier2022quantum](@cite), and [gu2022efficient](@cite).

### Arguments
- `mgraph`: A multigraph representing the underlying graph structure. For quantum Tanner codes, this is `𝒢₀□` or `𝒢₁□` from the left-right Cayley complex.
- `edge_q_index`: A dictionary mapping `(vertex, vertex, multiplicity)` tuples to qubit indices. This identifies which physical qubit (placed on squares/faces) corresponds to each edge in the multigraph.
- `edge_ab_index`: A dictionary mapping `(vertex, vertex, multiplicity)` tuples to local coordinate indices. This provides the identification of each edge with an element of `A×B` in the local view.
- `local_code`: A binary matrix representing the parity check matrix of the local code. For quantum Tanner codes, this is C₀ = C_A ⊗ C_B for Z-stabilizers or C₁ = C_A⊥ ⊗ C_B⊥ for X-stabilizers.
"""
function tanner_code(mgraph::Multigraphs.DiMultigraph{Int64}, edge_q_index::Dict{Tuple{Int64, Int64, Int64}, Int64}, edge_ab_index::Dict{Tuple{Int64, Int64, Int64}, Int64}, local_code)
    V = nv(mgraph)
    E = ne(mgraph, count_mul=true)÷2 # edges are double counted
    r, Δ = size(local_code)
    code = zeros(Bool, r*V, E)
    local_code = eltype(local_code) <: Bool ? local_code : _to_bool_matrix(local_code)
    for v in Graphs.vertices(mgraph)
        neigh = neighbors(mgraph,v)
        q_indices = rem.([edge_q_index[(v,v₂,m)] for v₂ in neigh for m in 1:Multigraphs.mul(mgraph,v,v₂)] .-1, E).+1
        ab_indices = rem.([edge_ab_index[(v,v₂,m)] for v₂ in neigh for m in 1:Multigraphs.mul(mgraph,v,v₂)] .-1, E).+1
        indices = q_indices[sortperm(ab_indices)] # crucial to ensure consistent local view
        @assert length(q_indices) == Δ
        @assert length(Set(ab_indices)) == Δ
        for row in 1:r
            code[(v-1)*r+row,indices] .= local_code[row,:]
        end
    end
    code
end

"""Construct the Tanner code for a given multigraph, edge numbering and local code.

The edge numbering is a map from (vertex, vertex, multiplicity) to index.
Most convenient when used with [`cayley_complex_square_graphs_quadripartite`](@ref).

As depicted in [dinur2022locally](@cite), [leverrier2022quantum](@cite), and [gu2022efficient](@cite)."""
function tanner_code_quadripartite(mgraph::Multigraphs.DiMultigraph{Int64}, edge_q_index::Dict{Tuple{Int64, Int64, Int64}, Int64}, edge_ab_index::Dict{Tuple{Int64, Int64, Int64}, Int64}, local_code)
    V = nv(mgraph)
    E = ne(mgraph, count_mul=true) # edges are not double counted here
    r, Δ = size(local_code)
    code = zeros(Bool, r*V÷2, E)
    local_code = eltype(local_code) <: Bool ? local_code : _to_bool_matrix(local_code)
    for v in sort!(Graphs.vertices(mgraph))[1:V÷2] # only first half of vertices have outgoing edges
        neigh = neighbors(mgraph,v)
        q_indices = rem.([edge_q_index[(v,v₂,m)] for v₂ in neigh for m in 1:Multigraphs.mul(mgraph,v,v₂)] .-1, E).+1
        ab_indices = rem.([edge_ab_index[(v,v₂,m)] for v₂ in neigh for m in 1:Multigraphs.mul(mgraph,v,v₂)] .-1, E).+1
        indices = q_indices[sortperm(ab_indices)] # crucial to ensure consistent local view
        @assert length(q_indices) == Δ
        @assert length(Set(ab_indices)) == Δ
        for row in 1:r 
            code[(v-1)*r+row,indices] .= local_code[row,:]
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
function is_nonconjugate(group::Group, genA::Vector{<:GroupElem}, genB::Vector{<:GroupElem})
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
is_symmetric_gen(gens::Vector{<:GroupElem}) = Set(Nemo.inv.(gens)) == Set(gens)
