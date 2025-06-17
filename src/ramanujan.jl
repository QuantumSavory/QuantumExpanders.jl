using Oscar
using LinearAlgebra
using Graphs
using Graphs: nv, add_edge!, adjacency_matrix, degree, vertices, edges, Edge, nv, ne, src, dst

"""Compute the Legendre symbol (a/p) for an odd prime p."""
function legendre_symbol(a::Int, p::Int)
    @assert is_prime(p) "p must be prime"
    ls = powermod(a, (p - 1) ÷ 2, p)
    return ls == p - 1 ? -1 : ls
end

"""Scalar matrices for GL(2,F): every nonzero scalar gives a valid scalar matrix."""
function scalar_matrices_GL(GL2)
    F = base_ring(GL2)
    return [GL2([x 0; 0 x]) for x in F if x != 0]
end

"""Scalar matrices for SL(2,F): only those with x^2 == 1 have determinant 1."""
function scalar_matrices_SL(SL2)
    F = base_ring(SL2)
    return [SL2([x 0; 0 x]) for x in F if x^2 == one(F)]
end

"""Solve p = a² + b² + c² + d² (over the integers)"""
function solve_four_squares(p::Int)
    solutions = Tuple{Int,Int,Int,Int}[]
    max_val = isqrt(p)
    for w in 0:max_val, x in 0:max_val, y in 0:max_val
        rem = p - (w^2 + x^2 + y^2)
        rem ≥ 0 || continue
        z = isqrt(rem)
        z^2 == rem || continue
        # Generate sign permutations for nonzero components.
        signs_w = w == 0 ? [0] : [1, -1]
        signs_x = x == 0 ? [0] : [1, -1]
        signs_y = y == 0 ? [0] : [1, -1]
        signs_z = z == 0 ? [0] : [1, -1]
        for sw in signs_w, sx in signs_x, sy in signs_y, sz in signs_z
            push!(solutions, (sw*w, sx*x, sy*y, sz*z))
        end
    end
    return unique(solutions)
end

"""Filter solutions: select those with a > 0 and b, c, d even."""
function process_solutions(solutions, p)
    filtered = filter(sol -> sol[1] > 0 && all(iseven, sol[2:4]), solutions)
    @assert length(filtered) == p + 1 "Incorrect number of solutions"
    return filtered
end

"""Create generator matrices over F from four–square solutions.
Each matrix has determinant equal to F(p)."""
function create_generators(solutions, F, p)
    u = sqrt(F(-1))  # Find u such that u² = -1 in F.
    generators = MatrixElem{typeof(F(0))}[]
    for (a, b, c, d) in solutions
        mat = matrix(F, [a + u*b   c + u*d;
                         -c + u*d   a - u*b])
        @assert det(mat) == F(p) "Generator matrix must have determinant p"
        push!(generators, mat)
    end
    return generators
end

""" Construct the Cayley graph from a given group quotient PG,
a list of generators (as group elements), and the canonical morphism."""
function construct_cayley_graph(PG, generators, morphism)
    element_dict = Dict{eltype(PG), Int}()
    elements = [one(PG)]
    element_dict[elements[1]] = 1
    queue = copy(elements)
    
    while !isempty(queue)
        current = popfirst!(queue)
        for gen in generators
            # Ensure gen is interpreted as a group element.
            gen_elem = parent(gen)(gen)
            new_elem = morphism(gen_elem) * current
            if !haskey(element_dict, new_elem)
                push!(elements, new_elem)
                element_dict[new_elem] = length(elements)
                push!(queue, new_elem)
            end
        end
    end

    g = SimpleGraph(length(elements))
    for (idx, elem) in enumerate(elements)
        for gen in generators
            gen_elem = parent(gen)(gen)
            neighbor = morphism(gen_elem) * elem
            neighbor_idx = element_dict[neighbor]
            add_edge!(g, idx, neighbor_idx)
        end
    end
    return g
end

"""
Check the Ramanujan property:
- For a (p+1)-regular graph, the trivial eigenvalue is p+1.
- All other eigenvalues should have absolute value ≤ 2√p.
"""
function is_ramanujan(g::SimpleGraph, p::Int)
    A = adjacency_matrix(g)
    λ = sort(eigvals(Matrix(A)), rev=true)
    bound = 2 * sqrt(p)
    non_trivial = filter(x -> abs(x - (p+1)) > 1e-6, λ)
    return all(v -> abs(v) ≤ bound + 1e-6, non_trivial)
end


"""
Construct the Ramanujan graph X^(p,q).
Chooses the group quotient based on the Legendre symbol (p/q):
- If (p/q) = -1, use PGL₂(F).
- If (p/q) = 1, use PSL₂(F) (scaling generators so their determinant becomes 1).
"""
function ramanujan_graph(p::Int, q::Int)
    @assert is_prime(p) && is_prime(q) "p and q must be primes."
    @assert p % 4 == 1 && q % 4 == 1 "p and q must be ≡ 1 mod 4."
    
    F = GF(q)
    symbol = legendre_symbol(p, q)  # Compute (p/q)
    
    if symbol == -1
        # Use GL(2,F) → PGL₂(F)
        GL2 = GL(2, F)
        center = scalar_matrices_GL(GL2)
        PG, morphism = quo(GL2, center)
        solutions = solve_four_squares(p)
        solutions_processed = process_solutions(solutions, p)
        generators = create_generators(solutions_processed, F, p)
        gl_gens = [GL2(mat) for mat in generators]
        return construct_cayley_graph(PG, gl_gens, morphism)
    elseif symbol == 1
        # Use SL(2,F) → PSL₂(F)
        SL2 = SL(2, F)
        center = scalar_matrices_SL(SL2)
        PG, morphism = quo(SL2, center)
        solutions = solve_four_squares(p)
        solutions_processed = process_solutions(solutions, p)
        generators = create_generators(solutions_processed, F, p)
        # In these generators, det(mat) == F(p). Since (p/q)=1, p is a square in F.
        # Let s be a square root of F(p) (note F(p) means the image of p in GF(q)).
        s = sqrt(F(p))
        s_inv = inv(s)
        # Scale each generator so that its determinant becomes 1:
        # (s_inv)^2 * p = 1.
        generators_scaled = [s_inv * mat for mat in generators]
        sl_gens = [SL2(mat) for mat in generators_scaled]
        return construct_cayley_graph(PG, sl_gens, morphism)
    else
        error("Unexpected Legendre symbol value.")
    end
end


"""
Constructs the edge-vertex incidence graph of a given graph G.
Let G = (V, E) be a graph with vertex set V and edge set E. The
edge-vertex incidence graph of G is a bipartite graph with vertex
set E ∪ V and edge set {(e, v) ∈ E × V : v is an endpoint of e}.

# Example
```jldoctest
julia> G = ramanujan_graph(5, 13);

julia> B = edge_vertex_incidence_graph(G);

julia> is_unbalanced_bipartite(B)
true
```
"""
function edge_vertex_incidence_graph(G)
    n = nv(G)
    m = ne(G)
    B = SimpleGraph(n + m)
    edge_to_index = Dict{Graphs.Edge, Int}()
    for (i, e) in enumerate(edges(G))
        edge_to_index[e] = n + i
    end
    for (e, idx) in edge_to_index
        add_edge!(B, src(e), idx)
        add_edge!(B, dst(e), idx)
    end
    return B
end

"""Verify that the edge-vertex incidence graph B is an unbalanced bipartite graph."""
function is_unbalanced_bipartite(B)
    # The first n vertices in B correspond to the vertices of G
    n = 0
    while n + 1 <= nv(B) && degree(B, n + 1) > 2
        n += 1
    end
    # If no vertices satisfy degree > 2, return false
    if n == 0
        return false
    end
    # The degree of any vertex in the first n vertices of B is d
    d = degree(B, 1)
    # The remaining vertices in B correspond to the edges of G
    m = nv(B) - n
    for v in 1:n
        if degree(B, v) != d
            return false
        end
    end
    for v in (n + 1):(n + m)
        if degree(B, v) != 2
            return false
        end
    end
    if n * d != m * 2
        return false
    end
    return true
end

"""
Test Alon–Chung Lemma, which states that for a `d`-regular graph `G` with `n`
vertices and second-largest eigenvalue `λ`, the number of edges in the subgraph
induced by any subset `X` of size `γ*n` is at most:
- (d * n / 2) * (γ² + (λ / d) * γ * (1 - γ))

# Example
```jldoctest
julia> G = ramanujan_graph(5, 13);

julia> B = edge_vertex_incidence_graph(G);

julia> alon_chung_lemma(B, 0.1)
(true, 143, 1520.6609615130378)
```
"""
function alon_chung_lemma(g, γ)
    n = nv(g)
    d = degree(g, 1)
    A = adjacency_matrix(g)
    Λ = eigvals(Matrix(A))
    λ = sort(Λ, rev=true)[2]
    # Select a random subset X of size γ * n
    subset_size = Int(round(γ * n))
    X = randperm(n)[1:subset_size]
    # Count the number of edges in the subgraph induced by X
    induced_subgraph_result = induced_subgraph(g, X)
    induced_subgraph_g = induced_subgraph_result[1]
    num_edges_in_subgraph = ne(induced_subgraph_g)
    bound = (d * n / 2) * (γ^2 + (λ / d) * γ * (1 - γ))
    satisfies_bound = num_edges_in_subgraph ≤ bound
    return satisfies_bound, num_edges_in_subgraph, bound
end

"""
    expander_code_parity_matrix

Construct the parity-check matrix for a classical expander code (a type of Tanner code).

An expander code is defined by:
- A **regular base graph** `g` (typically a good expander graph like a Ramanujan graph)
- A **short inner code** specified by its parity-check matrix `H`

The resulting code has:
- Variables corresponding to edges of `g`
- Constraints enforcing the inner code at each vertex of `g`

The construction follows Sipser-Spielman expander codes:
- Each edge in `g` becomes a variable in the code
- Each vertex enforces its incident edges to satisfy the inner code
- The code has rate ≥ 1 - (1 - r) * (n_edges/n_vertices) where r is the rate of the inner code.
"""
function expander_code_parity_matrix(g::AbstractGraph, H_inner::AbstractMatrix{<:Integer})
    H_inner = sparse(H_inner)
    n_vertices = nv(g)
    n_edges = ne(g)
    n_vertices == 0 && error("Graph must have vertices")
    degrees = degree(g)
    graph_degree = first(degrees)
    all(==(graph_degree), degrees) || error("Graph must be regular")
    _, inner_code_length = size(H_inner)
    inner_code_length == graph_degree || error("Inner code length must match graph degree")
    edge_list = collect(edges(g))
    edge_indices = Dict{Edge,Int}()
    for (idx, e) in enumerate(edge_list)
        edge_indices[e] = idx
    end
    nnz_est = n_vertices * nnz(H_inner)
    I = sizehint!(Int[], nnz_est)
    J = sizehint!(Int[], nnz_est)
    V = sizehint!(Int[], nnz_est)
    for v in 1:n_vertices
        incident_edges = [e for e in edge_list if src(e) == v || dst(e) == v]
        sort!(incident_edges, by=e->(min(src(e),dst(e)), max(src(e),dst(e))))
        for (check_row, row_vals) in enumerate(eachrow(H_inner))
            for (edge_pos, val) in enumerate(row_vals)
                val == 0 && continue
                edge_idx = edge_indices[incident_edges[edge_pos]]
                push!(I, (v-1)*size(H_inner,1) + check_row)
                push!(J, edge_idx)
                push!(V, val)
            end
        end
    end
    sparse(I, J, V, n_vertices*size(H_inner,1), n_edges)
end
