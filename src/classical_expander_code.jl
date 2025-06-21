"""Generate a random binary parity check matrix with specified dimensions and rank properties.
"""
function random_parity_check_matrix(rows::Int, cols::Int; min_rank::Int=1, max_attempts::Int=100)
    if rows <= 0 || cols <= 0
        error("Matrix dimensions must be positive")
    end
    if min_rank > min(rows, cols)
        error("Requested rank exceeds possible maximum")
    end
    for _ in 1:max_attempts
        H = zeros(Int, rows, cols)
        for c in 1:cols
            while true
                col = rand(0:1, rows)
                sum(col) > 0 && (H[:,c] = col; break)
            end
        end
        rk = rank(Matrix(H))
        rk >= min_rank && return H
    end
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
