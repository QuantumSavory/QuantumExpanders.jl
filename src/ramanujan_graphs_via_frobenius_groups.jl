"""A Cayley subset S is called normal if it is a union of conjugacy classes of G."""
function normal_cayley_subset(G)
    classes = conjugacy_classes(G)
    identity_class = nothing
    for cls in classes
        if length(cls) == 1 && first(cls) == one(G)
            identity_class = cls
            break
        end
    end
    non_trivial_classes = filter(cls -> cls != identity_class, classes)
    S = union(non_trivial_classes...)
    return collect(S)
end

function adjacency_matrix(g::SimpleGraph)
    n = nv(g)
    A = zeros(Int, n, n)
    for e in edges(g)
        i, j = src(e), dst(e)
        A[i, j] = 1
        A[j, i] = 1
    end
    return A
end

function is_ramanujan(g, k)
    A = adjacency_matrix(g)
    evals = real.(eigvals(A))
    sorted_evals = sort(evals, rev=true)
    λ1 = sorted_evals[1]
    non_trivial_evals = sorted_evals[2:end]
    μ = maximum(abs.(non_trivial_evals))
    ramanujan_bound = 2 * sqrt(k - 1)
    return μ <= ramanujan_bound + 1e-10
end
