@testitem "Test classical  expander graph properties" begin
    using Oscar
    using LinearAlgebra
    using Graphs
    using Graphs: nv, neighbors, AbstractGraph, degree, ne
    using QuantumExpanders
    using LogExpFunctions

    function generate_parity_matrix(rows::Int, cols::Int; min_rank::Int=1, max_attempts::Int=100)
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

    @testset "Classical Expander code using Ramanujan Graphs" begin
        test_pairs = [(5, 29),
                      (13, 17),
                      (29, 13)]
        for (p, q) in test_pairs
            g = ramanujan_graph(p, q)
            H_inner = generate_parity_matrix(2, degree(g, 1))
            @test size(H_inner, 2) == degree(g, 1)
            H_expander = QuantumExpanders.expander_code_parity_matrix(g, H_inner)
            m, n = size(H_expander)
            r = rank(H_expander)
            k = n - r
            @test n == ne(g)
            @test m == nv(g) * size(H_inner, 1)
            @test k == n - r
            d = degree(g,1)
            inner_rate = (d - rank((H_inner))) / d
            theoretical_min_rate = 2*inner_rate - 1
            actual_rate = k/n
            @test actual_rate ≥ theoretical_min_rate
        end
    end
end
