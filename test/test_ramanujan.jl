@testitem "Test Ramanujan Graph properties" begin
    using Oscar
    using LinearAlgebra
    using Graphs
    using Graphs: nv, neighbors, AbstractGraph, degree, ne
    using QuantumExpanders
    using LogExpFunctions

    function binary_entropy(α)
        if α == 0.0 || α == 1.0
            return 0.0
        else
            return -α * log(α) - (1 - α) * log(1 - α)
        end
    end

    function expansion_bound(α, c, d)
        term1 = (c / d) * (1 - (1 - α)^d)
        term2 = (2 * c * α * binary_entropy(α)) / log(2)
        return term1 - term2
    end

    function verify_expansion_property(B, α)
        @assert is_unbalanced_bipartite(B)
        n = div(nv(B), 2)
        m = nv(B) - n
        d = maximum(degree(B))
        max_S_size = floor(Int, α * m)
        for _ in 1:100
            S = rand((n + 1):(n + m), rand(1:max_S_size))
            N_S = Set{Int}()
            for v in S
                for u in Graphs.neighbors(B, v)
                    push!(N_S, u)
                end
            end
            bound = expansion_bound(α, 2, d) * length(S)
            if length(N_S) < bound
                return false
            end
        end
        return true
    end

    @testset "Ramanujan Graph Extended Tests" begin
        # Define a list of valid (p, q) pairs: both p and q are primes, p,q ≡ 1 (mod 4), and p ≠ q.
        test_pairs = [( 5, 29),
                      (13, 17),
                      (17, 13),
                      (29, 13)]
        for (p, q) in test_pairs
            @info "Testing with p = $p, q = $q"
            @test is_prime(p) && p % 4 == 1
            @test is_prime(q) && q % 4 == 1 && p != q
            # Build the finite field
            F = GF(q)
            # Compute the Legendre symbol (p/q)
            symbol = legendre_symbol(p, q)
            if symbol == -1
                # Use GL(2,F) for PGL₂(F)
                GL2 = GL(2, F)
                center = scalar_matrices_GL(GL2)
                PG, morphism = quo(GL2, center)
                # Compute the four-square solutions and process them.
                sols = solve_four_squares(p)
                processed = process_solutions(sols, p)
                @test length(processed) == p + 1
                # Create generators and verify determinants.
                gens = create_generators(processed, F, p)
                for g in gens
                    @test det(g) == F(p)
                end
                # Build the Cayley graph manually.
                gl_gens = [GL2(mat) for mat in gens]
                cayley_g = construct_cayley_graph(PG, gl_gens, morphism)
                # Test regularity: each vertex should have degree p+1.
                for v in Graphs.vertices(cayley_g)
                    @test Graphs.degree(cayley_g, v) == p + 1
                end
                # Test connectivity.
                @test Graphs.is_connected(cayley_g)
                # Expected order for PGL₂(F): q * (q+1) * (q-1)
                expected_order = q * (q + 1) * (q - 1)
                @test nv(cayley_g) == expected_order
                ram_graph = ramanujan_graph(p, q)
                @test is_ramanujan(ram_graph, p) == true
                @test nv(ram_graph) == expected_order
                B = edge_vertex_incidence_graph(ram_graph)
                @test is_unbalanced_bipartite(B)
                @test verify_expansion_property(B, 0.1)
            elseif symbol == 1
                # Use SL(2,F) for PSL₂(F)
                SL2 = SL(2, F)
                center = scalar_matrices_SL(SL2)
                PG, morphism = quo(SL2, center)
                sols = solve_four_squares(p)
                processed = process_solutions(sols, p)
                @test length(processed) == p + 1
                gens = create_generators(processed, F, p)
                for g in gens
                    @test det(g) == F(p)
                end
                # Scale generators so that their determinants become 1.
                s = sqrt(F(p))
                s_inv = inv(s)
                gens_scaled = [s_inv * g for g in gens]
                sl_gens = [SL2(mat) for mat in gens_scaled]
                cayley_g = construct_cayley_graph(PG, sl_gens, morphism)
                for v in Graphs.vertices(cayley_g)
                    @test Graphs.degree(cayley_g, v) == p + 1
                end
                @test Graphs.is_connected(cayley_g)
                # Expected order for PSL₂(F): (q*(q+1)*(q-1)) ÷ 2
                expected_order = (q * (q + 1) * (q - 1)) ÷ 2
                @test nv(cayley_g) == expected_order
                ram_graph = ramanujan_graph(p, q)
                @test is_ramanujan(ram_graph, p) == true
                @test nv(ram_graph) == expected_order
                B = edge_vertex_incidence_graph(ram_graph)
                @test is_unbalanced_bipartite(B)
                @test verify_expansion_property(B, 0.1)
            end
        end
    end

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
