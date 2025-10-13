@testitem "Test LPS Ramanujan Graph properties" begin
    using Oscar
    using Graphs
    using Graphs: degree, vertices, nv, ne, is_bipartite, adjacency_matrix, diameter, is_connected, independent_set, has_edge, MaximalIndependentSet, greedy_color, neighbors, AbstractGraph, degree
    using GraphsColoring: DSATUR, color, Greedy
    using IGraphs: IGraph, IGVectorInt, LibIGraph
    using LinearAlgebra
    using QuantumExpanders: legendre_symbol
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

    @testset "LPS Ramanujan Graph Tests" begin
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
                n = q*(q^2-1)
                cayley_g = LPS(p, q)
                # Test regularity: each vertex should have degree p+1.
                for v in Graphs.vertices(cayley_g)
                    @test Graphs.degree(cayley_g, v) == p + 1
                end
                # Test connectivity.
                @test Graphs.is_connected(cayley_g)
                # Expected order for PGL₂(F): q * (q+1) * (q-1)
                expected_order = q*(q+1)*(q-1)
                @test nv(cayley_g) == expected_order
                ram_graph = cayley_g
                @test is_ramanujan(ram_graph, p) == true
                @test nv(ram_graph) == expected_order
                B = edge_vertex_incidence_graph(ram_graph)
                @test is_unbalanced_bipartite(B)
                @test verify_expansion_property(B, 0.1)
                @test Graphs.is_bipartite(cayley_g)
                # https://igraph.org/c/doc/igraph-Structural.html#igraph_girth
                girth_Γ_g = floor(Int, 4*log(p, q)-log(p, 4))
                g_igraph = IGraph(cayley_g)
                girth_val = Ref{LibIGraph.igraph_real_t}(0.0)
                cycle = IGVectorInt()
                LibIGraph.igraph_girth(g_igraph.objref, girth_val, cycle.objref)
                actual_girth = isinf(girth_val[]) ? 0 : Int(girth_val[])
                @test actual_girth >= girth_Γ_g 
            elseif symbol == 1
                n = q*(q^2-1)/2
                cayley_g = LPS(p, q)
                for v in Graphs.vertices(cayley_g)
                    @test Graphs.degree(cayley_g, v) == p + 1
                end
                @test Graphs.is_connected(cayley_g)
                # Expected order for PSL₂(F): (q*(q+1)*(q-1)) ÷ 2
                expected_order = (q*(q+1)*(q-1))÷2
                @test nv(cayley_g) == expected_order
                ram_graph = cayley_g
                @test is_ramanujan(ram_graph, p) == true
                @test !Graphs.is_bipartite(cayley_g)
                @test nv(ram_graph) == expected_order
                B = edge_vertex_incidence_graph(ram_graph)
                @test is_unbalanced_bipartite(B)
                @test verify_expansion_property(B, 0.1)
                # https://igraph.org/c/doc/igraph-Structural.html#igraph_girth
                girth_Γ_g = floor(Int, 2*log(p, q))
                g_igraph = IGraph(cayley_g)
                girth_val = Ref{LibIGraph.igraph_real_t}(0.0)
                cycle = IGVectorInt()
                LibIGraph.igraph_girth(g_igraph.objref, girth_val, cycle.objref)
                actual_girth = isinf(girth_val[]) ? 0 : Int(girth_val[])
                @test actual_girth >= girth_Γ_g 
            end
        end
    end
end
