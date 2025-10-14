@testitem "Test LPS Ramanujan Graph properties" begin
    using Oscar
    using Oscar: center
    using Primes
    using Primes: primes
    using Graphs
    using Random
    using Graphs: degree, vertices, nv, ne, is_bipartite, adjacency_matrix, diameter, is_connected, independent_set, has_edge, MaximalIndependentSet, greedy_color, neighbors, AbstractGraph, degree
    using GraphsColoring: DSATUR, color, Greedy
    using IGraphs: IGraph, IGVectorInt, LibIGraph
    using LinearAlgebra
    using QuantumExpanders
    using QuantumExpanders: legendre_symbol, is_ramanujan
    using LogExpFunctions

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
                @test all(v -> Graphs.degree(cayley_g, v) == p + 1, Graphs.vertices(cayley_g))
                # Test connectivity.
                @test Graphs.is_connected(cayley_g)
                # Expected order for PGL₂(F): q * (q+1) * (q-1)
                expected_order = q*(q+1)*(q-1)
                @test nv(cayley_g) == expected_order
                ram_graph = cayley_g
                @test is_ramanujan(ram_graph, p)
                @test nv(ram_graph) == expected_order
                # When Legendre symbol = -1, the graph must be bipartite [lubotzky1988ramanujan](@cite).
                @test Graphs.is_bipartite(cayley_g)
                # https://igraph.org/c/doc/igraph-Structural.html#igraph_girth
                # Girth Property: g(Xᵖ,ᵗ) ≥ 4 logₚ q - logₚ 4
                # Page 263: Case i (a) of [lubotzky1988ramanujan](@cite)
                girth_Γ_g = floor(Int, 4*log(p, q)-log(p, 4))
                g_igraph = IGraph(cayley_g)
                girth_val = Ref{LibIGraph.igraph_real_t}(0.0)
                cycle = IGVectorInt()
                LibIGraph.igraph_girth(g_igraph.objref, girth_val, cycle.objref)
                actual_girth = isinf(girth_val[]) ? 0 : Int(girth_val[])
                @test actual_girth >= girth_Γ_g
                # Diameter property: diam(Xᵖ,ᵗ) ≤ 2logₚn+2logₚ2+1
                # Page 263: Case i (b) of [lubotzky1988ramanujan](@cite)
                diam = diameter(cayley_g)
                max_diameter = 2*log(p, n)+2*log(p, 2)+1
                @test diam ≤ ceil(Int, max_diameter)
            elseif symbol == 1
                n = q*(q^2-1)/2
                cayley_g = LPS(p, q)
                @test all(v -> Graphs.degree(cayley_g, v) == p + 1, Graphs.vertices(cayley_g))
                @test Graphs.is_connected(cayley_g)
                # Expected order for PSL₂(F): (q*(q+1)*(q-1)) ÷ 2
                expected_order = (q*(q+1)*(q-1))÷2
                @test nv(cayley_g) == expected_order
                ram_graph = cayley_g
                @test is_ramanujan(ram_graph, p)
                # When Legendre symbol = 1, the graph must be non-bipartite [lubotzky1988ramanujan](@cite).
                @test !Graphs.is_bipartite(cayley_g)
                @test nv(ram_graph) == expected_order
                # https://igraph.org/c/doc/igraph-Structural.html#igraph_girth
                # Girth property: g(Xᵖ,ᵗ) ≥ 2 logₚ q
                # Page 263: Case ii (a) of [lubotzky1988ramanujan](@cite)
                girth_Γ_g = floor(Int, 2*log(p, q))
                g_igraph = IGraph(cayley_g)
                girth_val = Ref{LibIGraph.igraph_real_t}(0.0)
                cycle = IGVectorInt()
                LibIGraph.igraph_girth(g_igraph.objref, girth_val, cycle.objref)
                actual_girth = isinf(girth_val[]) ? 0 : Int(girth_val[])
                @test actual_girth >= girth_Γ_g
                # Diameter property: diam(Xᵖ,ᵗ) ≤ 2logₚ n+2logₚ2+1
                # Page 263: Case ii (b) of [lubotzky1988ramanujan](@cite)
                diam = diameter(cayley_g)
                max_diameter = 2*log(p, n)+2*log(p, 2)+1
                @test diam ≤ ceil(Int, max_diameter)
                # Independence number: i(Xᵖ,ᵗ) ≤ (2√p)/(p + 1) n
                # Page 263: Case ii (c) of [lubotzky1988ramanujan](@cite)
                ind_set = independent_set(cayley_g, MaximalIndependentSet())
                independenceₙᵤₘ = length(ind_set)
                independenceₘₐₓₙᵤₘ = ((2*sqrt(p))/(p+1))*n
                @test independenceₙᵤₘ ≤ ceil(Int, independenceₘₐₓₙᵤₘ)
            end
        end
    end

    @testset "Correctness of center computation of SL/GL against Oscar.center" begin
        q_vals = filter(p -> mod(p, 4) == 1, primes(500))
        for q in q_vals
            F = GF(q)
            GL₂ = GL(2, F)
            SL₂ = SL(2, F)
            GL_center = scalar_matrices_GL(GL₂)
            GL_center_groupₒₛ, emb = center(GL₂)
            GL_centerₒₛ = collect(GL_center_groupₒₛ)
            @test length(GL_center) == length(GL_centerₒₛ)
            @test all(M -> emb\M in GL_centerₒₛ, GL_center)
            @test length(GL_center) == order(F) - 1
            SL_center = scalar_matrices_SL(SL₂)
            SL_center_groupₒₛ, emb = center(SL₂)
            SL_centerₒₛ = collect(SL_center_groupₒₛ)
            @test length(SL_center) == length(SL_centerₒₛ)
            @test all(M -> emb\M in SL_centerₒₛ, SL_center)
            expected_sl_size = gcd(2, q - 1)
            @test length(SL_center) == expected_sl_size
        end
    end
end
