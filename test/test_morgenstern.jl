@testitem "Test Morgensterm generators properties" begin
    using Oscar
    using Graphs
    using Graphs: degree, vertices, nv, ne, is_bipartite, adjacency_matrix, diameter, is_connected, independent_set, has_edge, MaximalIndependentSet, greedy_color
    using GraphsColoring: DSATUR, color, Greedy
    using IGraphs: IGraph, IGVectorInt, LibIGraph
    using NautyGraphs: NautyGraph, is_isomorphic
    using LinearAlgebra
    using QECCore
    using QuantumClifford: stab_looks_good, Stabilizer
    using QuantumClifford.ECC
    using QuantumExpanders
    using SimpleGraphConverter
    using SimpleGraphAlgorithms: chromatic_number, UG
    using QuantumExpanders: morgenstern_f, morgenstern_solutions, alternative_morgenstern_generators, FirstOnly, AllPairs, cayley_right

    @testset "Morgenstern Generators" begin
        # Test cases: (l, i) pairs where q=2^l and i is even
        test_cases = [
            (1, 2), # PSL(2,4)
            (1, 4), # PSL(2,16)
            (1, 6), # PSL(2,64)
            (2, 2), # PSL(2,16)
            (3, 2)  # PSL(2,64)
        ]
        for (l, i) in test_cases
            @testset "l=$l, i=$i (q=$(2^l)^$i=$(2^(l*i)))" begin
                SL₂, gens = morgenstern_generators(l, i)
                q = 2^l
                @test length(gens) == q + 1
                Fq = finite_field(2, l, :a)[1]
                R, x = polynomial_ring(Fq, :x)
                ε, sols = morgenstern_solutions(R)
                @test length(sols) == q + 1
                @test all(((γ,δ),) -> γ^2 + γ*δ + ε*δ^2 == one(Fq), sols)
                A_first = alternative_morgenstern_generators(gens, FirstOnly())
                @test length(A_first) == 2*q
                gensₐₗₗ = vcat(gens, A_first)
                H, _ = sub(SL₂, gensₐₗₗ)
                @test H == SL₂
                A_pairs = alternative_morgenstern_generators(gens, AllPairs())
                @test length(A_pairs) == q*(q+1)
                @test is_nonconjugate(SL₂, A_first, gens)
                @test is_nonconjugate(SL₂, A_pairs, gens)
                @test is_symmetric_gen(A_pairs)
                @test is_symmetric_gen(A_first)
                gensₐₗₗ = vcat(gens, A_pairs)
                H, _ = sub(SL₂, gensₐₗₗ)
                @test H == SL₂
            end
        end
    end

    @testset "Morgenstern Ramanujan Graph Properties" begin
        # Test cases: (l, i) pairs where q=2^l and i is even
        test_cases = [
            (1, 2), # PSL(2,4)
            (1, 4), # PSL(2,16)
            # (1, 6), # PSL(2,64) # take long time
            (2, 2), # PSL(2,16)
            # (3, 2) # PSL(2,64) # take long time
        ]
        for (l, i) in test_cases
            @testset "l=$l, i=$i (q=$(2^l), |Γ|=$(2^(3l*i) - 2^(l*i)))" begin
                q = 2^l
                r = q + 1
                expected_order = q^(3*i)-q^i
                G, gens = morgenstern_generators(l, i)
                graph = cayley_right(G, gens)
                # Property I: Regularity
                degrees = [degree(graph, v) for v in vertices(graph)]
                @test all(deg == q+1 for deg in degrees)
                @test maximum(degrees) == minimum(degrees)
                # Order
                @test nv(graph) == expected_order
                # Property II: Non-bipartiteness
                @test !is_bipartite(graph)
                # Property III: Ramanujan bound
                A = Matrix(adjacency_matrix(graph))
                eigenvalues = eigvals(A)
                sorted_evals = sort(real.(eigenvalues), rev=true)
                # The largest eigenvalue should be q+1.
                @test isapprox(sorted_evals[1], q+1, atol=1e-10)
                # All other eigenvalues should satisfy |μ| ≤ 2√q
                ramanujan_bound = 2*sqrt(q)
                @test all(abs(μ) ≤ ramanujan_bound+1e-10 for μ in sorted_evals[2:end])
                violating_eigenvalues = count(μ -> abs(μ) > ramanujan_bound+1e-10, sorted_evals[2:end])
                @test violating_eigenvalues == 0
                # Property III of Theorem 5.13: g(Γ) ≥ (2/3)⋅log_q(|Γ|)
                # https://igraph.org/c/doc/igraph-Structural.html#igraph_girth
                girth_Γ_g = floor(Int, (2/3)*log(q, expected_order))
                g_igraph = IGraph(graph)
                girth_val = Ref{LibIGraph.igraph_real_t}(0.0)
                cycle = IGVectorInt()
                LibIGraph.igraph_girth(g_igraph.objref, girth_val, cycle.objref)
                actual_girth = isinf(girth_val[]) ? 0 : Int(girth_val[])
                @test actual_girth >= girth_Γ_g 
                # Property IV: Diameter bound
                diam = diameter(graph)
                max_diameter = 2*log(q, expected_order)+2
                @test diam ≤ ceil(Int, max_diameter)
                # Property V: The chromatic number property
                # Theorem 5.13: χ(Γ_g) ≥ ((q+1)/2√q)+1
                coloring = greedy_color(graph; sort_degree=false, reps=1000)
                χ_greedy = coloring.num_colors
                @test χ_greedy >= (q+1)/(2*sqrt(q))+1
                colors_dsatur = color(graph; algorithm=DSATUR())
                χ_dsatur = length(colors_dsatur.colors) 
                @test χ_dsatur >= (q+1)/(2*sqrt(q))+1
                colors_greedy = color(graph; algorithm=Greedy())
                χ_greedy₂ = length(colors_greedy.colors) 
                @test χ_greedy₂ >= (q+1)/(2*sqrt(q))+1
                # Properties of generator set B
                @test length(gens) == q+1
                # All generators should have determinant 1.
                @test all(det(gen) == one(base_ring(gen)) for gen in gens)
                # All generators should have order 2.
                @test all(matrix(gen^2) == identity_matrix(base_ring(gen), 2) for gen in gens)
                @test is_connected(graph)
                # Property VI: Independence number
                ind_set = independent_set(graph, MaximalIndependentSet())
                independenceₙᵤₘ = length(ind_set)
                # Theorem: i(Γ_g) ≤ (2√q/(q+1)) |Γ_g|
                independenceₘₐₓₙᵤₘ = (2*sqrt(q)/(q+1))*nv(graph)
                @test independenceₙᵤₘ ≤ ceil(Int, independenceₘₐₓₙᵤₘ)
                @test all(u == v || !has_edge(graph, u, v) for u in ind_set, v in ind_set)
                # Additional expander properties from Theorem [A1, M1]
                n = nv(graph)
                λ = maximum(abs.(sorted_evals[2:end]))
                # Γ is an (n, r, 1 - λ²/r²)-expander
                d_expander = 1-(λ^2)/(r^2)
                @test d_expander > 0
                # λ ≤ r - d²/8r with d = d_expander
                @test λ ≤ r-(d_expander^2)/(8*r)+1e-10
                optimal_bound = 2*sqrt(r-1)
                @test λ ≤ optimal_bound+1e-10
            end
        end
    end

    @testset "Morgenstern Spectral Expansion Bounds" begin
        # Test cases: (l, i) pairs where q=2^l and i is even
        test_cases = [
            (1, 2), # PSL(2,4)
            (1, 4), # PSL(2,16)
            # (1, 6), # PSL(2,64) # take long time
            (2, 2), # PSL(2,16)
            # (3, 2) # PSL(2,64) # take long time
        ]
        for (l, i) in test_cases
            @testset "l=$l, i=$i (q=$(2^l))" begin
                q = 2^l
                G, B = morgenstern_generators(l, i)
                @testset "Alternative AllPairs Generators" begin
                    A = alternative_morgenstern_generators(B, AllPairs())
                    graph = cayley_right(G, A)
                    adj_mat = Matrix(adjacency_matrix(graph))
                    # Normalize by degree: k₁ = q²+q
                    adj_matₙₒᵣₘ = adj_mat/(q^2+q)
                    eigenvals = sort(real.(eigvals(adj_matₙₒᵣₘ)), rev=true)
                    λ = eigenvals[2]
                    tbound = (3q-1)/(q^2+q) # Claim 6.1 (ii) of [dinur2022locally](@cite)
                    @test λ < tbound
                end

                @testset "Alternative FirstOnly Generators" begin
                    A = alternative_morgenstern_generators(B, FirstOnly())
                    graph = cayley_right(G, A)
                    adj_matrix = Matrix(adjacency_matrix(graph))
                    # Normalize by degree: k₁ = 2q
                    adj_matₙₒᵣₘ = adj_matrix/2q
                    eigenvals = sort(real.(eigvals(adj_matₙₒᵣₘ)), rev=true)
                    λ = eigenvals[2]
                    tbound = (3*sqrt(2q-1))/(2q) # Claim 6.2 of [dinur2022locally](@cite)
                    @test λ < tbound
                end
            end
        end
    end

    function morgenstern_solutions_slow(R::FqPolyRing, ε)
        F = base_ring(R)
        sols = [(one(F), zero(F))]
        for δ in F
            iszero(δ) && continue
            for γ in F
                if γ^2+γ*δ+δ^2*ε == one(F)
                    push!(sols, (γ, δ))
                    other_γ = γ+δ
                    push!(sols, (other_γ, δ))
                    break
                end
            end
        end
        return unique(sols)
    end

    function morgenstern_solutions_fast(R::FqPolyRing, ε)
        F = base_ring(R)
        sols = [(one(F), zero(F))]
        x = gen(R)
        f = x^2+x+ε
        for s in F
            fs = f(s)
            sfs = sqrt(fs)
            α = s*inv(sfs)
            β = inv(sfs)
            push!(sols, (α, β))
        end
        return sols
    end

    @testset "Morgenstern solutions correctness" begin
        for exp in [4, 6, 8, 10]
            F = GF(2, exp)
            R, x = polynomial_ring(F, :x)
            f = morgenstern_f(R)
            ε = coeff(f, 0)
            sols_fast = morgenstern_solutions_fast(R, ε)
            sols_slow = morgenstern_solutions_slow(R, ε)
            @test Set(sols_fast) == Set(sols_slow)
            @test length(sols_fast) == length(sols_slow) == order(F) + 1
            @test all(((γ, δ),) -> γ^2+γ*δ+δ^2*ε == one(F), sols_fast)
            @test all(((γ, δ),) -> γ^2+γ*δ+δ^2*ε == one(F), sols_slow)
        end
    end

    @testset "Cayley Graph Isomorphism: Remark 3.2 of [dinur2022locally](@cite)" begin
        test_cases = [
            (1, 2), # PSL(2,4)
            (1, 4), # PSL(2,16)
            (2, 2)  # PSL(2,16)
        ]
        for (l, i) in test_cases
            @testset "l=$l, i=$i (q=$(2^l)^$i=$(2^(l*i)))" begin
                SL₂, gens = morgenstern_generators(l, i)
                q = 2^l
                cayleyᵣ = cayley_right(SL₂, gens)
                cayleyₗ = cayley_left(SL₂, gens)
                cgᵣ = NautyGraph(cayleyᵣ)
                cgₗ = NautyGraph(cayleyₗ)
                @test is_isomorphic(cgᵣ, cgₗ)
                A_first = alternative_morgenstern_generators(gens, FirstOnly())
                cayleyᵣ = cayley_right(SL₂, A_first)
                cayleyₗ = cayley_left(SL₂, A_first)
                cgᵣ = NautyGraph(cayleyᵣ)
                cgₗ = NautyGraph(cayleyₗ)
                @test is_isomorphic(cgᵣ, cgₗ)
                A_pairs = alternative_morgenstern_generators(gens, AllPairs())
                cayleyᵣ = cayley_right(SL₂, A_pairs)
                cayleyₗ = cayley_left(SL₂, A_pairs)
                cgᵣ = NautyGraph(cayleyᵣ)
                cgₗ = NautyGraph(cayleyₗ)
                @test is_isomorphic(cgᵣ, cgₗ)
            end
        end
    end

    @testset "Quantum Tanner codes based on Morgenstern Generators" begin
        test_cases = [
            (1, 2), # PSL(2,4)
            # (1, 4) # PSL(2,16) takes long time
        ]
        for (l, i) in test_cases
            @testset "l=$l, i=$i (q=$(2^l)^$i=$(2^(l*i)))" begin
                q = 2^l
                Δ = q+1
                SL₂, B = morgenstern_generators(l, i)
                A = alternative_morgenstern_generators(B, FirstOnly())
                for _ in 1:5
                    for rate in [0.4, 0.5, 0.6, 0.7]
                        hx, hz = gen_code(rate, SL₂, A, B) # Q = 360
                        c = Stabilizer(CSS(hx, hz))
                        @test stab_looks_good(c, remove_redundant_rows=true)
                        hx, hz = gen_good_code(rate, SL₂, A, B)
                        c = Stabilizer(CSS(hx, hz))
                        @test stab_looks_good(c, remove_redundant_rows=true)
                    end
                end
            end
            @testset "l=$l, i=$i (q=$(2^l)^$i=$(2^(l*i)))" begin
                q = 2^l
                Δ = q+1
                SL₂, B = morgenstern_generators(l, i)
                A = alternative_morgenstern_generators(B, FirstOnly())
                for _ in 1:5
                    for rate in [0.5, 0.6, 0.7]
                        hx, hz = gen_code(rate, SL₂, A, A, bipartite=false) # Q = 960
                        c = Stabilizer(CSS(hx, hz))
                        @test stab_looks_good(c, remove_redundant_rows=true)
                        hx, hz = gen_code(0.4, SL₂, B, B, use_same_local_code=true, bipartite=false)
                        c = Stabilizer(CSS(hx, hz))
                        @test stab_looks_good(c, remove_redundant_rows=true)
                        hx, hz = gen_good_code(rate, SL₂, A, B)
                        c = Stabilizer(CSS(hx, hz))
                        @test stab_looks_good(c, remove_redundant_rows=true)
                    end
                end
            end
            @testset "l=$l, i=$i (q=$(2^l)^$i=$(2^(l*i)))" begin
                q = 2^l
                Δ = q+1
                SL₂, B = morgenstern_generators(l, i)
                A = alternative_morgenstern_generators(B, FirstOnly())
                for _ in 1:5
                    hx, hz = gen_code(0.4, SL₂, B, B, bipartite=false) # Q = 540
                    c = Stabilizer(CSS(hx, hz))
                    @test stab_looks_good(c, remove_redundant_rows=true)
                    c = Stabilizer(CSS(hx, hz))
                    hx, hz = gen_code(0.4, SL₂, B, B, use_same_local_code=true, bipartite=false)
                    @test stab_looks_good(c, remove_redundant_rows=true)
                    hx, hz = gen_good_code(rate, SL₂, A, B)
                    c = Stabilizer(CSS(hx, hz))
                    @test stab_looks_good(c, remove_redundant_rows=true)
                end
            end
            @testset "l=$l, i=$i (q=$(2^l)^$i=$(2^(l*i)))" begin
                q = 2^l
                Δ = q+1
                SL₂, B = morgenstern_generators(l, i)
                A = alternative_morgenstern_generators(B, FirstOnly())
                for _ in 1:5
                    for rate in [0.4, 0.5, 0.6, 0.7]
                        hx, hz = gen_code(rate, SL₂, A, B, bipartite=false) # Q = 720
                        c = Stabilizer(CSS(hx, hz))
                        @test stab_looks_good(c, remove_redundant_rows=true)
                        hx, hz = gen_good_code(rate, SL₂, A, B)
                        c = Stabilizer(CSS(hx, hz))
                        @test stab_looks_good(c, remove_redundant_rows=true)
                    end
                end
            end
            @testset "l=$l, i=$i (q=$(2^l)^$i=$(2^(l*i)))" begin
                q = 2^l
                Δ = q+1
                SL₂, B = morgenstern_generators(l, i)
                A = alternative_morgenstern_generators(B, FirstOnly())
                for _ in 1:5
                    for rate in [0.4, 0.5, 0.6]
                        hx, hz = gen_code(rate, SL₂, B, A, bipartite=false) # Q = 720
                        c = Stabilizer(CSS(hx, hz))
                        @test stab_looks_good(c, remove_redundant_rows=true)
                        hx, hz = gen_good_code(rate, SL₂, A, B)
                        c = Stabilizer(CSS(hx, hz))
                        @test stab_looks_good(c, remove_redundant_rows=true)
                    end
                end
            end
        end
    end
end
