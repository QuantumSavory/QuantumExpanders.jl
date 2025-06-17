@testitem "Test Cayley graph properties" begin
    using Oscar
    using LinearAlgebra
    using Graphs
    using Graphs: nv, neighbors, AbstractGraph, degree, ne, add_edge!, SimpleDiGraph, adjacency_matrix
    using QuantumExpanders
    using QuantumExpanders: morgenstern_solutions, morgenstern_f
    using LogExpFunctions
    using LinearAlgebra
    using NautyGraphs
    using Multigraphs
    using Nemo
    import NautyGraphs: NautyGraph, is_isomorphic

    # Construct the Cayleyʳⁱᵍʰᵗ graph for a given group and set of generators.
    function cayley_right_nonlazy(group,generators)
        idx_to_mat = collect(group)
        mat_to_idx = Dict(mat=>i for (i,mat) in pairs(idx_to_mat))
        N = length(group)
        graph = SimpleDiGraph(N)
        for (i,g) in pairs(idx_to_mat)
            for b in generators
                j = mat_to_idx[g*b]
               add_edge!(graph,i,j)
            end
        end
       graph
    end

    # Construct the Cayleyˡᵉᶠᵗ graph for a given group and set of generators.
    function cayley_left_nonlazy(group, generators)
        idx_to_mat = collect(group)
        mat_to_idx = Dict(mat=>i for (i,mat) in pairs(idx_to_mat))
        N = length(group)
        graph = SimpleDiGraph(N)
        for (i,g) in pairs(idx_to_mat)
            for b in generators
                j = mat_to_idx[b*g]
                add_edge!(graph,i,j)
            end
        end
        graph
    end

    # Give all Morgenstern generators over PSL₂qⁱ, where i is even, q=2ˡ, and p is prime.
    function morgenstern_generators_nonlazy(l,i)
        @assert iseven(i)
        p = 2
        q = p^l
        qⁱ = q^i
        𝔽q , unit = finite_field(p,l)
        𝔽qⁱ, punit = finite_field(p,l*i)
        morph = embed(𝔽q,𝔽qⁱ)
        R𝔽q, x = polynomial_ring(𝔽q, "x")
        R𝔽qⁱ, y = polynomial_ring(𝔽qⁱ, "y")
        ε, Bsols = morgenstern_solutions(R𝔽q)
        𝕚s = roots(y^2+y+morph(ε))
        𝕚 = rand(𝕚s)
        SL₂qⁱ = special_linear_group(2,𝔽qⁱ)
        slunit = one(SL₂qⁱ)
        B = typeof(slunit)[]
        for sol in Bsols
            γ,δ = morph.(sol)
            _mat = 𝔽qⁱ[1 γ+δ*𝕚; (γ+δ*𝕚+δ)*punit 1]
            _matp = _mat * inv(sqrt(det(_mat)))
            @assert _mat * inv(sqrt(det(_mat))) == inv(sqrt(det(_mat))) * _mat
            b = SL₂qⁱ(_matp)
            @assert b^2==slunit
            push!(B,b)
        end
        SL₂qⁱ, B
    end

    @testset "Cayley graphs consistency" begin
        test_cases = [
            (1, 2),  # PSL(2,4)
            (1, 4),  # PSL(2,16)
            # (1, 6),# PSL(2,64)
            (2, 2),  # PSL(2,16)
            # (3, 2) # PSL(2,64)
        ]
        for (l, i) in test_cases
            @testset "Cayley Graph Tests: l=$l, i=$i" begin
                group_generators, gorder = morgenstern_generators(l, i)
                gᵣ = cayley_right(group_generators, gorder)
                groupᵣ, group_generatorsᵣ = morgenstern_generators_nonlazy(l, i)
                g1ᵣ = cayley_right_nonlazy(groupᵣ, group_generatorsᵣ)
                mᵣ = Matrix{Int}(adjacency_matrix(gᵣ))
                m1ᵣ = Matrix{Int}(adjacency_matrix(g1ᵣ))
                gnᵣ = NautyGraph(mᵣ)
                gn1ᵣ = NautyGraph(m1ᵣ)
                @test is_isomorphic(gn1ᵣ, gnᵣ)
                @test is_isomorphic(gnᵣ, gn1ᵣ)
                gₗ = cayley_left(group_generators, gorder)
                groupₗ, group_generatorsₗ = morgenstern_generators_nonlazy(l, i)
                g1ₗ = cayley_left_nonlazy(groupₗ, group_generatorsₗ)
                mₗ = Matrix{Int}(adjacency_matrix(gₗ))
                m1ₗ = Matrix{Int}(adjacency_matrix(g1ₗ))
                gnₗ = NautyGraph(mₗ)
                gn1ₗ = NautyGraph(m1ₗ)
                @test is_isomorphic(gn1ₗ, gnₗ)
                @test is_isomorphic(gnₗ, gn1ₗ)
            end
        end

    @testset "Morgenstern generators consistency check: l=$l, i=$i " begin
        # Test cases: (l, i) pairs where q=2^l and i is even
        test_cases = [
            (1, 2), # PSL(2,4)
            (1, 4), # PSL(2,16)
            (1, 6), # PSL(2,64)
            (2, 2), # PSL(2,16)
            (3, 2)  # PSL(2,64)
        ]
        for (l, i) in test_cases
            A, order_A = morgenstern_generators(l, i)
            groupᵣ, gensᵣ = morgenstern_generators_nonlazy(l, i)
            @test order_A == order(groupᵣ)
            F = base_ring(gensᵣ[1])
            G = GL(2, F)
            gens_A = [G(matrix(F, m)) for m in A]
            gens_B = [G(g) for g in gensᵣ]
            H_A = sub(G, gens_A)[1]
            H_B = sub(G, gens_B)[1]
            @test Oscar.is_isomorphic(H_A, H_B)
            SLG = SL(2, F)
            @test first(is_subgroup(H_A, SLG))
            @test first(is_subgroup(H_B, SLG))
            @test exponent(H_A) == exponent(H_B)
            @test length(conjugacy_classes(H_A)) == length(conjugacy_classes(H_B))
            gens_A = sort([G(matrix(F, m)) for m in A], by=x->string(x))
            gens_B = sort([G(g) for g in gensᵣ], by=x->string(x))
            @test gens_A == gens_B
        end
    end
end
