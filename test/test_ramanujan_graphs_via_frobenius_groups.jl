@testitem "Ramanujan Cayley graphs using Frobenius groups" begin
    using Oscar
    using LinearAlgebra
    using Graphs
    using Graphs: degree, is_connected, adjacency_matrix
    using QuantumExpanders
    using QuantumClifford
    using QuantumClifford.ECC

    function _is_ramanujanᵧ(g::SimpleGraph, k::Int)
        A = adjacency_matrix(g)
        λ = real.(eigvals(Matrix(A)))
        sorted_λ = sort(λ, rev=true)
        non_trivial = sorted_λ[2:end]
        bound = 2 * sqrt(k - 1)
        is_ram = all(λ_i -> abs(λ_i) <= bound + 1e-10, non_trivial)
        return is_ram
    end

    @testset "Cayley Graph Construction" begin
        for o in [2, 4, 6]
            for g in [dihedral_group, cyclic_group, symmetric_group]
                G = g(o)
                S = normal_cayley_subset(G)
                g = cayley_right(G, S)
                @test Graphs.nv(g) == order(G)
                @test all(Graphs.degree(g, i) == length(S) for i in 1:Graphs.nv(g))
                @test Graphs.is_connected(g)
                S = normal_cayley_subset(G)
                @test all(s^-1 in S for s in S)
                @test !(one(G) in S)
                @test length(S) > 0
            end
        end
    end

    @testset "Ramanujan graphs using Dihedral Groups" begin
        primes = [3, 5, 7, 11, 13, 17, 19, 23, 29]
        for p in primes
            G = dihedral_group(2*p)
            S = normal_cayley_subset(G)
            g = cayley_right(G, S)
            k = length(S)
            is_ram = _is_ramanujanᵧ(g, k)
            @test is_ram
        end
    end
    
    # Frobenius groups with r ≥ 4 of https://arxiv.org/pdf/1503.04075, see Theorem 3.3, page 6.
    @testset "Frobenius Groups r ≥ 4" begin
        large_primes = [11, 13, 17, 19]
        for p in large_primes
            G = dihedral_group(2*p)
            S = normal_cayley_subset(G)
            g = cayley_right(G, S)
            k = length(S)
            is_ram = _is_ramanujanᵧ(g, k)
            @test is_ram
            l = order(G) - k
            go = Float64(order(G))
            trivial_bound = 2 * (sqrt(go) - 1)
            # if  l ≤ 2(sqrt(|G|) - 1) of https://arxiv.org/pdf/1503.04075
            @test Float64(l) <= trivial_bound + 1e-10
        end
    end

    @testset "Quantum Tanner codes from Frobenius Groups" begin
        primes =  [3, 5, 7, 11, 13, 17, 19, 23, 29]
        for p in primes
            for rate in [0.5, 0.6, 0.7]
                G = dihedral_group(2*p)
                S = normal_cayley_subset(G)
                hx, hz = gen_code(rate, G, S, S, bipartite=false)
                c = Stabilizer(CSS(hx, hz))
                @test stab_looks_good(c, remove_redundant_rows=true)
                hx, hz = gen_good_code(rate, G, S, S, use_same_local_code=true, bipartite=false)
                c = Stabilizer(CSS(hx, hz))
                @test stab_looks_good(c, remove_redundant_rows=true)
            end
        end
    end
end
