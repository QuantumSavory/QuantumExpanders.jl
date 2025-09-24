using Test
using Oscar
using LinearAlgebra
import Graphs: degree, is_connected

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
        is_ram = is_ramanujan(g, k)
        @test is_ram == true
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
        is_ram = is_ramanujan(g, k)
        @test is_ram == true
        l = order(G) - k
        go = Float64(order(G))
        trivial_bound = 2 * (sqrt(go) - 1)
        # if  l ≤ 2(sqrt(|G|) - 1) of https://arxiv.org/pdf/1503.04075
        @test Float64(l) <= trivial_bound + 1e-10
    end
end
