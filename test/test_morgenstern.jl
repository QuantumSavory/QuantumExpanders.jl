@testitem "Test Morgensterm generators properties" begin
    using Oscar
    using LinearAlgebra
    using QuantumExpanders
    using QuantumExpanders: morgenstern_solutions, alternative_morgenstern_generators, FirstOnly, AllPairs

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
                _, gens = morgenstern_generators(l, i)
                q = 2^l
                @test length(gens) == q + 1
                Fq = finite_field(2, l, "a")[1]
                R, x = polynomial_ring(Fq, "x")
                ε, sols = morgenstern_solutions(R)
                @test length(sols) == q + 1
                @test all(((γ,δ),) -> γ^2 + γ*δ + ε*δ^2 == one(Fq), sols)
                A_first = alternative_morgenstern_generators(gens, FirstOnly())
                @test length(A_first) == 2*q
                A_pairs = alternative_morgenstern_generators(gens, AllPairs())
                @test length(A_pairs) == q*(q+1)
            end
        end
    end
end
