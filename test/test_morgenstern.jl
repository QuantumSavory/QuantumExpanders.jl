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
                group, gens = morgenstern_generators(l, i)
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

    @testset "Non-conjugate and symmetric condition verification" begin
        test_cases = [
            (1, 2), # PSL(2,4)
            (1, 4), # PSL(2,16)
          # (1, 6), # PSL(2,64)
            (2, 2), # PSL(2,16)
          # (3, 2)  # PSL(2,64)
        ]
        for (l, i) in test_cases
            @testset "l=$l, i=$i (q=$(2^l)^$i=$(2^(l*i)))" begin
                group, B = morgenstern_generators(l, i)
                A = alternative_morgenstern_generators(B, FirstOnly())
                Al = alternative_morgenstern_generators(B, AllPairs())
                @test is_symmetric_gen(B)
                @test is_symmetric_gen(A)
                @test is_symmetric_gen(Al)
                @test length(unique(A)) == length(A)
                @test length(unique(Al)) == length(Al)
                @test all(x -> x^2 != one(parent(x[1,1])), A)
                @test all(x -> x^2 != one(parent(x[1,1])), Al)
                @test is_nonconjugate(group, A, B)
                @test is_nonconjugate(group, Al, B)
                @test !is_nonconjugate(group, A, A)
                @test !is_nonconjugate(group, B, B)
                q = 2^(l*i)
                expected_order = q * (q^2 - 1) ÷ gcd(2, q - 1)
                @test order(group) == expected_order
            end
        end
    end
end
