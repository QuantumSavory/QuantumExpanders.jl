using Test
using Oscar
using QuantumExpanders
using QuantumExpanders: morgenstern_generators, morgenstern_solutions, FirstOnly, AllPairs, alternative_morgenstern_generators, is_ramanujan, cayley_right

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
                group_generators, group_order = morgenstern_generators(l, i)
                q = 2^l
                @test length(group_generators) == q + 1
                Fq = finite_field(2, l, "a")[1]
                R, x = polynomial_ring(Fq, "x")
                ε, sols = morgenstern_solutions(R)
                @test length(sols) == q + 1
                @test all(((γ,δ),) -> γ^2 + γ*δ + ε*δ^2 == one(Fq), sols)
                A_first = alternative_morgenstern_generators(group_generators, FirstOnly())
                @test length(A_first) == 2*q
                A_pairs = alternative_morgenstern_generators(group_generators, AllPairs())
                @test length(A_pairs) == q*(q+1)
            end
        end
        test_cases = [
            (1, 2), # PSL(2,4)
            (1, 4), # PSL(2,16)
            (2, 2), # PSL(2,16)
        ]
        for (l, i) in test_cases
            @testset "Ramanujan Graphs via Morgenstern generators: l=$l, i=$i (q=$(2^l)^$i=$(2^(l*i)))" begin
                SL, B = morgenstern_generators(l, i)
                g = QuantumExpanders.cayley_left(SL, B)
                q = 2^l
                @test is_ramanujan(g, q)
                SL, B = morgenstern_generators(l, i)
                g = QuantumExpanders.cayley_right(SL, B)
                @test is_ramanujan(g, q)
            end
        end
    end

    @testset "Non-conjugate condition verification" begin
        test_cases = [
            (1, 2), # PSL(2,4)
            (1, 4), # PSL(2,16)
          # (1, 6), # PSL(2,64)
            (2, 2), # PSL(2,16)
          # (3, 2)  # PSL(2,64)
        ]
        for (l, i) in test_cases
            @testset "l=$l, i=$i (q=$(2^l)^$i=$(2^(l*i)))" begin
                group_generators, group_order = morgenstern_generators(l, i)
                B = group_generators[1:end-1]
                A = alternative_morgenstern_generators(B, FirstOnly())
                Al = alternative_morgenstern_generators(B, AllPairs())
                @test length(unique(A)) == length(A)
                @test length(unique(Al)) == length(Al)
                @test all(x -> x^2 != one(parent(x[1,1])), A)
                @test all(x -> x^2 != one(parent(x[1,1])), Al)
                @test is_nonconjugate(group_generators, group_order, A, B)
                @test is_nonconjugate(group_generators, group_order, Al, B)
                @test !is_nonconjugate(group_generators, group_order, A, A)
                @test !is_nonconjugate(group_generators, group_order, B, B)
                q = 2^(l*i)
                expected_order = q * (q^2 - 1) ÷ gcd(2, q - 1)
                @test group_order == expected_order
            end
        end
    end
end
