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
                SL, B = morgenstern_generators(l, i)
                q = 2^l
                qⁱ = q^i
                @test order(SL) == qⁱ*(qⁱ^2-1)/gcd(2, qⁱ-1)
                @test length(B) == q + 1
                Fq = finite_field(2, l, "a")[1]
                R, x = polynomial_ring(Fq, "x")
                ε, sols = morgenstern_solutions(R)
                @test length(sols) == q + 1
                @test all(((γ,δ),) -> γ^2 + γ*δ + ε*δ^2 == one(Fq), sols)
                @test all(b -> b^2 == one(SL), B)
                A_first = alternative_morgenstern_generators(B, FirstOnly())
                @test length(A_first) == 2*q
                A_pairs = alternative_morgenstern_generators(B, AllPairs())
                @test length(A_pairs) == q*(q+1)
                G = sub(SL, B)[1]
                @test order(G) == order(SL)
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
end
