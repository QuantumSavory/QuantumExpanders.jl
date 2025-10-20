@testitem "Test Tensor Codes" begin
    using Test
    using Oscar
    using QuantumExpanders
    using QuantumExpanders: random_code, parity_matrix
    using Nemo: zero_matrix, base_ring, transpose, rank

    @testset "Morgenstern Generators" begin
        # Test cases: (l, i) pairs where q=2^l and i is even
        test_cases = [
            (1, 2), # PSL(2,4)
            (1, 4), # PSL(2,16)
            (2, 2), # PSL(2,16)
        ]
        for (l, i) in test_cases
            @testset "l=$l, i=$i (q=$(2^l)^$i=$(2^(l*i)))" begin
                SL₂, B = morgenstern_generators(l, i)
                A_first = alternative_morgenstern_generators(B, FirstOnly())
                A, G = A_first, SL₂
                ρ, Δ = 0.3, 4
                squares = enumerate_square_incidences(G, A, B)
                classical_code_pair = random_code_pair(ρ, Δ)
                hx, hz = parity_matrix(length(G), squares, classical_code_pair)
                iszero(mod.(hx*hz',2))
                A_pairs = alternative_morgenstern_generators(B, AllPairs())
                A, G = A_pairs, SL₂
                squares = enumerate_square_incidences(G, A, B)
                classical_code_pair = random_code_pair(rho, delta)
                hx, hz = parity_matrix(length(G), squares, classical_code_pair)
                iszero(mod.(hx*hz',2))
            end
        end
    end
end
