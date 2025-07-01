@testitem "Test Dual code properties" begin
    using Test
    using QuantumExpanders
    using QuantumExpanders: uniformly_random_code_checkmatrix, dual_code
    using Nemo: zero_matrix, base_ring, transpose, rank

    @testset "Dual code validity" begin
        for Δ in (5, 10, 20)
            for ρ in (0.2, 0.5, 0.8)
                H = uniformly_random_code_checkmatrix(ρ, Δ)
                G = dual_code(H)
                rH = rank(H)
                rG = rank(G)
                @test size(H, 2) == Δ
                @test size(G, 2) == Δ
                @test rH + rG == Δ
                Z1 = zero_matrix(base_ring(H), size(H, 1), size(G, 1))
                Z2 = zero_matrix(base_ring(G), size(G, 1), size(H, 1))
                @test H * transpose(G) == Z1
                @test G * transpose(H) == Z2
            end
        end
    end
end
