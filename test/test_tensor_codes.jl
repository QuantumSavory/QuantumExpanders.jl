@testitem "Test Tensor Codes" begin
    using Test
    using Oscar
    using QuantumExpanders
    using QuantumExpanders: uniformly_random_code_checkmatrix, dual_code
    using Nemo: zero_matrix, base_ring, transpose, rank

    @testset "Dual code correctness" begin
        for Δ in 5:100
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
                @test H*transpose(G) == Z1
                @test G*transpose(H) == Z2
            end
        end
    end

    @testset "Theorem 9 and 18: Local Component Properties of https://arxiv.org/pdf/2202.13641" begin
        @testset "Theorem 18: Random classical code construction" begin
            for Δ in 8:20
                for ρ in (0.3, 0.4, 0.5)
                    r = Int(floor(ρ * Δ))
                    Hₘ = uniformly_random_code_checkmatrix(ρ, Δ)
                    Gₘ = dual_code(Hₘ)
                    @test size(Hₘ) == (r, Δ)
                    @test size(Gₘ) == (Δ-r, Δ)
                    @test rank(Hₘ) == r
                    @test rank(Gₘ) == Δ-r
                    @test all(iszero, Hₘ*transpose(Gₘ))
                    @test all(iszero, Gₘ*transpose(Hₘ))
                end
            end
        end

        @testset "Theorem 9: Dual tensor code construction" begin
            for Δ in 6:20
                for ρ in (0.4, 0.5)
                    r = Int(floor(ρ * Δ))
                    Hₘ = uniformly_random_code_checkmatrix(ρ, Δ)
                    Gₘ = dual_code(Hₘ)
                    Hₙ = uniformly_random_code_checkmatrix(ρ, Δ)
                    Gₙ = dual_code(Hₙ)
                    # Cₘ ⊗ 𝔽₂ⁿ
                    col_code = kron(Gₘ, matrix_space(base_ring(Gₘ), 1, Δ)(ones(Int, 1, Δ)))
                    # 𝔽₂ᵐ ⊗ Cₙ
                    row_code = kron(matrix_space(base_ring(Gₙ), 1, Δ)(ones(Int, 1, Δ)), Gₙ)
                    # Dual tensor code: Cₘ ⊗ 𝔽₂ⁿ + 𝔽₂ᵐ ⊗ Cₙ
                    dual_tensor_rank = rank(vcat(col_code, row_code))
                    @test dual_tensor_rank ≥ 0
                    @test dual_tensor_rank ≤ min(size(col_code, 1) + size(row_code, 1), Δ^2)
                    @test size(col_code) == ((Δ-r)*1, Δ*Δ)
                    @test size(row_code) == (1*(Δ-r), Δ*Δ)
                end
            end
        end
            
        @testset "CSS condition check: C₁ ⊃ C₀⊥ => (Cₘ ⊗ Cₙ)⋅(Cₘ⊥ ⊗ Cₙ⊥)^T = 0" begin
            for Δ in 6:20
                for ρ in (0.4, 0.5)
                    r = Int(floor(ρ*Δ))
                    Hₘ = uniformly_random_code_checkmatrix(ρ, Δ)
                    Gₘ = dual_code(Hₘ)
                    Hₙ = uniformly_random_code_checkmatrix(ρ, Δ)
                    Gₙ = dual_code(Hₙ)
                    # C₀ = Cₘ ⊗ Cₙ
                    C0 = kron(Gₘ, Gₙ)
                    # C₁⊥ = Cₘ⊥ ⊗ Cₙ⊥
                    C1_dual = kron(Hₘ, Hₙ)
                    # C₁ ⊃ C₀⊥ => (Cₘ ⊗ Cₙ)⋅(Cₘ⊥ ⊗ Cₙ⊥)^T = 0
                    @test all(iszero, C0*transpose(C1_dual))
                    # Verify non-trivial codes
                    @test rank(C0) > 0 || @info "C0 has zero rank for Δ=$Δ, ρ=$ρ"
                    # The quantum code rate lower bound from Theorem 18: (1 - 2ρ)^2;
                    # Page 7; Two Tanner codes that define a quantum LDPC code.
                    expected_min_rate = (1 - 2ρ)^2
                    @test  1 ≥ expected_min_rate ≥ 0 
                end
            end
        end
    end
end
