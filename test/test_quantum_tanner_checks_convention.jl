@testitem "Quantum Tanner Checks conventions" begin
    using Test
    using Random: MersenneTwister
    using Oscar
    using QuantumExpanders
    binary_rank(H) = rank(matrix(GF(2), Int.(H)))
    F = free_group([:s, :r])
    s, r = gens(F)
    G, projection = quo(F, [s^2, r^4, s*r*s*r])
    s, r = projection(s), projection(r)
    A = [s, r, r^3]
    B = [s*r, s*r^3, r^2]
    ρ = 0.6
    seed = 64
    @testset "default matches QuantumTannerCode X/Z rank order" begin
        local_rng = MersenneTwister(seed)
        H_A = uniformly_random_code_checkmatrix(1 - ρ,length(A);rng=local_rng,)
        H_B = uniformly_random_code_checkmatrix(ρ,length(B);rng=local_rng,)
        G_A = dual_code(H_A)
        G_B = dual_code(H_B)
        classical_codes = (
            (Matrix{Int}(lift.(H_A)), Matrix{Int}(lift.(G_A)),),
            (Matrix{Int}(lift.(H_B)), Matrix{Int}(lift.(G_B)),),)
        reference = QuantumTannerCode(G,A,B,classical_codes,)
        hx_reference, hz_reference = parity_matrix_xz(reference)
        hx_random, hz_random = random_quantum_Tanner_code(ρ,G,A,B;rng=MersenneTwister(seed),)
        reference_ranks = (binary_rank(hx_reference),binary_rank(hz_reference),)
        random_ranks = (binary_rank(hx_random),binary_rank(hz_random),)
        @test reference_ranks[1] != reference_ranks[2]
        @test random_ranks == reference_ranks
        @test size(hx_random, 2) == size(hx_reference, 2)
        @test size(hz_random, 2) == size(hz_reference, 2)
        @test iszero(mod.(hx_random * hz_random', 2))
    end

    @testset "Gu convention explicitly exchanges X and Z" begin
        for bipartite in (true, false)
            hx_default, hz_default = random_quantum_Tanner_code(ρ,G,A,B;bipartite,rng=MersenneTwister(seed),)
            hx_gu, hz_gu = random_quantum_Tanner_code(ρ,G,A,B;bipartite,stabilizer_convention=:gu,rng=MersenneTwister(seed),)
            @test hx_gu == hz_default
            @test hz_gu == hx_default
            @test iszero(mod.(hx_gu * hz_gu', 2))
        end
        @test_throws ArgumentError random_quantum_Tanner_code(ρ,G,A,B;stabilizer_convention=:unknown,)
    end
end
