@testitem "Test Quantum Tanner Codes" begin
    using Test
    using Oscar
    using QECCore
    using HiGHS
    using JuMP
    using Random
    using Random: seed!
    using QuantumExpanders
    using QuantumClifford
    using QuantumClifford: stab_looks_good, Stabilizer
    using QuantumClifford.ECC
    using QuantumClifford.ECC: code_n, code_k
    using Nemo: zero_matrix, base_ring, transpose, rank

    @testset "Random Symmetric Generating Tests" begin
       for seed in 1:50
           for o in 4:5
               G = symmetric_group(o)
               rng = seed!(seed) 
               A, B = find_random_generating_sets(G, 3; rng=deepcopy(rng))
               H_A = [1 0 1; 1 1 0]
               G_A = [1 1 1]
               H_B = [1 1 1; 1 1 0]
               G_B = [1 1 0]
               classical_code_pair = ((H_A, G_A), (H_B, G_B))
               c = QuantumTannerCode(G, A, B, classical_code_pair)
               @test stab_looks_good(parity_checks(c), remove_redundant_rows=true)
           end
       end
    end

    @testset "Test Puncture" begin
       for seed in 1:50
           G = small_group(12,1)
           rng = MersenneTwister(seed)
           A, B = find_random_generating_sets(G, 6, 5; rng=rng)
           H_A = [1 0 0 0 1 1;
                  0 1 0 1 0 1;
                  0 0 1 1 1 0];
           G_A = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_A))))
           H_B = puncture(H_A, [6])
           G_B = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_B))))
           classical_code_pair = ((Matrix{Int}(H_A), G_A), (H_B, G_B))
           c = QuantumTannerCode(G, A, B, classical_code_pair)
           @test stab_looks_good(parity_checks(c), remove_redundant_rows=true)
       end
       
       # Here is an example of novel [[180, 2, 8]] code
       G = small_group(12,1)
       rng = MersenneTwister(1)
       A, B = find_random_generating_sets(G, 6, 5; rng=rng)
       H_A = [1 0 0 0 1 1;
              0 1 0 1 0 1;
              0 0 1 1 1 0];
       G_A = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_A))))
       H_B = puncture(H_A, [6])
       G_B = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_B))))
       classical_code_pair = ((Matrix{Int}(H_A), G_A), (H_B, G_B))
       c = QuantumTannerCode(G, A, B, classical_code_pair)
       @test stab_looks_good(parity_checks(c), remove_redundant_rows=true) 
    end
end
