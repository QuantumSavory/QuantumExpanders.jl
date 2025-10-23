@testitem "Test Quantum Tanner Codes" begin
    using Test
    using Oscar
    using QECCore
    using HiGHS
    using JuMP
    using QuantumExpanders
    using QuantumClifford
    using QuantumClifford.ECC
    using QuantumClifford.ECC: code_n, code_k
    using QuantumExpanders: random_code_pair, parity_matrix
    using Nemo: zero_matrix, base_ring, transpose, rank

    @testset "New instances of Quantum Tanner codes using standard dihedral group" begin
        # [[36, 9, 2]]
        F = free_group([:s, :r])
        s, r = gens(F)
        rels = [s^2, r^4, s*r*s*r]
        G, epimorphism = quo(F, rels)
        s = epimorphism(s)
        r = epimorphism(r)
        A = [s, r, r^3]
        B = [s*r, s*r^3, r^2]
        @test is_symmetric_gen(A)
        @test is_symmetric_gen(B)
        @test is_nonconjugate(G, B, A)
        H_A = [1 0 0; 1 1 1]
        G_A = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_A))))
        H_B = [1 1 1]
        G_B = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_B))))
        classical_code_pair = ((H_A, G_A), (H_B, G_B))
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 36 == ns && code_k(c) == 9 == ks
        distance(c, DistanceMIPAlgorithm(solver=HiGHS))
        classical_code_pair = random_code_pair(0.5, 3)
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        distance(c, DistanceMIPAlgorithm(solver=HiGHS))
    
        # [[54, 7, 4]]
        F = free_group([:s, :r])
        s, r = gens(F)
        rels = [s^2, r^6, s*r*s*r]
        G, epimorphism = quo(F, rels)
        s = epimorphism(s)
        r = epimorphism(r)
        A = [r, r^3, r^5]
        B = [s*r^2, s*r^4, s*r^5]
        @test is_symmetric_gen(A)
        @test is_symmetric_gen(B)
        @test is_nonconjugate(G, B, A)
        H_A = [1 0 0; 1 1 1]
        G_A = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_A))))
        H_B = [1 1 1]
        G_B = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_B))))
        classical_code_pair = ((H_A, G_A), (H_B, G_B))
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 54 == ns && code_k(c) == 7 == ks
        distance(c, DistanceMIPAlgorithm(solver=HiGHS))

        # [[72, 16, 4]]
        F = free_group([:s, :r])
        s, r = gens(F)
        rels = [s^2, r^8, s*r*s*r]
        G, epimorphism = quo(F, rels)
        s = epimorphism(s)
        r = epimorphism(r)
        A = [s, s*r^4, r^4]
        B = [s*r, s*r^3, s*r^7]
        @test is_symmetric_gen(A)
        @test is_symmetric_gen(B)
        @test is_nonconjugate(G, B, A)
        H_A = [1 0 0; 1 1 1]
        G_A = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_A))))
        H_B = [1 1 1]
        G_B = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_B))))
        classical_code_pair = ((H_A, G_A), (H_B, G_B))
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 72 == ns && code_k(c) == 16 == ks
        distance(c, DistanceMIPAlgorithm(solver=HiGHS))

        # [[200, 13, 4]]
        F = free_group([:s, :r])
        s, r = gens(F)
        rels = [s^2, r^8, s*r*s*r]
        G, epimorphism = quo(F, rels)
        s = epimorphism(s)
        r = epimorphism(r)
        A = [s*r^6, r, r^3, r^5, r^7]
        B = [s*r, s*r^3, s*r^7, r^2, r^6]
        @test is_symmetric_gen(A)
        @test is_symmetric_gen(B)
        @test is_nonconjugate(G, B, A)
        H_A = [1 0 1 0 1; 1 1 0 0 0; 1 0 0 0 1]
        G_A = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_A))))
        H_B = [1 1 1 1 1; 0 1 0 0 1]
        G_B = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_B))))
        classical_code_pair = ((H_A, G_A), (H_B, G_B))
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 200 == ns && code_k(c) == 13 == ks
        distance(c, DistanceMIPAlgorithm(solver=HiGHS))

        # [[250, 14, 6]]
        F = free_group([:s, :r])
        s, r = gens(F)
        rels = [s^2, r^10, s*r*s*r]
        G, epimorphism = quo(F, rels)
        s = epimorphism(s)
        r = epimorphism(r)
        A = [s*r, r, r^3, r^7, r^9]
        B = [s*r^6, r^2, r^4, r^6, r^8]
        @test is_symmetric_gen(A)
        @test is_symmetric_gen(B)
        @test is_nonconjugate(G, B, A)
        H_A = [1 1 1 0 1; 1 1 0 0 0; 1 0 0 0 1]
        G_A = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_A))))
        H_B = [1 1 1 0 0; 1 1 0 0 1]
        G_B = Matrix{Int}(lift.(dual_code(matrix(ZZ, H_B))))
        classical_code_pair = ((H_A, G_A), (H_B, G_B))
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 250 == ns && code_k(c) == 14 == ks
        distance(c, DistanceMIPAlgorithm(solver=HiGHS))

        # [[250, 14, 8]]
        A = [r^3, r^-3, r, r^-1, s*r^-2]
        B = [r^4, r^-4, r^2, r^-2, r^5]
        @test is_symmetric_gen(A)
        @test is_symmetric_gen(B)
        @test is_nonconjugate(G, B, A)
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 250 == ns && code_k(c) == 14 == ks
        distance(c, DistanceMIPAlgorithm(solver=HiGHS)) == 8
    end
end
