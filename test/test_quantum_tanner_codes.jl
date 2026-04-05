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

    @testset "New instances of Quantum Tanner codes using standard dihedral group" begin
        # [[36, 9, (2,1)]]
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
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 2
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # found via random search
        # [[36, 8, (3,3)]]
        H_A = [1  0  1; 0  1  1]
        G_A = [1  1  1]
        @test iszero(mod.(H_A*G_A', 2))
        H_B = [1  1  1]
        G_B = [1  1  0; 1  0  1]
        @test iszero(mod.(H_B*G_B', 2))
        classical_code_pair = ((H_A, G_A), (H_B, G_B))
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 36 == ns && code_k(c) == 8 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 3
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 3

        # [[54, 7, (4,3)]]
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
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 3

        # [[54, 11, (4,1)]]
        H_A = [1 1 1]
        G_A = [1 1 0; 1 0 1]
        @test iszero(mod.(H_A*G_A', 2))
        H_B = [1 1 1; 0 1 1]
        G_B = [0 1 1]
        @test iszero(mod.(H_B*G_B', 2))
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 54 == ns && code_k(c) == 11 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[54, 12, (3,1)]]
        H_A = [1 1 0]
        G_A = [1 1 0; 0 0 1]
        @test iszero(mod.(H_A*G_A', 2))
        H_B = [0 1 0; 0 0 1]
        G_B = [1 0 0]
        @test iszero(mod.(H_B*G_B', 2))
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 54 == ns && code_k(c) == 12 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 3
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[54, 8, (5,1)]]
        H_A = [1 0 1; 1 1 0]
        G_A = [1 1 1]
        @test iszero(mod.(H_A*G_A', 2))
        H_B = [1 1 1; 1 1 0]
        G_B = [1 1 0]
        @test iszero(mod.(H_B*G_B', 2))
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 54 == ns && code_k(c) == 8 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 5
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[72, 16, (4,1)]]
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
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[72, 10, (4,1)]]
        H_A = [1 1 1; 0 0 1]
        G_A = [1 1 0]
        @test iszero(mod.(H_A*G_A', 2))
        H_B = [1 1 1; 1 1 0]
        G_B = [1 1 0]
        @test iszero(mod.(H_B*G_B', 2))
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 72 == ns && code_k(c) == 10 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[72, 12, (4,1)]]
        H_A = [1 0 0; 0 1 1]
        G_A = [0 1 1]
        @test iszero(mod.(H_A*G_A', 2))
        H_B = [1 1 1; 0 1 1]
        G_B = [0 1 1]
        @test iszero(mod.(H_B*G_B', 2))
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab) 
        @test code_n(c) == 72 == ns && code_k(c) == 12 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[72, 14, (4,1)]]
        H_A = [1 1 1]
        G_A = [1 1 0; 1 0 1]
        @test iszero(mod.(H_A*G_A', 2))
        H_B = [1 1 1; 1 0 1]
        G_B = [1 0 1]
        @test iszero(mod.(H_B*G_B', 2))
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 72 == ns && code_k(c) == 14 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[200, 13, (4,4)]]
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
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 4

        # [[200, 19, (7,1)]] 
        H_A = [1  0  1  0  1;
               0  0  0  0  1;
               1  1  0  1  1;
               1  0  0  1  0];
        G_A = [1  0  1  1  0];
        H_B = [1  1  0  1  0;
               1  0  1  0  1;
               0  1  0  0  0]
        G_B = [1  0  1  1  0;
               0  0  1  0  1]
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 200 == ns && code_k(c) == 19 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 7
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[200, 36, (4,1)]]
        H_A = [0  1  0  0  0;
               1  1  1  0  1;
               0  1  0  1  1;
               0  0  1  1  0];
        G_A = [0  0  1  1  1];
        H_B = [0  1  0  1  1;
               1  1  1  0  0]
        G_B = [1  0  1  0  0;
               1  1  0  1  0;
               1  1  0  0  1]
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 200 == ns && code_k(c) == 36 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[200, 27, (6,1)]]
        H_A = [0  1  0  1  0;
               1  0  1  1  0];
        G_A = [1  0  1  0  0;
               1  1  0  1  0;
               0  0  0  0  1];
        H_B = [0  1  0  1  0;
               1  0  1  1  1;
               1  1  0  1  1;
               1  0  0  1  1];
        G_B = [1  0  0  0  1];
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 200 == ns && code_k(c) == 27 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[200, 12, 12]]
        H_A = [1  1  1  1  1;
               0  1  1  1  1;
               0  1  1  0  0;
               0  1  1  0  1]
        G_A = [0  1  1  0  0]
        H_B = [0  1  1  1  0;
               1  1  0  0  0;
               1  1  0  1  1]
        G_B = [1  1  1  0  0;
               1  1  0  1  1]
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 200 == ns && code_k(c) == 12 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 12
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[200, 38, 4]]
        H_A = [0  1  0  1  1;
               1  1  1  0  1]
        G_A = [1  0  1  0  0;
               1  1  0  1  0;
               0  1  0  0  1]
        H_B = [1  0  0  0  0;
               1  1  1  0  1;
               0  1  1  1  0;
               0  1  1  1  1]
        G_B = [0  1  1  0  0]
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 200 == ns && code_k(c) == 38 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[200, 75, 4]]
        H_A = [1  1  1  1  1]
        G_A = [1  1  0  0  0;
               1  0  1  0  0;
               1  0  0  1  0;
               1  0  0  0  1]
        H_B = [1  1  1  0  1;
               1  1  0  0  0;
               0  1  0  0  0;
               1  0  0  1  0]
        G_B = [0  0  1  0  1]
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 200 == ns && code_k(c) == 75 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 1

        # [[250, 14, (6,4)]]
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
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 4

        # [[250, 14, (8,4)]]
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
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z)) == 4

       # [[250, 40, (5,1)]]
        H_A = [0  0  0  1  0;
               0  0  1  0  1;
               0  1  0  0  0;
               0  0  0  1  1]
        G_A = [1  0  0  0  0]
        H_B = [0  1  0  1  0;
               1  0  1  0  0]
        G_B = [1  0  1  0  0;
               0  1  0  1  0;
               0  0  0  0  1]
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 250 == ns && code_k(c) == 40 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 1
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 1

        # [[250, 33, (5,1)]]
        H_A = [1  1  1  1  0;
               1  1  0  1  1;
               0  0  0  1  1;
               1  1  0  1  0];
        G_A = [1  1  0  0  0];
        H_B = [1  1  0  0  1;
               1  0  0  1  0];
        G_B = [0  0  1  0  0;
               1  1  0  1  0;
               0  1  0  0  1];
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 250 == ns && code_k(c) == 33 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 5
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 1
   
       # [[250, 13, (6,1)]]
       H_A =[0  1  1  0  1;
             1  0  0  0  0;
             1  1  1  1  0];
       G_A = [0  1  1  0  0;
              0  1  0  1  1];
       H_B = [1  1  1  0  0;
              1  0  0  1  1]
       G_B = [0  1  1  0  0;
              1  1  0  1  0;
              1  1  0  0  1];
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 250 == ns && code_k(c) == 13 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 1

        # [[250, 93, (3,3)]]
        H_A = [1  1  1  1  1]
        G_A = [1  1  0  0  0;
               1  0  1  0  0;
               1  0  0  1  0;
               1  0  0  0  1]
        H_B = [1  1  0  1  0;
               1  1  0  0  1;
               0  1  1  0  1;
               0  0  1  1  0]
        G_B = [1  0  1  1  1]
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 250 == ns && code_k(c) == 93 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 3
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 3

        # [[250, 93, (4,1)]]
        H_A = [1  0  1  1  1;
               0  0  1  1  1;
               1  0  1  1  0;
               1  0  1  0  0];
        G_A = [0  1  0  0  0];
        H_B = [1  1  1  0  1;
               0  0  1  1  0]
        G_B = [1  1  0  0  0;
               1  0  1  1  0;
               1  0  0  0  1]
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 250 == ns && code_k(c) == 39 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 1

        # [[250, 10, (14,1)]]
        H_A = [1  0  0  0  0;
               1  1  0  1  1];
        G_A = [0  0  1  0  0;
               0  1  0  1  0;
               0  1  0  0  1];
        H_B = [0  1  0  1  0;
               1  1  1  0  1];
        G_B = [1  0  1  0  0;
               1  1  0  1  0;
               1  0  0  0  1]
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 250 == ns && code_k(c) == 10 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 14
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 1

        # [[250, 39, (10,1)]]
        H_A = [1  0  0  0  1;
               1  0  1  0  0;
               1  1  0  0  1;
               0  1  1  0  0]
        G_A = [0  0  0  1  0];
        H_B = [0  1  1  1  1;
               1  1  0  0  1];
        G_B = [1  1  1  0  0;
               1  1  0  1  0;
               0  1  0  0  1];
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 250 == ns && code_k(c) == 39 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 1

        # [[250, 36, 10]]
        H_A = [0  1  0  1  1;
               0  1  0  1  0;
               1  0  0  1  0;
               0  0  0  1  0]
        G_A = [0  0  1  0  0]
        H_B = [1  1  0  0  1;
               1  0  1  1  1]
        G_B = [1  1  1  0  0;
               1  1  0  1  0;
               1  0  0  0  1];
        classical_code_pair = ((H_A, G_A), (H_B, G_B)) # found via random search
        c = QuantumTannerCode(G, A, B, classical_code_pair)
        hx, hz = parity_matrix_x(c), parity_matrix_z(c)
        @test iszero(mod.(hx*hz',2))
        stab = parity_checks(c)
        ns, ks = code_n(stab), code_k(stab)
        @test code_n(c) == 250 == ns && code_k(c) == 36 == ks
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 1

        @testset "New codes found via other non-abelian groups" begin
            # [[216, 30, (4,3)]]
            l = 12
            m = 4
            F = free_group(["r", "s"])
            r, s = gens(F)
            G, = quo(F, [s^m, r^l, s^(-1)*r*s*r])
            r, s = gens(G)
            A = [r^-4*s, s^-1*r^4, s^2]
            B = [r^-3*s, s^-1*r^3, r^6*s^2]
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
            @test code_n(c) == 216 == ns && code_k(c) == 30 == ks
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 4
            @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 3
            @test describe(G) == "C12 : C4"
        end
    end

    @testset "New instances of Quantum Tanner codes using Morgenstern generators" begin
        # [[360, 8, (10,3)]]
        l, i = 1, 2
        q = 2^l
        Δ = q+1
        SL₂, B = morgenstern_generators(l, i)
        A = alternative_morgenstern_generators(B, FirstOnly())
        classical_code_pair = (([0 0 0 1; 1 1 0 0], [1 1 0 0; 0 0 1 0]), ([1 0 1; 0 1 1], [1 1 1])) # found via random search
        c = QuantumTannerCode(SL₂, A, B, classical_code_pair)
        @test stab_looks_good(parity_checks(c), remove_redundant_rows=true)
        @test code_n(c) == 360
        @test code_k(c) == 8
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 3

        # [[360, 10, (4,1)]]
        classical_code_pair = (([0 0 1 1; 0 1 1 1], [1 0 0 0; 0 0 1 1]), ([1 0 1; 0 1 0], [1 0 1])) # found via random search
        c = QuantumTannerCode(SL₂, A, B, classical_code_pair)
        @test stab_looks_good(parity_checks(c), remove_redundant_rows=true)
        @test code_n(c) == 360
        @test code_k(c) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 1

        # [[360, 10, (6,1)]]
        classical_code_pair = (([1 0 1 0; 1 1 0 0], [1 1 1 0; 0 0 0 1]), ([1 1 0; 0 0 1], [1 1 0])) # found via random search
        c = QuantumTannerCode(SL₂, A, B, classical_code_pair)
        @test stab_looks_good(parity_checks(c), remove_redundant_rows=true)
        @test code_n(c) == 360
        @test code_k(c) == 10
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 1

        # [[360, 20, (4,1)]]
        classical_code_pair = (([1 1 0 0; 0 0 1 1], [1 1 0 0; 0 0 1 1]), ([0 0 1; 1 1 0], [1 1 0])) # found via random search
        c = QuantumTannerCode(SL₂, A, B, classical_code_pair)
        @test stab_looks_good(parity_checks(c), remove_redundant_rows=true)
        @test code_n(c) == 360
        @test code_k(c) == 20
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 4
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 1

        # [[360, 22, (6,4)]]
        classical_code_pair = (([1 0 0 1; 0 1 1 0], [0 1 1 0; 1 0 0 1]), ([1 0 1; 1 1 0], [1 1 1])) # found via random search
        c = QuantumTannerCode(SL₂, A, B, classical_code_pair)
        @test stab_looks_good(parity_checks(c), remove_redundant_rows=true)
        @test code_n(c) == 360
        @test code_k(c) == 22
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X, time_limit=900)) == 6
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:Z, time_limit=900)) == 4
    end

    @testset "New instances of Quantum Tanner codes using other Frobenius groups" begin
        # [[150, 48, 4]]
        G = symmetric_group(3)
        rng = MersenneTwister(43)
        S = normal_cayley_subset(G);
        Hᴬ = [1 1 0 1 1]
        Cᴬ = [1 1 0 0 0; 0 0 1 0 0; 1 0 0 1 0; 1 0 0 0 1]
        classical_code_pair = ((Hᴬ, Cᴬ), (Hᴬ, Cᴬ))
        c = GeneralizedQuantumTannerCode(G, S, S, classical_code_pair, bipartite=false, use_same_local_code=true);
        @test stab_looks_good(parity_checks(c), remove_redundant_rows=true)
        @test code_n(c) == 150
        @test code_k(c) == 48
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 4

        # [[36, 16, 3]]
        x = cperm([1,2,3,4])
        G = permutation_group(4, [x])
        rng = MersenneTwister(52)
        S = normal_cayley_subset(G)
        Hᴬ = [1 1 1]
        Cᴬ = [1 1 0; 1 0 1]
        classical_code_pair = ((Hᴬ, Cᴬ), (Hᴬ, Cᴬ))
        c = GeneralizedQuantumTannerCode(G, S, S, classical_code_pair, bipartite=false, use_same_local_code=true);
        @test stab_looks_good(parity_checks(c), remove_redundant_rows=true)
        @test code_n(c) == 36
        @test code_k(c) == 16
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 3

        # [[252, 70, 6]]
        rng = MersenneTwister(54)
        G = cyclic_group(7)
        S = normal_cayley_subset(G)
        Hᴬ = [1 1 1 1 1 1]
        Cᴬ = [1 1 0 0 0 0; 1 0 1 0 0 0; 1 0 0 1 0 0; 1 0 0 0 1 0; 1 0 0 0 0 1]
        classical_code_pair = ((Hᴬ, Cᴬ), (Hᴬ, Cᴬ))
        c = GeneralizedQuantumTannerCode(G, S, S, classical_code_pair, bipartite=false, use_same_local_code=true);
        @test stab_looks_good(parity_checks(c), remove_redundant_rows=true)
        @test code_n(c) == 252
        @test code_k(c) == 70
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 6

       # [[392, 96, 5]]
        rng = MersenneTwister(42)
        G = small_group(8, 4)
        S = normal_cayley_subset(G)
        Hᴬ = [0 1 0 1 1 1 1]
        Cᴬ = [1 0 0 0 0 0 0; 0 0 1 0 0 0 0; 0 1 0 1 0 0 0; 0 1 0 0 1 0 0; 0 1 0 0 0 1 0; 0 1 0 0 0 0 1]
        classical_code_pair = ((Hᴬ, Cᴬ), (Hᴬ, Cᴬ))
        c = GeneralizedQuantumTannerCode(G, S, S, classical_code_pair, bipartite=false, use_same_local_code=true);
        @test stab_looks_good(parity_checks(c), remove_redundant_rows=true)
        @test code_n(c) == 392
        @test code_k(c) == 96
        @test distance(c, DistanceMIPAlgorithm(solver=HiGHS, logical_operator_type=:X)) == 5
    end

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
