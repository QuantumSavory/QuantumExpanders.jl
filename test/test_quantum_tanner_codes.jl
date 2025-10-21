@testitem "Test Quantum Tanner Codes" begin
    using Test
    using Oscar
    using QuantumExpanders
    using QuantumExpanders: random_code_pair, parity_matrix, enumerate_square_incidences
    using Nemo: zero_matrix, base_ring, transpose, rank

    @testset "New instances of Quantum Tanner codes using standard dihedral group" begin
        # [[36, 9, d]]
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
        Q, Q_red = enumerate_squares(G, A, B)
        squares = QuantumExpanders.convert_squares_to_incidence_matrix(Q_red)
        hx, hz = parity_matrix(length(G), squares, classical_code_pair)
        iszero(mod.(hx*hz',2))

        # [[54, 7, d]]
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
        Q, Q_red = enumerate_squares(G, A, B)
        squares = QuantumExpanders.convert_squares_to_incidence_matrix(Q_red)
        hx, hz = parity_matrix(length(G), squares, classical_code_pair)
        iszero(mod.(hx*hz',2))

        # [[72, 16, d]]
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
        Q, Q_red = enumerate_squares(G, A, B)
        squares = QuantumExpanders.convert_squares_to_incidence_matrix(Q_red)
        hx, hz = parity_matrix(length(G), squares, classical_code_pair)
        iszero(mod.(hx*hz',2))

        # [[200, 13, d]]
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
        Q, Q_red = enumerate_squares(G, A, B)
        squares = QuantumExpanders.convert_squares_to_incidence_matrix(Q_red)
        hx, hz = parity_matrix(length(G), squares, classical_code_pair)
        iszero(mod.(hx*hz',2))

        # [[250, 14, d]]
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
        Q, Q_red = enumerate_squares(G, A, B)
        squares = QuantumExpanders.convert_squares_to_incidence_matrix(Q_red)
        hx, hz = parity_matrix(length(G), squares, classical_code_pair)
        iszero(mod.(hx*hz',2))
    end
end
