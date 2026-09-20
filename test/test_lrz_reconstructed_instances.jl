@testitem "Leverrier–Rozendaal–Zémor QT instances" begin
    using Test
    using Oscar
    using GAP
    using QECCore
    using QuantumExpanders
    using QuantumClifford
    using QuantumClifford.ECC
    import Nemo
    using Nemo: matrix, GF

    GAP.Globals.LoadPackage(GAP.GapObj("QDistRnd"))

    const QDIST_TRIALS = 10_000

    function julia_to_gap_gf2(M::AbstractMatrix)
        F2 = GAP.Globals.GF(GAP.Obj(2))
        oneF = GAP.Globals.One(F2)
        zeroF = GAP.Globals.Zero(F2)
        gap_rows = GAP.GapObj([
            GAP.GapObj([isodd(Int(M[i, j])) ? oneF : zeroF for j in axes(M, 2)])
            for i in axes(M, 1)
        ])
        return GAP.Globals.Matrix(gap_rows)
    end

    function compute_qdistrnd_distance(hx, hz; num=QDIST_TRIALS)
        @assert iszero(mod.(hx * hz', 2))
        GX = julia_to_gap_gf2(hx)
        GZ = julia_to_gap_gf2(hz)
        dz = GAP.Globals.DistRandCSS(GX, GZ, GAP.Obj(num), GAP.Obj(0), GAP.Obj(0))
        dx = GAP.Globals.DistRandCSS(GZ, GX, GAP.Obj(num), GAP.Obj(0), GAP.Obj(0))
        return Int(dx), Int(dz)
    end

    # shared shortened Hamming [6,3,3] local code
    H633 = [1 0 0 0 1 1;
            0 1 0 1 0 1;
            0 0 1 1 1 0]
    G633 = [0 1 1 1 0 0;
            1 0 1 0 1 0;
            1 1 0 0 0 1]
    @test iszero(mod.(H633 * G633', 2))
    p633 = [1, 2, 6, 4, 5, 3]

    V4 = codomain(isomorphism(PermGroup, small_group(4, 2)))
    x = cperm(V4, [1, 2]); y = cperm(V4, [3, 4]); e = one(V4)
    @test order(V4) == 4
    A_V4 = [e, e, e, y, x, x*y]

    @testset "LRZ [[144, 8, (≤12, ≤12)]]  (C2 x C2)" begin
        B = [e, e, y, y, x, x*y]
        code = QuantumTannerViaLeftRightActions(V4, A_V4, B, H633, G633, H633, G633; p1=p633, p2=p633)
        hx, hz = parity_matrix_xz(code)
        @test all(sum(hx, dims=2) == 9) && all(sum(hz, dims=2) == 9)
        @test code_n(code) == 144
        @test code_k(code) == 8
        stab = QuantumClifford.ECC.parity_checks(code)
        mat = matrix(GF(2), stab_to_gf2(stab))
        @test rank(mat) == code_n(code) - code_k(code)
        dx, dz = compute_qdistrnd_distance(hx, hz)
        @test min(dx, dz) == 12
    end

    @testset "LRZ [[144, 12, (≤11, ≤11)]]  (C2 x C2)" begin
        B = [e, e, y, y, x, x]
        code = QuantumTannerViaLeftRightActions(V4, A_V4, B, H633, G633, H633, G633; p1=p633, p2=p633)
        hx, hz = parity_matrix_xz(code)
        @test all(sum(hx, dims=2) == 9) && all(sum(hz, dims=2) == 9)
        @test code_n(code) == 144
        @test code_k(code) == 12
        stab = QuantumClifford.ECC.parity_checks(code)
        mat = matrix(GF(2), stab_to_gf2(stab))
        @test rank(mat) == code_n(code) - code_k(code)
        dx, dz = compute_qdistrnd_distance(hx, hz)
        @test min(dx, dz) == 11
    end

    @testset "LRZ [[288, 8, (≤19, ≤19)]]  (C8)" begin
        C8 = codomain(isomorphism(PermGroup, small_group(8, 1)))
        a    = cperm(C8, [1,7,5,3], [2,8,6,4])
        ainv = cperm(C8, [1,3,5,7], [2,4,6,8])
        b    = cperm(C8, [1,8,7,6,5,4,3,2])
        c    = cperm(C8, [1,4,7,2,5,8,3,6])
        d    = cperm(C8, [1,6,3,8,5,2,7,4])
        @test order(C8) == 8
        A = [a, a, ainv, b, c, d]
        B = [one(C8), a, b, b, c, d]
        code = QuantumTannerViaLeftRightActions(C8, A, B, H633, G633, H633, G633; p1=p633, p2=p633)
        hx, hz = parity_matrix_xz(code)
        @test all(sum(hx, dims=2) == 9) && all(sum(hz, dims=2) == 9)
        @test code_n(code) == 288
        @test code_k(code) == 8
        stab = QuantumClifford.ECC.parity_checks(code)
        mat = matrix(GF(2), stab_to_gf2(stab))
        @test rank(mat) == code_n(code) - code_k(code)
        dx, dz = compute_qdistrnd_distance(hx, hz; num=QDIST_TRIALS)
        @test min(dx, dz) == 19
    end
end
