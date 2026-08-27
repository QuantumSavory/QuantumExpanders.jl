@testitem "Morgenstern QT codes from appendix" begin
    using Test
    using Oscar
    using GAP
    using QECCore
    using QuantumExpanders
    using QuantumClifford
    using QuantumClifford.ECC
    import Nemo
    using Nemo: matrix, GF

    loaded = GAP.Globals.LoadPackage(GAP.GapObj("QDistRnd"))

    const QDIST_TRIALS2 = 50_000

    function dual_code_int(H)
        H_nemo = matrix(GF(2), H)
        null = Nemo.nullspace(H_nemo)[2]
        @assert all(iszero, H_nemo * null)
        null_t = Nemo.transpose(null)
        return Matrix{Int}([
            Int(lift(ZZ, null_t[i, j]))
            for i in 1:nrows(null_t), j in 1:ncols(null_t)
        ])
    end

    function julia_to_gap_gf2(M::AbstractMatrix)
        F2 = GAP.Globals.GF(GAP.Obj(2))
        oneF = GAP.Globals.One(F2)
        zeroF = GAP.Globals.Zero(F2)
        gap_rows = GAP.GapObj([
            GAP.GapObj([
                isodd(Int(M[i, j])) ? oneF : zeroF
                for j in axes(M, 2)
            ])
            for i in axes(M, 1)
        ])
        return GAP.Globals.Matrix(gap_rows)
    end

    function compute_qdistrnd_distance(hx, hz; num=QDIST_TRIALS)
        @assert iszero(mod.(hx * hz', 2))
        GX = julia_to_gap_gf2(hx)
        GZ = julia_to_gap_gf2(hz)
        dz = GAP.Globals.DistRandCSS(
            GX,
            GZ,
            GAP.Obj(num),
            GAP.Obj(0),
            GAP.Obj(0),
        )
        dx = GAP.Globals.DistRandCSS(
            GZ,
            GX,
            GAP.Obj(num),
            GAP.Obj(0),
            GAP.Obj(0),
        )
        return Int(dx), Int(dz)
    end

    l = 1
    i = 2
    SL2, B = morgenstern_generators(l, i)
    A = alternative_morgenstern_generators(B, FirstOnly())
    @test order(SL2) == 60
    @test length(A) == 4
    @test length(B) == 3

    @testset "Morgenstern [[360, 8, (≤ 20, ≤ 3)]]" begin
        H_A = [
            1 1 0 1;
            1 1 1 1
        ]
        H_B = [1 1 1]
        G_A = dual_code_int(H_A)
        G_B = dual_code_int(H_B)
        c = QuantumTannerCode(
            SL2,
            A,
            B,
            ((H_A, G_A), (H_B, G_B)),
        )
        hx, hz = parity_matrix_xz(c)
        @test code_n(c) == 360
        @test code_k(c) == 8
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        @test rank(mat) == code_n(c) - code_k(c)
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS2,
        )
        println()
        println("Morgenstern [[360, 8, (≤ 20, ≤ 3)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  paper table  : (dx, dz) = (≤20, ≤3)")
    end

    @testset "Morgenstern [[360, 61, (≤ 10, ≤ 3)]]" begin
        H_A = [1 1 1 0]
        H_B = [
            0 1 1;
            1 1 0
        ]
        G_A = dual_code_int(H_A)
        G_B = dual_code_int(H_B)
        c = QuantumTannerCode(
            SL2,
            A,
            B,
            ((H_A, G_A), (H_B, G_B)),
        )
        hx, hz = parity_matrix_xz(c)
        @test code_n(c) == 360
        @test code_k(c) == 61
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        @test rank(mat) == code_n(c) - code_k(c)
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println()
        println("Morgenstern [[360, 61, (≤ 10, ≤ 3)]]")
        println("  QDistRnd 50K          : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  paper table           : (dx, dz) = (≤10, ≤3)")
        println("  appendix example 50K  : (dx, dz) = (10, 3)")
    end
end
