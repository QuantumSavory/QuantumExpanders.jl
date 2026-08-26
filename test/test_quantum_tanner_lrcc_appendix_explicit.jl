@testitem "Quantum Tanner Codes. All LRCC appendix constructions" begin
    using Test
    using Oscar
    using QECCore
    using QuantumExpanders
    using QuantumClifford
    using QuantumClifford.ECC
    using GAP

    H633 = [
        1 0 0 0 1 1;
        0 1 0 1 0 1;
        0 0 1 1 1 0
    ]
    G633 = Matrix{Int}(lift.(dual_code(H633)))

    H844 = [
        1 0 0 0 0 1 1 1;
        0 1 0 0 1 0 1 1;
        0 0 1 0 1 1 0 1;
        0 0 0 1 1 1 1 0
    ]
    G844 = Matrix{Int}(lift.(dual_code(H844)))

    G734 = [
        1  0  1  1  1  0  0;
        1  1  1  0  0  1  0;
        0  1  1  1  0  0  1
    ];

    H734 = Matrix{Int}(lift.(dual_code(G734)));

    G743 = [
        1  0  0  0  1  1  0;
        0  1  0  0  1  0  1;
        0  0  1  0  0  1  1;
        0  0  0  1  1  1  1
        ];

    H743 = Matrix{Int}(lift.(dual_code(G743)))


    loaded = GAP.Globals.LoadPackage(GAP.GapObj("QDistRnd"))

    const QDIST_TRIALS = 20_000

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

    @testset "C3 x S3 [[324, 8, (≤ 17, ≤ 14)]]" begin
        G = small_group(18, 3) # C3 x S3
        r, s, t = Oscar.gens(G)
        A = [s, s^2, t, t^2, r*t^2, r]
        B = [r*s, r*s^2, s*t, s^2*t^2, r*s*t^2, r*s^2*t^2]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H633, G633)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 324
        @test code_k(c) == 8
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 9
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 x S3 [[324, 8, (≤ 17, ≤ 14)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤17, ≤14)")

    end

    @testset "C3 x S3 [[576, 22, (≤ 16, ≤ 16)]]" begin
        G = small_group(18, 3) # C3 x S3
        r, s, t = Oscar.gens(G)
        A = [r*s*t^2, r*s^2*t^2, r*s^2*t, r*s*t, t^2, t, r*s, r*s^2]
        B = [r*t^2, s^2, s, s*t, s^2*t^2, s^2*t, s*t^2, r]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 22
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 x S3 [[576, 22, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C2 x A4 [[432, 16, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 13) # C2 x A4
        r, s, t, u = Oscar.gens(G)
        A = [s^2*t*u, s*t, t, s^2*u, s*t*u, r*t]
        B = [r*s^2*t*u, r*s*t, r*s^2*u, r*s*t*u, r*s^2, r*s]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H633, G633)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 432
        @test code_k(c) == 16
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 9
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x A4 [[432, 16, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "S4 [[768, 32, (≤ 16, ≤ 16)]]" begin
        G = small_group(24, 12) # S4
        r, s, t, u = Oscar.gens(G)
        A = [s^2*t, s*u, s^2, s, s*t*u, s^2*u, s*t, s^2*t*u]
        B = [r*s^2*u, r*s^2*t*u, r*t, r*u, r*s^2, r*s*t, r*s*t*u, r]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 768
        @test code_k(c) == 32
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("S4 [[768, 32, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C2 x A4 [[768, 42, (≤ 16, ≤ 16)]]" begin
        G = small_group(24, 13) # C2 x A4
        r, s, t, u = Oscar.gens(G)
        A = [r, r*s^2*t*u, r*s*t, r*s^2*t, r*s*u, r*s^2*u, r*s*t*u, t]
        B = [s, s^2, s*t, s^2*t*u, s*t*u, s^2*u, r*u, r*t*u]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 768
        @test code_k(c) == 42
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x A4 [[768, 42, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C2 x A4 [[768, 20, (≤ 16, ≤ 16)]]" begin
        G = small_group(24, 13) # C2 x A4
        r, s, t, u = Oscar.gens(G)
        A = [s, s^2, s^2*t*u, s*t, r*t, s*t*u, s^2*u, r]
        B = [t, r*s*t*u, r*s^2*u, r*s, r*s^2, r*s*t, r*s^2*t*u, t*u]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 768
        @test code_k(c) == 20
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x A4 [[768, 20, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C12 x C2 [[768, 40, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 9) # C12 x C2
        r, s, t, u = Oscar.gens(G)
        A = [s*u, r*t^2, r*t*u, r*s*t^2, r*s*t*u, t^2*u, t*u, u]
        B = [s*t^2*u, s*t*u, r*s*u, r*s, t^2, t, r*u, r]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 768
        @test code_k(c) == 40
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C12 x C2 [[768, 40, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "(C6 x C2) : C2 [[768, 38, (≤ 16, ≤ 16)]]" begin
        G = small_group(24, 8) # (C6 x C2) : C2
        r, s, t, u = Oscar.gens(G)
        A = [r*u, u^2, u, s*t, r*u^2, t*u, t*u^2, r*t]
        B = [r*s*t*u^2, r*s*u^2, r*s*u, r*s*t*u, r*s*t, r*s, s*t*u^2, s*t*u]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 768
        @test code_k(c) == 38
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C6 x C2) : C2 [[768, 38, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C2 x (C3 : C4) [[768, 34, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 7) # C2 x (C3 : C4)
        r, s, t, u = Oscar.gens(G)
        A = [r*u^2, r*t*u^2, t*u, t*u^2, s*t*u^2, s*t*u, s*u, s*u^2]
        B = [u^2, u, r*s*u^2, r*s*t*u^2, t, s, r*s*t, r*s]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 768
        @test code_k(c) == 34
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x (C3 : C4) [[768, 34, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C4 x S3 [[768, 48, (≤ 16, ≤ 16)]]" begin
        G = small_group(24, 5) # C4 x S3
        r, s, t, u = Oscar.gens(G)
        A = [r*s*t*u^2, r*s*u^2, r*s*t, r*s, s, s*t, r*s*t*u, r*s*u]
        B = [r*t, t*u^2, t*u, s*u^2, s*t*u, r*u^2, t, r*t*u^2]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 768
        @test code_k(c) == 48
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C4 x S3 [[768, 48, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C3 : Q8 [[768, 28, (≤ 16, ≤ 16)]]" begin
        G = small_group(24, 4) # C3 : Q8
        r, s, t, u = Oscar.gens(G)
        A = [r*s*t*u^2, r*s*u^2, r*s*t, r*s, s, s*t, r*s*t*u, r*s*u]
        B = [r*t, r, t*u^2, t*u, s*u^2, s*t*u, r*u^2, r*t*u^2]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 768
        @test code_k(c) == 28
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 : Q8 [[768, 28, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "SL(2,3) [[768, 32, (≤ 16, ≤ 16)]]" begin
        G = small_group(24, 3) # SL(2,3)
        r, s, t, u = Oscar.gens(G)
        A = [r^2*s*t*u, r*s, r^2, r, r^2*s*u, r*t, s, s*u]
        B = [r*s*u, r^2*s*t, r*s*t*u, r^2*t, r*u, r^2*u, r*t*u, r^2*s]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 768
        @test code_k(c) == 32
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("SL(2,3) [[768, 32, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C3 : C8 [[768, 32, (≤ 16, ≤ 16)]]" begin
        G = small_group(24, 1) # C3 : C8
        r, s, t, u = Oscar.gens(G)
        A = [r*u, r*s*t*u, r*u^2, r*s*t*u^2, s, s*t, t*u^2, t*u]
        B = [r*s, r*t, s*u^2, s*t*u, r*s*u^2, r*t*u^2, u, u^2]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 768
        @test code_k(c) == 32
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 : C8 [[768, 32, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C3 : C8 [[768, 24, (≤ 16, ≤ 16)]]" begin
        G = small_group(24, 1) # C3 : C8
        r, s, t, u = Oscar.gens(G)
        A = [r*s*t, r, s*t*u, s*u^2, r*s*t*u, r*u, s*t*u^2, s*u]
        B = [r*s, r*t, t*u^2, t*u, s, s*t, r*t*u, r*s*u]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 768
        @test code_k(c) == 24
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 : C8 [[768, 24, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C5 x C5 [[800, 24, (≤ 20, ≤ 20)]]" begin
        G = small_group(25, 2) # C5 x C5
        r, s = Oscar.gens(G)
        A = [r*s^2, r^4*s^3, r^3*s, r^2*s^4, r^4*s^2, r*s^3, r, r^4]
        B = [r^3, r^2, r^4*s, r*s^4, s, s^4, r^2*s, r^3*s^4]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 800
        @test code_k(c) == 24
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16 
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C5 x C5 [[800, 24, (≤ 20, ≤ 20)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤20, ≤20)")
    end

    @testset "C9 x C3 [[864, 32, (≤ 12, ≤ 12)]]" begin
        G = small_group(27, 2) # C9 x C3
        r, s, t = Oscar.gens(G)
        A = [r^2*t^2, r, t, t^2, r^2, r*t^2, r*s^2, r^2*s*t^2]
        B = [r^2*s^2, r*s*t^2, r^2*s*t, r*s^2*t, s^2*t, s*t^2, r^2*s, r*s^2*t^2]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 864
        @test code_k(c) == 32
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C9 x C3 [[864, 32, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "(C3 x C3) : C3 [[864, 20, (≤ 24, ≤ 24)]]" begin
        G = small_group(27, 3) # (C3 x C3) : C3
        r, s, t = Oscar.gens(G)
        A = [r^2*t, r*t^2, t^2, t, r*t, r^2*t^2, r, r^2]
        B = [r^2*s*t, r*s^2*t, r*s*t^2, r^2*s^2*t^2, r^2*s^2*t, r*s, s*t^2, s^2*t]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 864
        @test code_k(c) == 20
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C3 x C3) : C3 [[864, 20, (≤ 24, ≤ 24)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤24, ≤24)")
    end

    @testset "(C3 x C3) : C3 [[864, 16, (≤ 24, ≤ 24)]]" begin
        G = small_group(27, 3) # (C3 x C3) : C3
        r, s, t = Oscar.gens(G)
        A = [r^2*s, r*s^2*t^2, r^2*s*t, r*s^2*t, s^2, s, s*t, s^2*t^2]
        B = [r*s, r^2*s^2*t, r^2*s^2*t^2, r*s*t^2, r*t^2, r^2*t, r^2*s^2, r*s*t]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 864
        @test code_k(c) == 16
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C3 x C3) : C3 [[864, 16, (≤ 24, ≤ 24)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤24, ≤24)")
    end

    @testset "C9 : C3 [[864, 28, (≤ 18, ≤ 18)]]" begin
        G = small_group(27, 4) # C9 : C3
        r, s, t = Oscar.gens(G)
        A = [r, r^2*t^2, r^2, r*t^2, r^2*t, r*t, r^2*s^2, r*s]
        B = [s^2, s, s*t^2, s^2*t, s^2*t^2, s*t, r^2*s, r*s^2*t]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 864
        @test code_k(c) == 28
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C9 : C3 [[864, 28, (≤ 18, ≤ 18)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤18, ≤18)")
    end

    @testset "C3 x C3 x C3 [[864, 24, (≤ 24, ≤ 24)]]" begin
        G = small_group(27, 5) # C3 x C3 x C3
        r, s, t = Oscar.gens(G)
        A = [s^2*t, s*t^2, r*t, r^2*t^2, s, s^2, r, r^2]
        B = [s*t, s^2*t^2, t^2, t, r*s, r^2*s^2, r^2*t, r*t^2]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 864
        @test code_k(c) == 24
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 x C3 x C3 [[864, 24, (≤ 24, ≤ 24)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤24, ≤24)")

    end

    @testset "C7 : C4 [[896, 40, (≤ 16, ≤ 16)]]" begin
        G = small_group(28, 1) # C7 : C4
        r, s, t = Oscar.gens(G)
        A = [r*s*t^3, r*t^3, r, r*s, r*s*t^5, r*t^5, r*t^4, r*s*t^4]
        B = [s*t, s*t^6, t, t^6, s*t^2, s*t^5, t^2, t^5]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 896
        @test code_k(c) == 40
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C7 : C4 [[896, 40, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C14 x C2 [[896, 16, (≤ 16, ≤ 16)]]" begin
        G = small_group(28, 4) # C14 x C2
        r, s, t = Oscar.gens(G)
        A = [s*t^2, s*t^5, t^6, t, r*s*t^4, r*s*t^3, r*t, r*t^6]
        B = [r, r*s, t^5, t^2, r*s*t^2, r*s*t^5, r*t^2, r*t^5]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 896
        @test code_k(c) == 16
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C14 x C2 [[896, 16, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C3 x D10 [[960, 24, (≤ 16, ≤ 16)]]" begin
        G = small_group(30, 2) # C3 x D10
        r, s, t = Oscar.gens(G)
        A = [r*s^2*t, r*s*t, s^2*t^3, s*t^2, r*s*t^4, r*s^2*t^4, r*s^2, r*s]
        B = [t^3, t^2, s*t, s^2*t^4, t^4, t, s, s^2]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 960
        @test code_k(c) == 24
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 x D10 [[960, 24, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    @testset "C3 x D10 [[960, 32, (≤ 16, ≤ 16)]]" begin
        G = small_group(30, 2) # C3 x D10
        r, s, t = Oscar.gens(G)
        A = [r, s^2*t, s*t^4, s, s^2, s^2*t^3, s*t^2, r*t^4]
        B = [t^2, t^3, r*s*t^4, r*s^2*t^4, t^4, t, r*s^2*t^3, r*s*t^3]
        c = QuantumTannerCode(G, A, B, ((H844, G844), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 960
        @test code_k(c) == 32
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 x D10 [[960, 32, (≤ 16, ≤ 16)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤16)")
    end

    # =====================
    # LRCC appendix table 2
    # =====================

    @testset "(C4 x C2) : C2 [[384, 24, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 13) # (C4 x C2) : C2
        r, s, t, u = Oscar.gens(G)
        A = [r, u, s*t*u, s*t, r*s, r*s*u]
        B = [s, r*t*u, r*t, r*s*t, r*s*t*u, t, t*u, s*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 384
        @test code_k(c) == 24
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C4 x C2) : C2 [[384, 24, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "(C4 x C2) : C2 [[384, 8, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 13) # (C4 x C2) : C2
        r, s, t, u = Oscar.gens(G)
        A = [s*t, s*t*u, r*t, r*t*u, u, s*u]
        B = [r*s, r*s*u, r, r*u, r*s*t*u, t*u, t, r*s*t]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 384
        @test code_k(c) == 8
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C4 x C2) : C2 [[384, 8, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C2 x D8 [[384, 44, (≤ 10, ≤ 10)]]" begin
        G = small_group(16, 11) # C2 x D8
        r, s, t, u = Oscar.gens(G)
        A = [r*s, r*s*u, s, r*s*t*u, r*s*t, s*u]
        B = [r*u, t, r*t*u, r, t*u, s*t*u, r*t, u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 384
        @test code_k(c) == 44
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x D8 [[384, 44, (≤ 10, ≤ 10)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤10, ≤10)")
    end

    @testset "C2 x D8 [[384, 38, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 11) # C2 x D8
        r, s, t, u = Oscar.gens(G)
        A = [s, r*s, r*s*u, r*u, s*u, u]
        B = [r*t, t*u, t, s*t*u, r*s*t*u, r*s*t, s*t, r*t*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 384
        @test code_k(c) == 38
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x D8 [[384, 38, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "D16 [[384, 32, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 7) # D16
        r, s, t, u = Oscar.gens(G)
        A = [s, s*t*u, t, t*u, s*t, s*u]
        B = [r*s*t*u, r*s*u, r*t, r*u, r*s, r*s*t, r, r*t*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 384
        @test code_k(c) == 32
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("D16 [[384, 32, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C4 x C2 x C2 [[384, 16, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 10) # C4 x C2 x C2
        r, s, t, u = Oscar.gens(G)
        A = [s, r*t*u, r*t, r*s*t, r*s*t*u, u]
        B = [r*s, r*s*u, s*t, r*u, r, s*t*u, t, t*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 384
        @test code_k(c) == 16
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C4 x C2 x C2 [[384, 16, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C8 : C2 [[384, 32, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 6) # C8 : C2
        r, s, t, u = Oscar.gens(G)
        A = [s*t*u, s*t, u, t*u, t, s*u]
        B = [r*s*u, r*s*t*u, r*s*t, r*s, r, r*t*u, r*u, r*t]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 384
        @test code_k(c) == 32
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C8 : C2 [[384, 32, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C8 : C2 [[384, 8, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 6) # C8 : C2
        r, s, t, u = Oscar.gens(G)
        A = [r*s*t, r*s, r*s*t*u, r*s*u, s*u, s]
        B = [r*t, r*u, t, t*u, s*t, s*t*u, r, r*t*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 384
        @test code_k(c) == 8
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C8 : C2 [[384, 8, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C8 x C2 [[384, 18, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 5) # C8 x C2
        r, s, t, u = Oscar.gens(G)
        A = [s*u, r*t*u, r, t*u, t, s]
        B = [r*s*t*u, r*s, r*t, r*u, s*t*u, s*t, r*s*u, r*s*t]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 384
        @test code_k(c) == 18
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C8 x C2 [[384, 18, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C4 : C4 [[384, 16, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 4) # C4 : C4
        r, s, t, u = Oscar.gens(G)
        A = [r*s*u, r*s, r*s*t, r*s*t*u, t*u, u]
        B = [s*u, s*t*u, s*t, s, r, r*u, r*t, r*t*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 384
        @test code_k(c) == 16
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C4 : C4 [[384, 16, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C4 : C4 [[384, 48, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 4) # C4 : C4
        r, s, t, u = Oscar.gens(G)
        A = [s*u, s*t*u, t, s, s*t, u]
        B = [r*s*t*u, r*s*t, r*t*u, r*t, r*s*u, r*s, r*u, r]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 384
        @test code_k(c) == 48
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C4 : C4 [[384, 48, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C4 x C4 [[392, 12, (≤ 13, ≤ 12)]]" begin
        G = small_group(16, 2) # C4 x C4
        r, s, t, u = Oscar.gens(G)
        A = [s*t*u, s*t, r*t, r, t*u, r*u, r*t*u]
        B = [s, s*u, t, r*s, r*s*t*u, r*s*u, r*s*t]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 392
        @test code_k(c) == 12
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C4 x C4 [[392, 12, (≤ 13, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤13, ≤12)")
    end

    @testset "C4 : C4 [[392, 15, (≤ 14, ≤ 10)]]" begin
        G = small_group(16, 4) # C4 : C4
        r, s, t, u = Oscar.gens(G)
        A = [t, s, s*t, r*u, r, r*t, r*t*u]
        B = [s*u, s*t*u, r*s*t*u, r*s*t, t*u, r*s, r*s*u]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 392
        @test code_k(c) == 15
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C4 : C4 [[392, 15, (≤ 14, ≤ 10)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤14, ≤10)")
    end

    @testset "C2 x D8 [[392, 36, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 11) # C2 x D8
        r, s, t, u = Oscar.gens(G)
        A = [r, s, r*u, r*s*u, r*s, u, s*u]
        B = [r*s*t*u, r*s*t, t, s*t, r*t*u, s*t*u, r*t]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 392
        @test code_k(c) == 36
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x D8 [[392, 36, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C2 x Q8 [[392, 8, (≤ 19, ≤ 20)]]" begin
        G = small_group(16, 12) # C2 x Q8
        r, s, t, u = Oscar.gens(G)
        A = [r*s*t*u, r*s*t, s, s*u, s*t, s*t*u, u]
        B = [r*s*u, r*s, r, r*u, t*u, r*t, r*t*u]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 392
        @test code_k(c) == 8
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x Q8 [[392, 8, (≤ 19, ≤ 20)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤19, ≤20)")
    end

    @testset "C2 x Q8 [[392, 15, (≤ 20, ≤ 10)]]" begin
        G = small_group(16, 12) # C2 x Q8
        r, s, t, u = Oscar.gens(G)
        A = [s*t*u, s*t, r*t, r*t*u, t, s*u, s]
        B = [r*s*t*u, r*s*t, r*s, r*s*u, u, r, r*u]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 392
        @test code_k(c) == 15
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x Q8 [[392, 15, (≤ 20, ≤ 10)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤20, ≤10)")
    end

    @testset "(C4 x C2) : C2 [[392, 33, (≤ 12, ≤ 12)]]" begin
        G = small_group(16, 13) # (C4 x C2) : C2
        r, s, t, u = Oscar.gens(G)
        A = [u, r*t*u, r*t, r*s, r*s*u, s*t*u, s*t]
        B = [t, t*u, s, r*s*t, r, r*s*t*u, r*u]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 392
        @test code_k(c) == 33
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C4 x C2) : C2 [[392, 33, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "(C4 x C2) : C2 [[392, 36, (≤ 12, ≤ 10)]]" begin
        G = small_group(16, 13) # (C4 x C2) : C2
        r, s, t, u = Oscar.gens(G)
        A = [s*t*u, s*t, u, t*u, t, r*s*u, r*s]
        B = [r*s*t*u, r*s*t, r, r*u, r*t, r*t*u, s*u]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 392
        @test code_k(c) == 36
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C4 x C2) : C2 [[392, 36, (≤ 12, ≤ 10)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤10)")
    end

    @testset "(C3 x C3) : C2 [[432, 4, (≤ 24, ≤ 24)]]" begin
        G = small_group(18, 4) # (C3 x C3) : C2
        r, s, t = Oscar.gens(G)
        A = [t, t^2, s^2, s, s*t^2, s^2*t]
        B = [r*s, r*s^2, r, r*t, r*s*t^2, r*s*t, r*t^2, r*s^2*t^2]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 432
        @test code_k(c) == 4
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C3 x C3) : C2 [[432, 4, (≤ 24, ≤ 24)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤24, ≤24)")
    end

    @testset "C6 x C3 [[432, 12, (≤ 20, ≤ 20)]]" begin
        G = small_group(18, 5) # C6 x C3
        r, s, t = Oscar.gens(G)
        A = [s^2*t, s*t^2, r*s^2*t, r*s*t^2, t^2, t]
        B = [r*s^2, r*s, r*s^2*t^2, r*s*t, r*t^2, r*t, s^2*t^2, s*t]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 432
        @test code_k(c) == 12
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C6 x C3 [[432, 12, (≤ 20, ≤ 20)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤20, ≤20)")
    end

    # =====================
    # LRCC appendix table 3
    # =====================

    @testset "C3 : C8 [[576, 16, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 1) # C3 : C8
        r, s, t, u = Oscar.gens(G)
        A = [s*u^2, s*t*u, s*t*u^2, s*u, t*u^2, t*u]
        B = [r*s*t, r, r*t*u^2, r*s*u^2, u, u^2, r*s*t*u^2, r*u^2]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 16
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 : C8 [[576, 16, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "SL(2,3) [[576, 12, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 3) # SL(2,3)
        r, s, t, u = Oscar.gens(G)
        A = [r*s, r^2*s*t*u, r*s*t, r^2*t*u, r*t, r^2*s*u]
        B = [r^2*s*t, r*s*u, r^2*s, r*t*u, t, t*u, r^2*t, r*s*t*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 12
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("SL(2,3) [[576, 12, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C3 : Q8 [[576, 8], (≤ 12, ≤ 12)]" begin
        G = small_group(24, 4) # C3 : Q8
        r, s, t, u = Oscar.gens(G)
        A = [s*t*u, s*u^2, r*s*t*u^2, r*s*u^2, r*s*u, r*s*t*u]
        B = [r*u, r*t*u, r*t*u^2, r*u^2, t*u, t*u^2, r*t, r]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 8
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 : Q8 [[576, 8], (≤ 12, ≤ 12)]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C4 x S3 [[576, 28, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 5) # C4 x S3
        r, s, t, u = Oscar.gens(G)
        A = [r*t*u^2, r*u^2, r*t, t*u^2, t*u, r*u]
        B = [r*s*u^2, r*s*t*u^2, s*u, s*t*u^2, r*s*t*u, r*s*u, r*s*t, r*s]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 28
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C4 x S3 [[576, 28, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "D24 [[576, 26, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 6) # D24
        r, s, t, u = Oscar.gens(G)
        A = [r*t*u^2, r*s*t, t*u^2, t*u, r*s*t*u^2, t]
        B = [u^2, u, s*u^2, s*t*u, s, s*t, s*u, s*t*u^2]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 26
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("D24 [[576, 26, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C2 x (C3 : C4) [[576, 16, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 7) # C2 x (C3 : C4)
        r, s, t, u = Oscar.gens(G)
        A = [t*u^2, t*u, s*t*u^2, s*t*u, t, s]
        B = [r*t, r, u^2, u, r*s*u^2, r*s*t*u^2, r*s*u, r*s*t*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 16
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x (C3 : C4) [[576, 16, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C2 x (C3 : C4) [[576, 40, (≤ 12, ≤ 9)]]" begin
        G = small_group(24, 7) # C2 x (C3 : C4)
        r, s, t, u = Oscar.gens(G)
        A = [s*t*u, s*t*u^2, t, u^2, u, s*t]
        B = [r*s*t*u^2, r*s*u^2, s*u, s*u^2, r*t*u^2, r*u^2, r*u, r*t*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 40
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x (C3 : C4) [[576, 40, (≤ 12, ≤ 9)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤9)")
    end

    @testset "(C6 x C2) : C2 [[576, 8, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 8) # (C6 x C2) : C2
        r, s, t, u = Oscar.gens(G)
        A = [t*u^2, t*u, u^2, u, s*t*u^2, s*t*u]
        B = [r*t, r*u^2, r*t*u, r*u, r*s*u, r*s*t*u, r*s*t*u^2, r*s*u^2]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 8
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C6 x C2) : C2 [[576, 8, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")

    end

    @testset "(C6 x C2) : C2 [[576, 28, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 8) # (C6 x C2) : C2
        r, s, t, u = Oscar.gens(G)
        A = [r*u, u^2, u, s*t, r*u^2, r*t]
        B = [t*u, t*u^2, r*s*t*u^2, r*s*u^2, r*s*u, r*s*t*u, r*s*t, r*s]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 28
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C6 x C2) : C2 [[576, 28, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C12 x C2 [[576, 8, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 9) # C12 x C2
        r, s, t, u = Oscar.gens(G)
        A = [s*u, s*t^2, s*t, t*u, t^2*u, s]
        B = [s*t*u, s*t^2*u, t^2, t, r, r*u, r*s*u, r*s]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 8
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C12 x C2 [[576, 8, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C12 x C2 [[576, 30, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 9) # C12 x C2
        r, s, t, u = Oscar.gens(G)
        A = [r*s*t, r*s*t^2*u, u, r*s, r*s*u, s*u]
        B = [r*t*u, r*t^2, s*t*u, s*t^2*u, r*s*t^2, r*s*t*u, r*u, r]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 30
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C12 x C2 [[576, 30, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C3 x D8 [[576, 16, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 10) # C3 x D8
        r, s, t, u = Oscar.gens(G)
        A = [r*t*u, r*t^2*u, s*u, r*s, r*s*u, s]
        B = [r*s*t, r*s*t^2*u, t, t^2, s*t, s*t^2, r*s*t^2, r*s*t*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 16
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 x D8 [[576, 16, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C3 x D8 [[576, 32, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 10) # C3 x D8
        r, s, t, u = Oscar.gens(G)
        A = [r, r*s*u, r*s, s*u, r*u, s]
        B = [r*t^2*u, r*t*u, r*s*t^2, r*s*t*u, s*t^2, s*t, s*t^2*u, s*t*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 32
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C3 x D8 [[576, 32, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "S4 [[576, 23, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 12) # S4
        r, s, t, u = Oscar.gens(G)
        A = [r*s^2*t*u, r*s^2*u, r, r*s, r*s*t, r*s*t*u]
        B = [s, s^2, s*u, s^2*t, u, t*u, s*t, s^2*t*u]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 23
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("S4 [[576, 23, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C2 x A4 [[576, 6, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 13) # C2 x A4
        r, s, t, u = Oscar.gens(G)
        A = [s^2*t, s*u, s^2*u, s*t*u, r*t*u, r*u]
        B = [r*s*u, r*s^2*t, r*s^2*u, r*s*t*u, r, r*s^2*t*u, r*s*t, t]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 6
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x A4 [[576, 6, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C2 x A4 [[576, 28, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 13) # C2 x A4
        r, s, t, u = Oscar.gens(G)
        A = [t*u, r*s, r*s^2, r*s*u, r*s^2*t, t]
        B = [s*u, s^2*t, s^2, s, s*t*u, s^2*u, s^2*t*u, s*t]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 28
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x A4 [[576, 28, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C2 x C2 x S3 [[576, 20, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 14) # C2 x C2 x S3
        r, s, t, u = Oscar.gens(G)
        A = [s*t*u^2, s*t*u, u^2, u, s*t, r*u^2]
        B = [r*s*t*u^2, r*s*t, s*u, s*u^2, r*s*u^2, r*t, r*s*t*u, r*s]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 20
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x C2 x S3 [[576, 20, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C2 x C2 x S3 [[576, 16, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 14) # C2 x C2 x S3
        r, s, t, u = Oscar.gens(G)
        A = [u^2, u, r*t*u^2, t, s*u^2, s*u]
        B = [r, r*s*t*u, r*s*t, r*s, t*u, t*u^2, s*t, r*s*u^2]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 16
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x C2 x S3 [[576, 16, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "C6 x C2 x C2 [[576, 26, (≤ 12, ≤ 12)]]" begin
        G = small_group(24, 15) # C6 x C2 x C2
        r, s, t, u = Oscar.gens(G)
        A = [r*s*u^2, r*s*u, u, u^2, r*s, r*t]
        B = [s*t, s*t*u^2, s*t*u, r*s*t*u, r*s*t*u^2, s*u, s*u^2, r*s*t]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 576
        @test code_k(c) == 26
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C6 x C2 x C2 [[576, 26, (≤ 12, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤12)")
    end

    @testset "(C3 x C3) : C3 [[648, 8, (≤ 20, ≤ 20)]]" begin
        G = small_group(27, 3) # (C3 x C3) : C3
        r, s, t = Oscar.gens(G)
        A = [r*s^2*t^2, r^2*s, r*t^2, r^2*t, t, t^2]
        B = [r*s*t^2, r^2*s^2*t^2, s^2*t, s*t^2, s, s^2, s^2*t^2, s*t]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 648
        @test code_k(c) == 8
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C3 x C3) : C3 [[648, 8, (≤ 20, ≤ 20)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤20, ≤20)")
    end

    # ======================
    # LRCC appendix table 4
    # ======================

    @testset "C7 : C4 [[672, 18, (≤ 24, ≤ 24)]]" begin
        G = small_group(28, 1) # C7 : C4
        r, s, t = Oscar.gens(G)
        A = [s*t^2, s*t^5, s*t^4, s*t^3, t, t^6]
        B = [r*t^5, r*s*t^5, r*s*t, r*t, s*t, s*t^6, r*t^6, r*s*t^6]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 672
        @test code_k(c) == 18
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C7 : C4 [[672, 18, (≤ 24, ≤ 24)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤24, ≤24)")
    end

    @testset "D28 [[672, 12, (≤ 22, ≤ 22)]]" begin
        G = small_group(28, 3) # D28
        r, s, t = Oscar.gens(G)
        A = [r*s*t, s*t^4, s*t^3, t^3, t^4, r*s*t^5]
        B = [t^2, t^5, s*t^2, s*t^5, t^6, t, s*t, s*t^6]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 672
        @test code_k(c) == 12
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("D28 [[672, 12, (≤ 22, ≤ 22)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤22, ≤22)")
    end

    @testset "C14 x C2 [[672, 20, (≤ 21, ≤ 21)]]" begin
        G = small_group(28, 4) # C14 x C2
        r, s, t = Oscar.gens(G)
        A = [s, s*t^3, s*t^4, t^3, t^4, r*s]
        B = [s*t^2, s*t^5, r*s*t^6, r*s*t, t^2, t^5, t^6, t]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 672
        @test code_k(c) == 20
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C14 x C2 [[672, 20, (≤ 21, ≤ 21)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤21, ≤21)")
    end

    @testset "C5 x S3 [[720, 32, (≤ 15, ≤ 15)]]" begin
        G = small_group(30, 1) # C5 x S3
        r, s, t = Oscar.gens(G)
        A = [s^2*t, s^3*t^2, r*t, r*s^2*t, r*s^3*t, r]
        B = [s*t^2, s^4*t, r*s^4, r*s, t, t^2, r*s*t, r*s^4*t]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 720
        @test code_k(c) == 32
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C5 x S3 [[720, 32, (≤ 15, ≤ 15)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤15, ≤15)")
    end

    @testset "D30 [[720, 16, (≤ 20, ≤ 20)]]" begin
        G = small_group(30, 3) # D30
        r, s, t = Oscar.gens(G)
        A = [s*t^2, s^2*t^3, r*s*t^2, s^2*t, s*t^4, r*s^2*t^3]
        B = [s^2, s, s^2*t^4, s*t, s^2*t^2, s*t^3, t^2, t^3]
        c = QuantumTannerCode(G, A, B, ((H633, G633), (H844, G844)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 720
        @test code_k(c) == 16
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 12
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("D30 [[720, 16, (≤ 20, ≤ 20)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤20, ≤20)")
    end

    @testset "C2 x A4 [[588, 33, (≤ 19, ≤ 18)]]" begin
        G = small_group(24, 13) # C2 x A4
        r, s, t, u = Oscar.gens(G)
        A = [r*s^2*t*u, r*s*t, r*t, r*s^2, r*s, r*s^2*u, r*s*t*u]
        B = [s^2*u, s*t*u, s, s^2, s*t, s^2*t*u, t]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 588
        @test code_k(c) == 33
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x A4 [[588, 33, (≤ 19, ≤ 18)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤19, ≤18)")
    end

    @testset "C2 x C2 x S3 [[588, 29, (≤ 18, ≤ 12)]]" begin
        G = small_group(24, 14) # C2 x C2 x S3
        r, s, t, u = Oscar.gens(G)
        A = [s*u, s*u^2, s*t*u, s*t*u^2, t, s*t, r]
        B = [t*u, t*u^2, u^2, u, s, r*s*u, r*t*u]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 588
        @test code_k(c) == 29
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x C2 x S3 [[588, 29, (≤ 18, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤18, ≤12)")
    end

    @testset "C2 x C2 x S3 [[588, 39, (≤ 12, ≤ 9)]]" begin
        G = small_group(24, 14) # C2 x C2 x S3
        r, s, t, u = Oscar.gens(G)
        A = [s*t, r*s*t*u, r*t*u, s*u^2, s*u, r*s*t*u^2, r*s*t]
        B = [t*u, t*u^2, u^2, u, r*s*u^2, r, s]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 588
        @test code_k(c) == 39
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C2 x C2 x S3 [[588, 39, (≤ 12, ≤ 9)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤9)")
    end

    @testset "C6 x C2 x C2 [[588, 40, (≤ 9, ≤ 18)]]" begin
        G = small_group(24, 15) # C6 x C2 x C2
        r, s, t, u = Oscar.gens(G)
        A = [r*s*u^2, r*s*u, u, u^2, r*s, r*t, s*t]
        B = [s*t*u^2, s*t*u, r*s*t*u, r*s*t*u^2, s*u, s*u^2, r*s*t]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 588
        @test code_k(c) == 40
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C6 x C2 x C2 [[588, 40, (≤ 9, ≤ 18)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤9, ≤18)")
    end

    @testset "C6 x C2 x C2 [[588, 38, (≤ 14, ≤ 10)]]" begin
        G = small_group(24, 15) # C6 x C2 x C2
        r, s, t, u = Oscar.gens(G)
        A = [s*u, s*u^2, r*t*u, r*t*u^2, r*s*t*u, r*s*t*u^2, r*t]
        B = [t*u^2, t*u, r*s*t, u, u^2, s*t*u^2, s*t*u]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 588
        @test code_k(c) == 38
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C6 x C2 x C2 [[588, 38, (≤ 14, ≤ 10)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤14, ≤10)")
    end

    @testset "C6 x C2 x C2 [[588, 21, (≤ 21, ≤ 12)]]" begin
        G = small_group(24, 15) # C6 x C2 x C2
        r, s, t, u = Oscar.gens(G)
        A = [r*t, s*u, s*u^2, r*s, r*s*u^2, r*s*u, s]
        B = [u^2, u, r*s*t*u, r*s*t*u^2, s*t*u^2, s*t*u, s*t]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 588
        @test code_k(c) == 21
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C6 x C2 x C2 [[588, 21, (≤ 21, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤21, ≤12)")
    end

    @testset "C8 x C4 [[784, 16, (≤ 19, ≤ 25)]]" begin
        G = small_group(32, 3) # C8 x C4
        r, s, t, u, v = Oscar.gens(G)
        A = [s*t*u*v, s*t, t*u*v, t*u, r*s, r*s*t*u*v, v]
        B = [r*t, r*v, r*s*v, r*s*t*u, r*t*u, r*u*v, u*v]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 784
        @test code_k(c) == 16
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("C8 x C4 [[784, 16, (≤ 19, ≤ 25)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤19, ≤25)")
    end

    @testset "(C8 x C2) : C2 [[784, 37, (≤ 14, ≤ 12)]]" begin
        G = small_group(32, 5) # (C8 x C2) : C2
        r, s, t, u, v = Oscar.gens(G)
        A = [r*t*u, r*t*v, t, t*u, t*u*v, r*s*t*u, r*s*v]
        B = [s*t*v, v, r*s*t*u*v, r*s, r*t, r*t*u*v, s*v]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 784
        @test code_k(c) == 37
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C8 x C2) : C2 [[784, 37, (≤ 14, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤14, ≤12)")
    end

    @testset "(C2 x C2 x C2) : C4 [[784, 40, (≤ 14, ≤ 12)]]" begin
        G = small_group(32, 6) # (C2 x C2 x C2) : C4
        r, s, t, u, v = Oscar.gens(G)
        A = [r*s*t, r*s*u*v, v, r*s, r*s*t*u, r*s*t*v, r*s*u]
        B = [r*t, r*t*u*v, r*u*v, r*v, u*v, t, s*t*v]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 784
        @test code_k(c) == 40
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C2 x C2 x C2) : C4 [[784, 40, (≤ 14, ≤ 12)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤14, ≤12)")
    end

    @testset "(C8 : C2) : C2 [[784, 24, (≤ 16, ≤ 13)]]" begin
        G = small_group(32, 7) # (C8 : C2) : C2
        r, s, t, u, v = Oscar.gens(G)
        A = [r*s*t, r*s*u, v, r*s, r*s*t*u*v, r*s*t*v, r*s*u*v]
        B = [r*t, r*t*u, r*u*v, r, u*v, u, t]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 784
        @test code_k(c) == 24
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C8 : C2) : C2 [[784, 24, (≤ 16, ≤ 13)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤16, ≤13)")
    end

    @testset "(C8 : C2) : C2 [[784, 58, (≤ 12, ≤ 10)]]" begin
        G = small_group(32, 7) # (C8 : C2) : C2
        r, s, t, u, v = Oscar.gens(G)
        A = [r*t*u, r*t, r*v, r*u, s, r*u*v, r]
        B = [s*u*v, t*v, r*s, r*s*t*u*v, s*t*u, t, v]
        c = QuantumTannerCode(G, A, B, ((H743, G743), (H734, G734)))
        hx, hz = parity_matrix_xz(c)
        stab = QuantumClifford.ECC.parity_checks(c)
        mat = matrix(GF(2), stab_to_gf2(stab))
        computed_rank = rank(mat)
        @test computed_rank == code_n(c) - code_k(c)
        @test code_n(c) == 784
        @test code_k(c) == 58
        wx = maximum(vec(sum(hx, dims=2)))
        wz = maximum(vec(sum(hz, dims=2)))
        @test max(wx, wz) == 16
        dx_qdist, dz_qdist = compute_qdistrnd_distance(
            hx,
            hz;
            num=QDIST_TRIALS,
        )
        println("(C8 : C2) : C2 [[784, 58, (≤ 12, ≤ 10)]]")
        println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
        println("  sQetch paper : (dx, dz) = (≤12, ≤10)")
    end
end
