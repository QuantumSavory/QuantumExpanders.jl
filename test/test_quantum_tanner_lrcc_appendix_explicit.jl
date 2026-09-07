@testitem "Quantum Tanner Codes. All LRCC appendix constructions" begin
    using Test
    using Oscar
    using QECCore: parity_matrix_x, parity_matrix_z, code_n, code_k
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

    const QDIST_TRIALS = 1000

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

    function test_lrcc_appendix(
            name,
            G,
            A,
            B,
            local_codes;
            n,
            k,
            max_weight,
            dx_bound,
            dz_bound,
            distance_slack=2,
        )
        title = "$name [[$n, $k, (≤ $dx_bound, ≤ $dz_bound)]]"

        @testset "$title" begin
            c = QuantumTannerCode(G, A, B, local_codes)
            hx, hz = parity_matrix_x(c), parity_matrix_z(c)
            stab = QuantumClifford.ECC.parity_checks(c)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == n
            @test code_k(c) == k
            wx = maximum(vec(sum(hx, dims=2)))
            wz = maximum(vec(sum(hz, dims=2)))
            @test max(wx, wz) == max_weight 
        end
    end

    G = small_group(18, 3) # C3 x S3
    r, s, t = Oscar.gens(G)
    A = [s, s^2, t, t^2, r*t^2, r]
    B = [r*s, r*s^2, s*t, s^2*t^2, r*s*t^2, r*s^2*t^2]

    test_lrcc_appendix(
        "C3 x S3",
        G,
        A,
        B,
        ((H633, G633), (H633, G633));
        n=324,
        k=8,
        max_weight=9,
        dx_bound=17,
        dz_bound=14,
    )

    G = small_group(18, 3) # C3 x S3
    r, s, t = Oscar.gens(G)
    A = [r*s*t^2, r*s^2*t^2, r*s^2*t, r*s*t, t^2, t, r*s, r*s^2]
    B = [r*t^2, s^2, s, s*t, s^2*t^2, s^2*t, s*t^2, r]

    test_lrcc_appendix(
        "C3 x S3",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=576,
        k=22,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(24, 13) # C2 x A4
    r, s, t, u = Oscar.gens(G)
    A = [s^2*t*u, s*t, t, s^2*u, s*t*u, r*t]
    B = [r*s^2*t*u, r*s*t, r*s^2*u, r*s*t*u, r*s^2, r*s]

    test_lrcc_appendix(
        "C2 x A4",
        G,
        A,
        B,
        ((H633, G633), (H633, G633));
        n=432,
        k=16,
        max_weight=9,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 12) # S4
    r, s, t, u = Oscar.gens(G)
    A = [s^2*t, s*u, s^2, s, s*t*u, s^2*u, s*t, s^2*t*u]
    B = [r*s^2*u, r*s^2*t*u, r*t, r*u, r*s^2, r*s*t, r*s*t*u, r]

    test_lrcc_appendix(
        "S4",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=768,
        k=32,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(24, 13) # C2 x A4
    r, s, t, u = Oscar.gens(G)
    A = [r, r*s^2*t*u, r*s*t, r*s^2*t, r*s*u, r*s^2*u, r*s*t*u, t]
    B = [s, s^2, s*t, s^2*t*u, s*t*u, s^2*u, r*u, r*t*u]

    test_lrcc_appendix(
        "C2 x A4",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=768,
        k=42,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(24, 13) # C2 x A4
    r, s, t, u = Oscar.gens(G)
    A = [s, s^2, s^2*t*u, s*t, r*t, s*t*u, s^2*u, r]
    B = [t, r*s*t*u, r*s^2*u, r*s, r*s^2, r*s*t, r*s^2*t*u, t*u]

    test_lrcc_appendix(
        "C2 x A4",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=768,
        k=20,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(24, 9) # C12 x C2
    r, s, t, u = Oscar.gens(G)
    A = [s*u, r*t^2, r*t*u, r*s*t^2, r*s*t*u, t^2*u, t*u, u]
    B = [s*t^2*u, s*t*u, r*s*u, r*s, t^2, t, r*u, r]

    test_lrcc_appendix(
        "C12 x C2",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=768,
        k=40,
        max_weight=16,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 8) # (C6 x C2) : C2
    r, s, t, u = Oscar.gens(G)
    A = [r*u, u^2, u, s*t, r*u^2, t*u, t*u^2, r*t]
    B = [r*s*t*u^2, r*s*u^2, r*s*u, r*s*t*u, r*s*t, r*s, s*t*u^2, s*t*u]

    test_lrcc_appendix(
        "(C6 x C2) : C2",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=768,
        k=38,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(24, 7) # C2 x (C3 : C4)
    r, s, t, u = Oscar.gens(G)
    A = [r*u^2, r*t*u^2, t*u, t*u^2, s*t*u^2, s*t*u, s*u, s*u^2]
    B = [u^2, u, r*s*u^2, r*s*t*u^2, t, s, r*s*t, r*s]

    test_lrcc_appendix(
        "C2 x (C3 : C4)",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=768,
        k=34,
        max_weight=16,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 5) # C4 x S3
    r, s, t, u = Oscar.gens(G)
    A = [r*s*t*u^2, r*s*u^2, r*s*t, r*s, s, s*t, r*s*t*u, r*s*u]
    B = [r*t, t*u^2, t*u, s*u^2, s*t*u, r*u^2, t, r*t*u^2]

    test_lrcc_appendix(
        "C4 x S3",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=768,
        k=48,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(24, 4) # C3 : Q8
    r, s, t, u = Oscar.gens(G)
    A = [r*s*t*u^2, r*s*u^2, r*s*t, r*s, s, s*t, r*s*t*u, r*s*u]
    B = [r*t, r, t*u^2, t*u, s*u^2, s*t*u, r*u^2, r*t*u^2]

    test_lrcc_appendix(
        "C3 : Q8",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=768,
        k=28,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(24, 3) # SL(2,3)
    r, s, t, u = Oscar.gens(G)
    A = [r^2*s*t*u, r*s, r^2, r, r^2*s*u, r*t, s, s*u]
    B = [r*s*u, r^2*s*t, r*s*t*u, r^2*t, r*u, r^2*u, r*t*u, r^2*s]

    test_lrcc_appendix(
        "SL(2,3)",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=768,
        k=32,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(24, 1) # C3 : C8
    r, s, t, u = Oscar.gens(G)
    A = [r*u, r*s*t*u, r*u^2, r*s*t*u^2, s, s*t, t*u^2, t*u]
    B = [r*s, r*t, s*u^2, s*t*u, r*s*u^2, r*t*u^2, u, u^2]

    test_lrcc_appendix(
        "C3 : C8",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=768,
        k=32,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(24, 1) # C3 : C8
    r, s, t, u = Oscar.gens(G)
    A = [r*s*t, r, s*t*u, s*u^2, r*s*t*u, r*u, s*t*u^2, s*u]
    B = [r*s, r*t, t*u^2, t*u, s, s*t, r*t*u, r*s*u]

    test_lrcc_appendix(
        "C3 : C8",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=768,
        k=24,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(25, 2) # C5 x C5
    r, s = Oscar.gens(G)
    A = [r*s^2, r^4*s^3, r^3*s, r^2*s^4, r^4*s^2, r*s^3, r, r^4]
    B = [r^3, r^2, r^4*s, r*s^4, s, s^4, r^2*s, r^3*s^4]

    test_lrcc_appendix(
        "C5 x C5",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=800,
        k=24,
        max_weight=16,
        dx_bound=20,
        dz_bound=20,
    )

    G = small_group(27, 2) # C9 x C3
    r, s, t = Oscar.gens(G)
    A = [r^2*t^2, r, t, t^2, r^2, r*t^2, r*s^2, r^2*s*t^2]
    B = [r^2*s^2, r*s*t^2, r^2*s*t, r*s^2*t, s^2*t, s*t^2, r^2*s, r*s^2*t^2]

    test_lrcc_appendix(
        "C9 x C3",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=864,
        k=32,
        max_weight=16,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(27, 3) # (C3 x C3) : C3
    r, s, t = Oscar.gens(G)
    A = [r^2*t, r*t^2, t^2, t, r*t, r^2*t^2, r, r^2]
    B = [r^2*s*t, r*s^2*t, r*s*t^2, r^2*s^2*t^2, r^2*s^2*t, r*s, s*t^2, s^2*t]

    test_lrcc_appendix(
        "(C3 x C3) : C3",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=864,
        k=20,
        max_weight=16,
        dx_bound=24,
        dz_bound=24,
    )

    G = small_group(27, 3) # (C3 x C3) : C3
    r, s, t = Oscar.gens(G)
    A = [r^2*s, r*s^2*t^2, r^2*s*t, r*s^2*t, s^2, s, s*t, s^2*t^2]
    B = [r*s, r^2*s^2*t, r^2*s^2*t^2, r*s*t^2, r*t^2, r^2*t, r^2*s^2, r*s*t]

    test_lrcc_appendix(
        "(C3 x C3) : C3",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=864,
        k=16,
        max_weight=16,
        dx_bound=24,
        dz_bound=24,
    )

    G = small_group(27, 4) # C9 : C3
    r, s, t = Oscar.gens(G)
    A = [r, r^2*t^2, r^2, r*t^2, r^2*t, r*t, r^2*s^2, r*s]
    B = [s^2, s, s*t^2, s^2*t, s^2*t^2, s*t, r^2*s, r*s^2*t]

    test_lrcc_appendix(
        "C9 : C3",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=864,
        k=28,
        max_weight=16,
        dx_bound=18,
        dz_bound=18,
    )

    G = small_group(27, 5) # C3 x C3 x C3
    r, s, t = Oscar.gens(G)
    A = [s^2*t, s*t^2, r*t, r^2*t^2, s, s^2, r, r^2]
    B = [s*t, s^2*t^2, t^2, t, r*s, r^2*s^2, r^2*t, r*t^2]

    test_lrcc_appendix(
        "C3 x C3 x C3",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=864,
        k=24,
        max_weight=16,
        dx_bound=24,
        dz_bound=24,
    )

    G = small_group(28, 1) # C7 : C4
    r, s, t = Oscar.gens(G)
    A = [r*s*t^3, r*t^3, r, r*s, r*s*t^5, r*t^5, r*t^4, r*s*t^4]
    B = [s*t, s*t^6, t, t^6, s*t^2, s*t^5, t^2, t^5]

    test_lrcc_appendix(
        "C7 : C4",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=896,
        k=40,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(28, 4) # C14 x C2
    r, s, t = Oscar.gens(G)
    A = [s*t^2, s*t^5, t^6, t, r*s*t^4, r*s*t^3, r*t, r*t^6]
    B = [r, r*s, t^5, t^2, r*s*t^2, r*s*t^5, r*t^2, r*t^5]

    test_lrcc_appendix(
        "C14 x C2",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=896,
        k=16,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(30, 2) # C3 x D10
    r, s, t = Oscar.gens(G)
    A = [r*s^2*t, r*s*t, s^2*t^3, s*t^2, r*s*t^4, r*s^2*t^4, r*s^2, r*s]
    B = [t^3, t^2, s*t, s^2*t^4, t^4, t, s, s^2]

    test_lrcc_appendix(
        "C3 x D10",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=960,
        k=24,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    G = small_group(30, 2) # C3 x D10
    r, s, t = Oscar.gens(G)
    A = [r, s^2*t, s*t^4, s, s^2, s^2*t^3, s*t^2, r*t^4]
    B = [t^2, t^3, r*s*t^4, r*s^2*t^4, t^4, t, r*s^2*t^3, r*s*t^3]

    test_lrcc_appendix(
        "C3 x D10",
        G,
        A,
        B,
        ((H844, G844), (H844, G844));
        n=960,
        k=32,
        max_weight=16,
        dx_bound=16,
        dz_bound=16,
    )

    # =====================
    # LRCC appendix table 2
    # =====================

    G = small_group(16, 13) # (C4 x C2) : C2
    r, s, t, u = Oscar.gens(G)
    A = [r, u, s*t*u, s*t, r*s, r*s*u]
    B = [s, r*t*u, r*t, r*s*t, r*s*t*u, t, t*u, s*u]

    test_lrcc_appendix(
        "(C4 x C2) : C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=384,
        k=24,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 13) # (C4 x C2) : C2
    r, s, t, u = Oscar.gens(G)
    A = [s*t, s*t*u, r*t, r*t*u, u, s*u]
    B = [r*s, r*s*u, r, r*u, r*s*t*u, t*u, t, r*s*t]

    test_lrcc_appendix(
        "(C4 x C2) : C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=384,
        k=8,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 11) # C2 x D8
    r, s, t, u = Oscar.gens(G)
    A = [r*s, r*s*u, s, r*s*t*u, r*s*t, s*u]
    B = [r*u, t, r*t*u, r, t*u, s*t*u, r*t, u]

    test_lrcc_appendix(
        "C2 x D8",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=384,
        k=44,
        max_weight=12,
        dx_bound=10,
        dz_bound=10,
    )

    G = small_group(16, 11) # C2 x D8
    r, s, t, u = Oscar.gens(G)
    A = [s, r*s, r*s*u, r*u, s*u, u]
    B = [r*t, t*u, t, s*t*u, r*s*t*u, r*s*t, s*t, r*t*u]

    test_lrcc_appendix(
        "C2 x D8",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=384,
        k=38,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 7) # D16
    r, s, t, u = Oscar.gens(G)
    A = [s, s*t*u, t, t*u, s*t, s*u]
    B = [r*s*t*u, r*s*u, r*t, r*u, r*s, r*s*t, r, r*t*u]

    test_lrcc_appendix(
        "D16",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=384,
        k=32,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 10) # C4 x C2 x C2
    r, s, t, u = Oscar.gens(G)
    A = [s, r*t*u, r*t, r*s*t, r*s*t*u, u]
    B = [r*s, r*s*u, s*t, r*u, r, s*t*u, t, t*u]

    test_lrcc_appendix(
        "C4 x C2 x C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=384,
        k=16,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 6) # C8 : C2
    r, s, t, u = Oscar.gens(G)
    A = [s*t*u, s*t, u, t*u, t, s*u]
    B = [r*s*u, r*s*t*u, r*s*t, r*s, r, r*t*u, r*u, r*t]

    test_lrcc_appendix(
        "C8 : C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=384,
        k=32,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 6) # C8 : C2
    r, s, t, u = Oscar.gens(G)
    A = [r*s*t, r*s, r*s*t*u, r*s*u, s*u, s]
    B = [r*t, r*u, t, t*u, s*t, s*t*u, r, r*t*u]

    test_lrcc_appendix(
        "C8 : C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=384,
        k=8,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 5) # C8 x C2
    r, s, t, u = Oscar.gens(G)
    A = [s*u, r*t*u, r, t*u, t, s]
    B = [r*s*t*u, r*s, r*t, r*u, s*t*u, s*t, r*s*u, r*s*t]

    test_lrcc_appendix(
        "C8 x C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=384,
        k=18,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 4) # C4 : C4
    r, s, t, u = Oscar.gens(G)
    A = [r*s*u, r*s, r*s*t, r*s*t*u, t*u, u]
    B = [s*u, s*t*u, s*t, s, r, r*u, r*t, r*t*u]

    test_lrcc_appendix(
        "C4 : C4",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=384,
        k=16,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 4) # C4 : C4
    r, s, t, u = Oscar.gens(G)
    A = [s*u, s*t*u, t, s, s*t, u]
    B = [r*s*t*u, r*s*t, r*t*u, r*t, r*s*u, r*s, r*u, r]

    test_lrcc_appendix(
        "C4 : C4",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=384,
        k=48,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 2) # C4 x C4
    r, s, t, u = Oscar.gens(G)
    A = [s*t*u, s*t, r*t, r, t*u, r*u, r*t*u]
    B = [s, s*u, t, r*s, r*s*t*u, r*s*u, r*s*t]

    test_lrcc_appendix(
        "C4 x C4",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=392,
        k=12,
        max_weight=16,
        dx_bound=13,
        dz_bound=12,
    )

    G = small_group(16, 4) # C4 : C4
    r, s, t, u = Oscar.gens(G)
    A = [t, s, s*t, r*u, r, r*t, r*t*u]
    B = [s*u, s*t*u, r*s*t*u, r*s*t, t*u, r*s, r*s*u]

    test_lrcc_appendix(
        "C4 : C4",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=392,
        k=15,
        max_weight=16,
        dx_bound=14,
        dz_bound=10,
    )

    G = small_group(16, 11) # C2 x D8
    r, s, t, u = Oscar.gens(G)
    A = [r, s, r*u, r*s*u, r*s, u, s*u]
    B = [r*s*t*u, r*s*t, t, s*t, r*t*u, s*t*u, r*t]

    test_lrcc_appendix(
        "C2 x D8",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=392,
        k=36,
        max_weight=16,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 12) # C2 x Q8
    r, s, t, u = Oscar.gens(G)
    A = [r*s*t*u, r*s*t, s, s*u, s*t, s*t*u, u]
    B = [r*s*u, r*s, r, r*u, t*u, r*t, r*t*u]

    test_lrcc_appendix(
        "C2 x Q8",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=392,
        k=8,
        max_weight=16,
        dx_bound=19,
        dz_bound=20,
    )

    G = small_group(16, 12) # C2 x Q8
    r, s, t, u = Oscar.gens(G)
    A = [s*t*u, s*t, r*t, r*t*u, t, s*u, s]
    B = [r*s*t*u, r*s*t, r*s, r*s*u, u, r, r*u]

    test_lrcc_appendix(
        "C2 x Q8",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=392,
        k=15,
        max_weight=16,
        dx_bound=20,
        dz_bound=10,
    )

    G = small_group(16, 13) # (C4 x C2) : C2
    r, s, t, u = Oscar.gens(G)
    A = [u, r*t*u, r*t, r*s, r*s*u, s*t*u, s*t]
    B = [t, t*u, s, r*s*t, r, r*s*t*u, r*u]

    test_lrcc_appendix(
        "(C4 x C2) : C2",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=392,
        k=33,
        max_weight=16,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(16, 13) # (C4 x C2) : C2
    r, s, t, u = Oscar.gens(G)
    A = [s*t*u, s*t, u, t*u, t, r*s*u, r*s]
    B = [r*s*t*u, r*s*t, r, r*u, r*t, r*t*u, s*u]

    test_lrcc_appendix(
        "(C4 x C2) : C2",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=392,
        k=36,
        max_weight=16,
        dx_bound=12,
        dz_bound=10,
    )

    G = small_group(18, 4) # (C3 x C3) : C2
    r, s, t = Oscar.gens(G)
    A = [t, t^2, s^2, s, s*t^2, s^2*t]
    B = [r*s, r*s^2, r, r*t, r*s*t^2, r*s*t, r*t^2, r*s^2*t^2]

    test_lrcc_appendix(
        "(C3 x C3) : C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=432,
        k=4,
        max_weight=12,
        dx_bound=24,
        dz_bound=24,
    )

    G = small_group(18, 5) # C6 x C3
    r, s, t = Oscar.gens(G)
    A = [s^2*t, s*t^2, r*s^2*t, r*s*t^2, t^2, t]
    B = [r*s^2, r*s, r*s^2*t^2, r*s*t, r*t^2, r*t, s^2*t^2, s*t]

    test_lrcc_appendix(
        "C6 x C3",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=432,
        k=12,
        max_weight=12,
        dx_bound=20,
        dz_bound=20,
    )

    # =====================
    # LRCC appendix table 3
    # =====================

    G = small_group(24, 1) # C3 : C8
    r, s, t, u = Oscar.gens(G)
    A = [s*u^2, s*t*u, s*t*u^2, s*u, t*u^2, t*u]
    B = [r*s*t, r, r*t*u^2, r*s*u^2, u, u^2, r*s*t*u^2, r*u^2]

    test_lrcc_appendix(
        "C3 : C8",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=16,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 3) # SL(2,3)
    r, s, t, u = Oscar.gens(G)
    A = [r*s, r^2*s*t*u, r*s*t, r^2*t*u, r*t, r^2*s*u]
    B = [r^2*s*t, r*s*u, r^2*s, r*t*u, t, t*u, r^2*t, r*s*t*u]

    test_lrcc_appendix(
        "SL(2,3)",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=12,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 4) # C3 : Q8
    r, s, t, u = Oscar.gens(G)
    A = [s*t*u, s*u^2, r*s*t*u^2, r*s*u^2, r*s*u, r*s*t*u]
    B = [r*u, r*t*u, r*t*u^2, r*u^2, t*u, t*u^2, r*t, r]

    test_lrcc_appendix(
        "C3 : Q8",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=8,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 5) # C4 x S3
    r, s, t, u = Oscar.gens(G)
    A = [r*t*u^2, r*u^2, r*t, t*u^2, t*u, r*u]
    B = [r*s*u^2, r*s*t*u^2, s*u, s*t*u^2, r*s*t*u, r*s*u, r*s*t, r*s]

    test_lrcc_appendix(
        "C4 x S3",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=28,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 6) # D24
    r, s, t, u = Oscar.gens(G)
    A = [r*t*u^2, r*s*t, t*u^2, t*u, r*s*t*u^2, t]
    B = [u^2, u, s*u^2, s*t*u, s, s*t, s*u, s*t*u^2]

    test_lrcc_appendix(
        "D24",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=26,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 7) # C2 x (C3 : C4)
    r, s, t, u = Oscar.gens(G)
    A = [t*u^2, t*u, s*t*u^2, s*t*u, t, s]
    B = [r*t, r, u^2, u, r*s*u^2, r*s*t*u^2, r*s*u, r*s*t*u]

    test_lrcc_appendix(
        "C2 x (C3 : C4)",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=16,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 7) # C2 x (C3 : C4)
    r, s, t, u = Oscar.gens(G)
    A = [s*t*u, s*t*u^2, t, u^2, u, s*t]
    B = [r*s*t*u^2, r*s*u^2, s*u, s*u^2, r*t*u^2, r*u^2, r*u, r*t*u]

    test_lrcc_appendix(
        "C2 x (C3 : C4)",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=40,
        max_weight=12,
        dx_bound=12,
        dz_bound=9,
    )

    G = small_group(24, 8) # (C6 x C2) : C2
    r, s, t, u = Oscar.gens(G)
    A = [t*u^2, t*u, u^2, u, s*t*u^2, s*t*u]
    B = [r*t, r*u^2, r*t*u, r*u, r*s*u, r*s*t*u, r*s*t*u^2, r*s*u^2]

    test_lrcc_appendix(
        "(C6 x C2) : C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=8,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 8) # (C6 x C2) : C2
    r, s, t, u = Oscar.gens(G)
    A = [r*u, u^2, u, s*t, r*u^2, r*t]
    B = [t*u, t*u^2, r*s*t*u^2, r*s*u^2, r*s*u, r*s*t*u, r*s*t, r*s]

    test_lrcc_appendix(
        "(C6 x C2) : C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=28,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 9) # C12 x C2
    r, s, t, u = Oscar.gens(G)
    A = [s*u, s*t^2, s*t, t*u, t^2*u, s]
    B = [s*t*u, s*t^2*u, t^2, t, r, r*u, r*s*u, r*s]

    test_lrcc_appendix(
        "C12 x C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=8,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 9) # C12 x C2
    r, s, t, u = Oscar.gens(G)
    A = [r*s*t, r*s*t^2*u, u, r*s, r*s*u, s*u]
    B = [r*t*u, r*t^2, s*t*u, s*t^2*u, r*s*t^2, r*s*t*u, r*u, r]

    test_lrcc_appendix(
        "C12 x C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=30,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 10) # C3 x D8
    r, s, t, u = Oscar.gens(G)
    A = [r*t*u, r*t^2*u, s*u, r*s, r*s*u, s]
    B = [r*s*t, r*s*t^2*u, t, t^2, s*t, s*t^2, r*s*t^2, r*s*t*u]

    test_lrcc_appendix(
        "C3 x D8",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=16,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 10) # C3 x D8
    r, s, t, u = Oscar.gens(G)
    A = [r, r*s*u, r*s, s*u, r*u, s]
    B = [r*t^2*u, r*t*u, r*s*t^2, r*s*t*u, s*t^2, s*t, s*t^2*u, s*t*u]

    test_lrcc_appendix(
        "C3 x D8",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=32,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 12) # S4
    r, s, t, u = Oscar.gens(G)
    A = [r*s^2*t*u, r*s^2*u, r, r*s, r*s*t, r*s*t*u]
    B = [s, s^2, s*u, s^2*t, u, t*u, s*t, s^2*t*u]

    test_lrcc_appendix(
        "S4",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=23,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 13) # C2 x A4
    r, s, t, u = Oscar.gens(G)
    A = [s^2*t, s*u, s^2*u, s*t*u, r*t*u, r*u]
    B = [r*s*u, r*s^2*t, r*s^2*u, r*s*t*u, r, r*s^2*t*u, r*s*t, t]

    test_lrcc_appendix(
        "C2 x A4",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=6,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 13) # C2 x A4
    r, s, t, u = Oscar.gens(G)
    A = [t*u, r*s, r*s^2, r*s*u, r*s^2*t, t]
    B = [s*u, s^2*t, s^2, s, s*t*u, s^2*u, s^2*t*u, s*t]

    test_lrcc_appendix(
        "C2 x A4",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=28,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 14) # C2 x C2 x S3
    r, s, t, u = Oscar.gens(G)
    A = [s*t*u^2, s*t*u, u^2, u, s*t, r*u^2]
    B = [r*s*t*u^2, r*s*t, s*u, s*u^2, r*s*u^2, r*t, r*s*t*u, r*s]

    test_lrcc_appendix(
        "C2 x C2 x S3",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=20,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 14) # C2 x C2 x S3
    r, s, t, u = Oscar.gens(G)
    A = [u^2, u, r*t*u^2, t, s*u^2, s*u]
    B = [r, r*s*t*u, r*s*t, r*s, t*u, t*u^2, s*t, r*s*u^2]

    test_lrcc_appendix(
        "C2 x C2 x S3",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=16,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(24, 15) # C6 x C2 x C2
    r, s, t, u = Oscar.gens(G)
    A = [r*s*u^2, r*s*u, u, u^2, r*s, r*t]
    B = [s*t, s*t*u^2, s*t*u, r*s*t*u, r*s*t*u^2, s*u, s*u^2, r*s*t]

    test_lrcc_appendix(
        "C6 x C2 x C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=576,
        k=26,
        max_weight=12,
        dx_bound=12,
        dz_bound=12,
    )

    G = small_group(27, 3) # (C3 x C3) : C3
    r, s, t = Oscar.gens(G)
    A = [r*s^2*t^2, r^2*s, r*t^2, r^2*t, t, t^2]
    B = [r*s*t^2, r^2*s^2*t^2, s^2*t, s*t^2, s, s^2, s^2*t^2, s*t]

    test_lrcc_appendix(
        "(C3 x C3) : C3",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=648,
        k=8,
        max_weight=12,
        dx_bound=20,
        dz_bound=20,
    )

    # ======================
    # LRCC appendix table 4
    # ======================

    G = small_group(28, 1) # C7 : C4
    r, s, t = Oscar.gens(G)
    A = [s*t^2, s*t^5, s*t^4, s*t^3, t, t^6]
    B = [r*t^5, r*s*t^5, r*s*t, r*t, s*t, s*t^6, r*t^6, r*s*t^6]

    test_lrcc_appendix(
        "C7 : C4",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=672,
        k=18,
        max_weight=12,
        dx_bound=24,
        dz_bound=24,
    )

    G = small_group(28, 3) # D28
    r, s, t = Oscar.gens(G)
    A = [r*s*t, s*t^4, s*t^3, t^3, t^4, r*s*t^5]
    B = [t^2, t^5, s*t^2, s*t^5, t^6, t, s*t, s*t^6]

    test_lrcc_appendix(
        "D28",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=672,
        k=12,
        max_weight=12,
        dx_bound=22,
        dz_bound=22,
    )

    G = small_group(28, 4) # C14 x C2
    r, s, t = Oscar.gens(G)
    A = [s, s*t^3, s*t^4, t^3, t^4, r*s]
    B = [s*t^2, s*t^5, r*s*t^6, r*s*t, t^2, t^5, t^6, t]

    test_lrcc_appendix(
        "C14 x C2",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=672,
        k=20,
        max_weight=12,
        dx_bound=21,
        dz_bound=21,
    )

    G = small_group(30, 1) # C5 x S3
    r, s, t = Oscar.gens(G)
    A = [s^2*t, s^3*t^2, r*t, r*s^2*t, r*s^3*t, r]
    B = [s*t^2, s^4*t, r*s^4, r*s, t, t^2, r*s*t, r*s^4*t]

    test_lrcc_appendix(
        "C5 x S3",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=720,
        k=32,
        max_weight=12,
        dx_bound=15,
        dz_bound=15,
    )

    G = small_group(30, 3) # D30
    r, s, t = Oscar.gens(G)
    A = [s*t^2, s^2*t^3, r*s*t^2, s^2*t, s*t^4, r*s^2*t^3]
    B = [s^2, s, s^2*t^4, s*t, s^2*t^2, s*t^3, t^2, t^3]

    test_lrcc_appendix(
        "D30",
        G,
        A,
        B,
        ((H633, G633), (H844, G844));
        n=720,
        k=16,
        max_weight=12,
        dx_bound=20,
        dz_bound=20,
    )

    G = small_group(24, 13) # C2 x A4
    r, s, t, u = Oscar.gens(G)
    A = [r*s^2*t*u, r*s*t, r*t, r*s^2, r*s, r*s^2*u, r*s*t*u]
    B = [s^2*u, s*t*u, s, s^2, s*t, s^2*t*u, t]

    test_lrcc_appendix(
        "C2 x A4",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=588,
        k=33,
        max_weight=16,
        dx_bound=19,
        dz_bound=18,
    )

    G = small_group(24, 14) # C2 x C2 x S3
    r, s, t, u = Oscar.gens(G)
    A = [s*u, s*u^2, s*t*u, s*t*u^2, t, s*t, r]
    B = [t*u, t*u^2, u^2, u, s, r*s*u, r*t*u]

    test_lrcc_appendix(
        "C2 x C2 x S3",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=588,
        k=29,
        max_weight=16,
        dx_bound=18,
        dz_bound=12,
    )

    G = small_group(24, 14) # C2 x C2 x S3
    r, s, t, u = Oscar.gens(G)
    A = [s*t, r*s*t*u, r*t*u, s*u^2, s*u, r*s*t*u^2, r*s*t]
    B = [t*u, t*u^2, u^2, u, r*s*u^2, r, s]

    test_lrcc_appendix(
        "C2 x C2 x S3",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=588,
        k=39,
        max_weight=16,
        dx_bound=12,
        dz_bound=9,
    )

    G = small_group(24, 15) # C6 x C2 x C2
    r, s, t, u = Oscar.gens(G)
    A = [r*s*u^2, r*s*u, u, u^2, r*s, r*t, s*t]
    B = [s*t*u^2, s*t*u, r*s*t*u, r*s*t*u^2, s*u, s*u^2, r*s*t]

    test_lrcc_appendix(
        "C6 x C2 x C2",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=588,
        k=40,
        max_weight=16,
        dx_bound=9,
        dz_bound=18,
    )

    G = small_group(24, 15) # C6 x C2 x C2
    r, s, t, u = Oscar.gens(G)
    A = [s*u, s*u^2, r*t*u, r*t*u^2, r*s*t*u, r*s*t*u^2, r*t]
    B = [t*u^2, t*u, r*s*t, u, u^2, s*t*u^2, s*t*u]

    test_lrcc_appendix(
        "C6 x C2 x C2",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=588,
        k=38,
        max_weight=16,
        dx_bound=14,
        dz_bound=10,
    )

    G = small_group(24, 15) # C6 x C2 x C2
    r, s, t, u = Oscar.gens(G)
    A = [r*t, s*u, s*u^2, r*s, r*s*u^2, r*s*u, s]
    B = [u^2, u, r*s*t*u, r*s*t*u^2, s*t*u^2, s*t*u, s*t]

    test_lrcc_appendix(
        "C6 x C2 x C2",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=588,
        k=21,
        max_weight=16,
        dx_bound=21,
        dz_bound=12,
    )

    G = small_group(32, 3) # C8 x C4
    r, s, t, u, v = Oscar.gens(G)
    A = [s*t*u*v, s*t, t*u*v, t*u, r*s, r*s*t*u*v, v]
    B = [r*t, r*v, r*s*v, r*s*t*u, r*t*u, r*u*v, u*v]

    test_lrcc_appendix(
        "C8 x C4",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=784,
        k=16,
        max_weight=16,
        dx_bound=19,
        dz_bound=25,
    )

    G = small_group(32, 5) # (C8 x C2) : C2
    r, s, t, u, v = Oscar.gens(G)
    A = [r*t*u, r*t*v, t, t*u, t*u*v, r*s*t*u, r*s*v]
    B = [s*t*v, v, r*s*t*u*v, r*s, r*t, r*t*u*v, s*v]

    test_lrcc_appendix(
        "(C8 x C2) : C2",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=784,
        k=37,
        max_weight=16,
        dx_bound=14,
        dz_bound=12,
    )

    G = small_group(32, 6) # (C2 x C2 x C2) : C4
    r, s, t, u, v = Oscar.gens(G)
    A = [r*s*t, r*s*u*v, v, r*s, r*s*t*u, r*s*t*v, r*s*u]
    B = [r*t, r*t*u*v, r*u*v, r*v, u*v, t, s*t*v]

    test_lrcc_appendix(
        "(C2 x C2 x C2) : C4",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=784,
        k=40,
        max_weight=16,
        dx_bound=14,
        dz_bound=12,
    )

    G = small_group(32, 7) # (C8 : C2) : C2
    r, s, t, u, v = Oscar.gens(G)
    A = [r*s*t, r*s*u, v, r*s, r*s*t*u*v, r*s*t*v, r*s*u*v]
    B = [r*t, r*t*u, r*u*v, r, u*v, u, t]

    test_lrcc_appendix(
        "(C8 : C2) : C2",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=784,
        k=24,
        max_weight=16,
        dx_bound=16,
        dz_bound=13,
    )

    G = small_group(32, 7) # (C8 : C2) : C2
    r, s, t, u, v = Oscar.gens(G)
    A = [r*t*u, r*t, r*v, r*u, s, r*u*v, r]
    B = [s*u*v, t*v, r*s, r*s*t*u*v, s*t*u, t, v]

    test_lrcc_appendix(
        "(C8 : C2) : C2",
        G,
        A,
        B,
        ((H743, G743), (H734, G734));
        n=784,
        k=58,
        max_weight=16,
        dx_bound=12,
        dz_bound=10,
    )
end
