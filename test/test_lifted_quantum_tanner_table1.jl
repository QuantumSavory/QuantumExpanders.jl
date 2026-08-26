@testitem "Lifted Quantum Tanner Codes. Paper Table 1" begin
    using Test
    using Oscar
    using GAP
    using QECCore
    using QuantumExpanders
    import Nemo
    using Nemo: matrix, GF
    using QuantumClifford
    using QuantumClifford.ECC

    loaded = GAP.Globals.LoadPackage(GAP.GapObj("QDistRnd"))

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

    function compute_distance(hx, hz; num=50_000)
        @assert iszero(mod.(hx * hz', 2))
        GX = julia_to_gap_gf2(hx)
        GZ = julia_to_gap_gf2(hz)
        dz = GAP.Globals.DistRandCSS(
            GX, GZ, GAP.Obj(num), GAP.Obj(0), GAP.Obj(0)
        )
        dx = GAP.Globals.DistRandCSS(
            GZ, GX, GAP.Obj(num), GAP.Obj(0), GAP.Obj(0)
        )
        return Int(dx), Int(dz)
    end

    function dual_code(H)
        H_nemo = matrix(GF(2), H)
        null = Nemo.nullspace(H_nemo)[2]
        @assert all(iszero, H_nemo*null)
        null_t = Nemo.transpose(null)
        Matrix{Int}([Int(lift(ZZ, null_t[i, j]))
                     for i in 1:nrows(null_t), j in 1:ncols(null_t)])
    end

    H633 = [1 0 0 0 1 1;
            0 1 0 1 0 1;
            0 0 1 1 1 0]

    G633 = [0 1 1 1 0 0;
            1 0 1 0 1 0;
            1 1 0 0 0 1]

    H844 = [1 0 0 0 0 1 1 1;
            0 1 0 0 1 0 1 1;
            0 0 1 0 1 1 0 1;
            0 0 0 1 1 1 1 0]

    G844 = [0 1 1 1 1 0 0 0;
            1 0 1 1 0 1 0 0;
            1 1 0 1 0 0 1 0;
            1 1 1 0 0 0 0 1]

    G734 = [1 0 1 1 1 0 0;
            1 1 1 0 0 1 0;
            0 1 1 1 0 0 1]
    H734 = dual_code(G734)

    G743 = [1 0 0 0 1 1 0;
            0 1 0 0 1 0 1;
            0 0 1 0 0 1 1;
            0 0 0 1 1 1 1]
    H743 = dual_code(G743)

    G953 = [1 0 0 0 0 1 1 1 1;
            0 1 0 0 0 1 1 1 0;
            0 0 1 0 0 1 1 0 1;
            0 0 0 1 0 1 0 1 1;
            0 0 0 0 1 0 1 1 1]

    H953 = dual_code(G953)
    local_codes = Dict(
        :c633 => (H633, G633),
        :c844 => (H844, G844),
        :c734 => (H734, G734),
        :c743 => (H743, G743),
        :c953 => (H953, G953),
    )

    function paper_perm(S, spec)
        isempty(spec) && return one(S)
        cperm(S, [collect(cycle) for cycle in spec]...)
    end

    GROUP_IDS = Dict(
        :S3      => (6, 1),  # describe(small_group(6,1))
        :D10     => (10, 1), # describe(small_group(10,1))
        :C3xC3   => (9, 2),  # describe(small_group(9,2))
        :C3sdC4  => (12, 1), # describe(small_group(12,1))
        :D14     => (14, 1), # describe(small_group(14,1))
        :C3xS3   => (18, 3), # describe(small_group(18,3))
        :A4      => (12, 3), # describe(small_group(12,3))
        :C4sdC4  => (16, 4), # describe(small_group(16,4))
        :C5sdC4  => (20, 1), # describe(small_group(20,1))
    )

    function build_lift(spec)
        n, id = GROUP_IDS[spec.group]
        G = codomain(isomorphism(PermGroup, small_group(n, id)))
        A = [paper_perm(G, p) for p in spec.A]
        B = [paper_perm(G, p) for p in spec.B]
        return G, A, B
    end

    A_s3_8 = [
        (), (), ((2,3),), ((2,3),),
        ((1,2,3),), ((1,2),), ((1,2),), ((1,3,2),)
    ]

    B_s3_6 = [
        (), (), ((2,3),), ((1,2,3),), ((1,2),), ((1,3,2),)
    ]

    B_s3_7 = [
        (), (), ((2,3),), ((2,3),),
        ((1,2,3),), ((1,2,3),), ((1,2),)
    ]

    A_d10_6 = [
        (), ((2,5),(3,4)), ((1,5,4,3,2),),
        ((1,5),(2,4)), ((1,4),(2,3)), ((1,2,3,4,5),)
    ]
    A_d10_8 = [
        (), (), ((2,5),(3,4)), ((2,5),(3,4)),
        ((1,5,4,3,2),), ((1,5),(2,4)), ((1,5),(2,4)),
        ((1,4,2,5,3),)
    ]

    A_c3xc3_8 = [
        (), (), ((4,6,5),), ((4,6,5),), ((1,3,2),),
        ((1,3,2),(4,6,5)), ((1,2,3),), ((1,2,3),(4,5,6))
    ]
    B_c3xc3_6 = [
        (), ((4,6,5),), ((4,5,6),), ((1,3,2),),
        ((1,3,2),(4,6,5)), ((1,2,3),)
    ]

    A_c3sc4_7 = [
        (), (), ((5,6,7),), ((5,6,7),),
        ((1,4,3,2),(6,7)), ((1,4,3,2),(6,7)),
        ((1,4,3,2),(5,7))
    ]
    B_c3sc4_6 = [
        (), ((1,3),(2,4),(5,6,7)), ((1,4,3,2),(6,7)),
        ((1,4,3,2),(5,7)), ((1,4,3,2),(5,6)),
        ((1,2,3,4),(6,7))
    ]

    A_d14_7 = [
        (), ((2,7),(3,6),(4,5)), ((1,5,2,6,3,7,4),),
        ((1,5),(2,4),(6,7)), ((1,2),(3,7),(4,6)),
        ((1,6,4,2,7,5,3),), ((1,3),(4,7),(5,6))
    ]
    B_d14_6 = A_d14_7[1:6]

    A_c3xs3_6 = [
        (), ((5,6),), ((4,5,6),), ((1,3,2),(4,5)),
        ((1,2,3),(4,5)), ((1,3,2),(4,6))
    ]

    A_756_1 = A_c3sc4_7
    A_756_2 = [
        (), (), ((1,3),(2,4),(5,6,7)), ((1,3),(2,4),(5,7,6)),
        ((1,4,3,2),(6,7)), ((1,4,3,2),(6,7)),
        ((1,4,3,2),(5,7))
    ]
    B_756_9 = [
        (), (), ((5,6,7),), ((5,6,7),),
        ((1,4,3,2),(6,7)), ((1,4,3,2),(5,7)),
        ((1,4,3,2),(5,6)), ((1,2,3,4),(6,7)),
        ((1,2,3,4),(5,7))
    ]

    A_c3sc4_8 = [
        (), (), ((5,6,7),), ((5,6,7),),
        ((1,4,3,2),(6,7)), ((1,4,3,2),(5,7)),
        ((1,4,3,2),(5,6)), ((1,2,3,4),(6,7))
    ]

    A_a4_8 = [
        (), ((2,4,3),), ((2,4,3),), ((1,4),(2,3)),
        ((1,4,2),), ((1,4,3),), ((1,2,3),), ((1,2,3),)
    ]
    B_a4_8 = [
        (), ((2,4,3),), ((2,3,4),), ((1,4),(2,3)),
        ((1,4,2),), ((1,4,3),), ((1,2),(3,4)), ((1,3,2),)
    ]

    A_c3xs3_8 = [
        (), (), ((5,6),), ((5,6),), ((1,3,2),(4,5,6)),
        ((1,3,2),(4,5)), ((1,3,2),(4,5)), ((1,3,2),(4,6,5))
    ]

    A_c4sc4_8 = [
        ((5,8,7,6),), ((1,3),(2,4)), ((1,3),(2,4),(5,7),(6,8)),
        ((1,3),(2,4),(5,6,7,8)), ((1,4,3,2),(6,8)),
        ((1,4,3,2),(5,6),(7,8)), ((1,2,3,4),(6,8)),
        ((1,2,3,4),(5,6),(7,8))
    ]
    B_c4sc4_7 = [
        (), (), ((5,8,7,6),), ((5,8,7,6),),
        ((1,4,3,2),(6,8)), ((1,4,3,2),(5,6),(7,8)),
        ((1,2,3,4),(5,6),(7,8))
    ]
    B_c4sc4_8 = [
        (), (), ((5,8,7,6),), ((5,6,7,8),),
        ((1,4,3,2),(6,8)), ((1,4,3,2),(5,6),(7,8)),
        ((1,2,3,4),(6,8)), ((1,2,3,4),(5,6),(7,8))
    ]

    A_c5sc4_8 = [
        (), (), ((5,9,8,7,6),), ((5,9,8,7,6),),
        ((1,4,3,2),(6,9),(7,8)), ((1,4,3,2),(5,6),(7,9)),
        ((1,4,3,2),(5,7),(8,9)), ((1,4,3,2),(5,8),(6,7))
    ]

    lift_specs = Dict(
        :s3_8x6       => (group=:S3,      A=A_s3_8,      B=B_s3_6),
        :s3_8x7       => (group=:S3,      A=A_s3_8,      B=B_s3_7),
        :d10_6x6      => (group=:D10,     A=A_d10_6,     B=A_d10_6),
        :s3_8x8       => (group=:S3,      A=A_s3_8,      B=A_s3_8),
        :c3xc3_8x6    => (group=:C3xC3,   A=A_c3xc3_8,   B=B_c3xc3_6),
        :d10_8x6      => (group=:D10,     A=A_d10_8,     B=A_d10_6),
        :c3sc4_7x6    => (group=:C3sdC4,  A=A_c3sc4_7,   B=B_c3sc4_6),
        :d14_7x6      => (group=:D14,     A=A_d14_7,     B=B_d14_6),
        :d10_8x8      => (group=:D10,     A=A_d10_8,     B=A_d10_8),
        :c3xs3_6x6    => (group=:C3xS3,   A=A_c3xs3_6,   B=A_c3xs3_6),
        :c3sc4_7x9_a1 => (group=:C3sdC4,  A=A_756_1,     B=B_756_9),
        :c3sc4_7x9_a2 => (group=:C3sdC4,  A=A_756_2,     B=B_756_9),
        :c3sc4_8x8    => (group=:C3sdC4,  A=A_c3sc4_8,   B=A_c3sc4_8),
        :a4_8x8       => (group=:A4,      A=A_a4_8,      B=B_a4_8),
        :c3xs3_8x6    => (group=:C3xS3,   A=A_c3xs3_8,   B=A_c3xs3_6),
        :c4sc4_8x7    => (group=:C4sdC4,  A=A_c4sc4_8,   B=B_c4sc4_7),
        :c4sc4_8x8    => (group=:C4sdC4,  A=A_c4sc4_8,   B=B_c4sc4_8),
        :c5sc4_8x8    => (group=:C5sdC4,  A=A_c5sc4_8,   B=A_c5sc4_8),
    )

    lift_cache = Dict{Symbol, Any}()
    function get_lift(key)
        get!(lift_cache, key) do
            build_lift(lift_specs[key])
        end
    end
    cases = [
        (name="[[288, 8, (≤ 15, ≤ 15)]]",    lift=:s3_8x6,       ca=:c844, cb=:c633, p1=1:8, p2=[1,2,4,3,6,5],       n=288,  k=8,  wx=12, wz=12, dx=15, dz=15),
        (name="[[336, 4, (≤ 24, ≤ 15)]]",    lift=:s3_8x7,       ca=:c844, cb=:c734, p1=1:8, p2=[1,2,3,5,6,7,4],     n=336,  k=4,  wx=16, wz=16, dx=24, dz=15),
        (name="[[336, 12, (≤ 14, ≤ 14)]]",   lift=:s3_8x7,       ca=:c844, cb=:c734, p1=1:8, p2=[1,2,3,5,4,7,6],     n=336,  k=12, wx=16, wz=16, dx=14, dz=14),
        (name="[[360, 6, (≤ 18, ≤ 18)]]",    lift=:d10_6x6,      ca=:c633, cb=:c633, p1=1:6, p2=[1,2,4,5,3,6],       n=360,  k=6,  wx=9,  wz=9, dx=18, dz=18),
        (name="[[384, 24, (≤ 14, ≤ 14)]]",   lift=:s3_8x8,       ca=:c844, cb=:c844, p1=1:8, p2=[1,2,3,4,5,6,8,7],   n=384,  k=24, wx=16, wz=16, dx=14, dz=14),
        (name="[[384, 8, (≤ 16, ≤ 16)]]",    lift=:s3_8x8,       ca=:c844, cb=:c844, p1=1:8, p2=[1,2,3,4,6,7,8,5],   n=384,  k=8,  wx=16, wz=16, dx=16, dz=16),
        (name="[[432, 16, (≤ 15, ≤ 15)]]",   lift=:c3xc3_8x6,    ca=:c844, cb=:c633, p1=1:8, p2=[1,2,4,3,5,6],       n=432,  k=16, wx=12, wz=12, dx=15, dz=15),
        (name="[[432, 8, (≤ 18, ≤ 18)]]",    lift=:c3xc3_8x6,    ca=:c844, cb=:c633, p1=1:8, p2=[1,2,4,6,5,3],       n=432,  k=8,  wx=12, wz=12, dx=18, dz=18),
        (name="[[480, 8, (≤ 21, ≤ 18)]]",    lift=:d10_8x6,      ca=:c844, cb=:c633, p1=1:8, p2=[1,2,3,6,4,5],       n=480,  k=8,  wx=12, wz=12, dx=21, dz=18),
        (name="[[504, 4, (≤ 34, ≤ 12)]]",    lift=:c3sc4_7x6,    ca=:c734, cb=:c633, p1=1:7, p2=[1,2,6,3,4,5],       n=504,  k=4,  wx=12, wz=12, dx=34, dz=12),
        (name="[[588, 6, (≤ 28, ≤ 23)]]",    lift=:d14_7x6,      ca=:c743, cb=:c633, p1=1:7, p2=[1,2,4,5,6,3],       n=588,  k=6,  wx=12, wz=12, dx=28, dz=23),
        (name="[[588, 14, (≤ 20, ≤ 18)]]",   lift=:d14_7x6,      ca=:c743, cb=:c633, p1=1:7, p2=[1,2,3,6,4,5],       n=588,  k=14, wx=12, wz=12, dx=20, dz=18),
        (name="[[640, 24, (≤ 14, ≤ 14)]]",   lift=:d10_8x8,      ca=:c844, cb=:c844, p1=1:8, p2=[1,2,3,4,5,6,8,7],   n=640,  k=24, wx=16, wz=16, dx=14, dz=14),
        (name="[[648, 26, (≤ 12, ≤ 15)]]",   lift=:c3xs3_6x6,    ca=:c633, cb=:c633, p1=1:6, p2=[1,2,3,4,6,5],       n=648,  k=26, wx=9,  wz=9, dx=12, dz=15),
        (name="[[648, 14, (≤ 16, ≤ 15)]]",   lift=:c3xs3_6x6,    ca=:c633, cb=:c633, p1=1:6, p2=[1,2,4,3,5,6],       n=648,  k=14, wx=9,  wz=9, dx=16, dz=15),
        (name="[[648, 2, (≤ 18, ≤ 32)]]",    lift=:c3xs3_6x6,    ca=:c633, cb=:c633, p1=1:6, p2=[1,2,4,6,3,5],       n=648,  k=2,  wx=9,  wz=9, dx=18, dz=32),
        (name="[[648, 4, (≤ 26, ≤ 21)]]",    lift=:c3xs3_6x6,    ca=:c633, cb=:c633, p1=1:6, p2=[1,2,6,3,4,5],       n=648,  k=4,  wx=9,  wz=9, dx=26, dz=21),
        (name="[[648, 10, (≤ 18, ≤ 20)]]",   lift=:c3xs3_6x6,    ca=:c633, cb=:c633, p1=1:6, p2=[1,2,5,6,4,3],       n=648,  k=10, wx=9,  wz=9, dx=18, dz=20),
        (name="[[756, 17, (≤ 9, ≤ 21]]",   lift=:c3sc4_7x9_a1, ca=:c734, cb=:c953, p1=1:7, p2=[1,2,3,4,5,6,8,9,7], n=756, k=17, wx=20, wz=20, dx=9, dz=21),
        (name="[[756, 10, (≤ 9, ≤ 42]]",   lift=:c3sc4_7x9_a1, ca=:c734, cb=:c953, p1=1:7, p2=[1,2,3,4,5,7,8,9,6], n=756, k=10, wx=20, wz=20, dx=9, dz=42),
        (name="[[756, 43, (≤ 12, ≤ 12)]]",   lift=:c3sc4_7x9_a2, ca=:c743, cb=:c953, p1=1:7, p2=1:9,                   n=756,  k=43, wx=20, wz=20, dx=12, dz=12),
        (name="[[756, 30, (≤ 21, ≤ 12)]]",   lift=:c3sc4_7x9_a2, ca=:c743, cb=:c953, p1=1:7, p2=[1,2,3,4,5,7,6,9,8], n=756, k=30, wx=20, wz=20, dx=21, dz=12),
        (name="[[768, 24, (≤ 16, ≤ 16)]]",   lift=:c3sc4_8x8,    ca=:c844, cb=:c844, p1=1:8, p2=[1,2,3,4,5,6,8,7],   n=768,  k=24, wx=16, wz=16, dx=18, dz=18),
        (name="[[768, 32, (≤ 16, ≤ 16)]]",   lift=:a4_8x8,       ca=:c844, cb=:c844, p1=1:8, p2=1:8,                   n=768,  k=32, wx=16, wz=16, dx=16, dz=16),
        (name="[[864, 24, (≤ 21, ≤ 18)]]",   lift=:c3xs3_8x6,    ca=:c844, cb=:c633, p1=1:8, p2=1:6,                   n=864,  k=24, wx=12, wz=12, dx=21, dz=18),
        (name="[[864, 16, (≤ 27, ≤ 18)]]",   lift=:c3xs3_8x6,    ca=:c844, cb=:c633, p1=1:8, p2=[1,2,3,5,4,6],       n=864,  k=16, wx=12, wz=12, dx=27, dz=18),
        (name="[[864, 8, (≤ 27, ≤ 32)]] A",  lift=:c3xs3_8x6,    ca=:c844, cb=:c633, p1=1:8, p2=[1,2,3,5,6,4],       n=864,  k=8,  wx=12, wz=12, dx=27, dz=32),
        (name="[[864, 8, (≤ 27, ≤ 47)]] B",  lift=:c3xs3_8x6,    ca=:c844, cb=:c633, p1=1:8, p2=[1,2,3,6,4,5],       n=864,  k=8,  wx=12, wz=12, dx=27, dz=47),
        (name="[[896, 40, (≤ 16, ≤ 12)]]",   lift=:c4sc4_8x7,    ca=:c844, cb=:c734, p1=1:8, p2=1:7,                   n=896,  k=40, wx=16, wz=16, dx=16, dz=12),
        (name="[[1280, 24, (≤ 16, ≤ 16)]]",  lift=:c5sc4_8x8,    ca=:c844, cb=:c844, p1=1:8, p2=[1,2,3,4,5,6,8,7],   n=1280, k=24, wx=16, wz=16, dx=16, dz=16),
    ]

    for case in cases
        @testset "$(case.name)" begin
            G, A, B = get_lift(case.lift)
            H_A, G_A = local_codes[case.ca]
            H_B, G_B = local_codes[case.cb]
            c = QuantumTannerViaLeftRightActions(G, A, B, H_A, G_A, H_B, G_B; p1 = collect(case.p1), p2 = collect(case.p2),)
            hx, hz = parity_matrix_xz(c)
            stab = QuantumClifford.ECC.parity_checks(c)
            mat = matrix(GF(2), stab_to_gf2(stab))
            computed_rank = rank(mat)
            @test computed_rank == code_n(c) - code_k(c)
            @test code_n(c) == case.n
            @test code_k(c) == case.k
            @test maximum(vec(sum(hx, dims=2))) == case.wx
            @test maximum(vec(sum(hz, dims=2))) == case.wz
        end
    end
end
