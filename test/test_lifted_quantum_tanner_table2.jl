@testitem "Lifted Quantum Tanner Codes — paper Table 2" begin
    using Test
    using Oscar
    using GAP
    using QECCore
    using QuantumExpanders
    import Nemo
    using Nemo: matrix, GF

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
        @assert all(iszero, H_nemo * null)
        null_t = Nemo.transpose(null)
        return Matrix{Int}([
            Int(lift(ZZ, null_t[i, j]))
            for i in 1:nrows(null_t), j in 1:ncols(null_t)
        ])
    end

    H633 = [
        1 0 0 0 1 1;
        0 1 0 1 0 1;
        0 0 1 1 1 0
    ]
    G633 = [
        0 1 1 1 0 0;
        1 0 1 0 1 0;
        1 1 0 0 0 1
    ]

    H844 = [
        1 0 0 0 0 1 1 1;
        0 1 0 0 1 0 1 1;
        0 0 1 0 1 1 0 1;
        0 0 0 1 1 1 1 0
    ]
    G844 = [
        0 1 1 1 1 0 0 0;
        1 0 1 1 0 1 0 0;
        1 1 0 1 0 0 1 0;
        1 1 1 0 0 0 0 1
    ]

    G734 = [
        1 0 1 1 1 0 0;
        1 1 1 0 0 1 0;
        0 1 1 1 0 0 1
    ]
    H734 = dual_code(G734)

    local_codes = Dict(
        :c633 => (H633, G633),
        :c844 => (H844, G844),
        :c734 => (H734, G734),
    )

    function paper_perm(G, spec)
        isempty(spec) && return one(G)
        cperm(G, [collect(cycle) for cycle in spec]...)
    end

    cases = [
        (
            name = "[[480, 8, (≤ 21, ≤ 21)]]",
            group = (10, 1), # D10
            ca = :c844,
            cb = :c633,
            A = [
                (), (),
                ((2,5),(3,4)), ((2,5),(3,4)),
                ((1,5,4,3,2),),
                ((1,5),(2,4)), ((1,5),(2,4)),
                ((1,4,2,5,3),),
            ],
            B = [
                (),
                ((2,5),(3,4)),
                ((1,5,4,3,2),),
                ((1,5),(2,4)),
                ((1,4),(2,3)),
                ((1,2,3,4,5),),
            ],
            p1 = 1:8,
            p2 = [1,2,5,6,3,4],
            n = 480,
            k = 8,
            wx = 12,
            wz = 12,
            dx = 21,
            dz = 21,
        ),
        (
            name = "[[504, 4, (≤ 36, ≤ 27)]]",
            group = (9, 2), # C3 x C3
            ca = :c844,
            cb = :c734,
            A = [
                (), (),
                ((4,6,5),), ((4,6,5),),
                ((1,3,2),),
                ((1,3,2),(4,6,5)),
                ((1,2,3),),
                ((1,2,3),(4,5,6)),
            ],
            B = [
                (), (),
                ((4,6,5),), ((4,6,5),),
                ((1,3,2),), ((1,3,2),),
                ((1,3,2),(4,6,5)),
            ],
            p1 = 1:8,
            p2 = [1,2,3,5,6,7,4],
            n = 504,
            k = 4,
            wx = 16,
            wz = 16,
            dx = 36,
            dz = 27,
        ),
        (
            name = "[[672, 4, (≤ 48, ≤ 28)]]",
            group = (12, 1), # C3 : C4
            ca = :c844,
            cb = :c734,
            A = [
                (), (),
                ((5,6,7),), ((5,6,7),),
                ((1,4,3,2),(6,7)),
                ((1,4,3,2),(5,7)),
                ((1,4,3,2),(5,6)),
                ((1,2,3,4),(6,7)),
            ],
            B = [
                (), (),
                ((5,6,7),), ((5,6,7),),
                ((1,4,3,2),(6,7)),
                ((1,4,3,2),(6,7)),
                ((1,4,3,2),(5,7)),
            ],
            p1 = 1:8,
            p2 = [1,2,3,5,6,7,4],
            n = 672,
            k = 4,
            wx = 16,
            wz = 16,
            dx = 48,
            dz = 28,
        ),
        (
            name = "[[720, 6, (≤ 30, ≤ 30)]]",
            group = (20, 3), # C5 : C4
            ca = :c633,
            cb = :c633,
            A = [
                (),
                ((2,4,5,3),),
                ((1,5,4,3,2),),
                ((1,5,2,3),),
                ((1,3,4,2),),
                ((1,2,3,4,5),),
            ],
            B = [
                (),
                ((2,4,5,3),),
                ((1,5,4,3,2),),
                ((1,5,2,3),),
                ((1,3,4,2),),
                ((1,2,3,4,5),),
            ],
            p1 = 1:6,
            p2 = [1,2,5,4,6,3],
            n = 720,
            k = 6,
            wx = 9,
            wz = 9,
            dx = 30,
            dz = 30,
        ),
        (
            name = "[[864, 8, (≤ 39, ≤ 31)]]",
            group = (18, 3), # C3 x S3
            ca = :c844,
            cb = :c633,
            A = [
                (), (),
                ((5,6),), ((5,6),),
                ((1,3,2),(4,5,6)),
                ((1,3,2),(4,5)),
                ((1,3,2),(4,5)),
                ((1,3,2),(4,6,5)),
            ],
            B = [
                (),
                ((5,6),),
                ((4,5,6),),
                ((1,3,2),(4,5)),
                ((1,2,3),(4,5)),
                ((1,3,2),(4,6)),
            ],
            p1 = 1:8,
            p2 = [1,2,4,5,3,6],
            n = 864,
            k = 8,
            wx = 12,
            wz = 12,
            dx = 39,
            dz = 31,
        ),
        (
            name = "[[864, 16, (≤ 36, ≤ 32)]]",
            group = (18, 3), # C3 x S3
            ca = :c844,
            cb = :c633,
            A = [
                (), (),
                ((5,6),), ((5,6),),
                ((1,3,2),(4,5,6)),
                ((1,3,2),(4,5)),
                ((1,3,2),(4,5)),
                ((1,3,2),(4,6,5)),
            ],
            B = [
                (),
                ((5,6),),
                ((4,5,6),),
                ((1,3,2),(4,5)),
                ((1,2,3),(4,5)),
                ((1,3,2),(4,6)),
            ],
            p1 = 1:8,
            p2 = [1,2,6,4,5,3],
            n = 864,
            k = 16,
            wx = 12,
            wz = 12,
            dx = 36,
            dz = 32,
        ),
        (
            name = "[[864, 8, (≤ 42, ≤ 40)]]",
            group = (18, 3), # C3 x S3
            ca = :c844,
            cb = :c633,
            A = [
                (), (),
                ((5,6),), ((5,6),),
                ((1,3,2),(4,5,6)),
                ((1,3,2),(4,5)),
                ((1,3,2),(4,5)),
                ((1,3,2),(4,6,5)),
            ],
            B = [
                (),
                ((5,6),),
                ((4,5,6),),
                ((1,3,2),(4,5)),
                ((1,2,3),(4,5)),
                ((1,3,2),(4,6)),
            ],
            p1 = 1:8,
            p2 = [1,2,6,5,4,3],
            n = 864,
            k = 8,
            wx = 12,
            wz = 12,
            dx = 42,
            dz = 40,
        ),
        (
            name = "[[1120, 4, (≤ 80, ≤ 35)]]",
            group = (20, 1), # C5 : C4
            ca = :c844,
            cb = :c734,
            A = [
                (), (),
                ((5,9,8,7,6),), ((5,9,8,7,6),),
                ((1,4,3,2),(6,9),(7,8)),
                ((1,4,3,2),(5,6),(7,9)),
                ((1,4,3,2),(5,7),(8,9)),
                ((1,4,3,2),(5,8),(6,7)),
            ],
            B = [
                (), (),
                ((5,9,8,7,6),), ((5,9,8,7,6),),
                ((1,4,3,2),(6,9),(7,8)),
                ((1,4,3,2),(6,9),(7,8)),
                ((1,4,3,2),(5,6),(7,9)),
            ],
            p1 = 1:8,
            p2 = [1,2,3,6,4,7,5],
            n = 1120,
            k = 4,
            wx = 16,
            wz = 16,
            dx = 80,
            dz = 35,
        ),
    ]
    group_cache = Dict{Tuple{Int,Int}, Any}()
    for case in cases
        @testset "$(case.name)" begin
            group_order, group_id = case.group
            G = get!(group_cache, case.group) do
                codomain(isomorphism(
                    PermGroup,
                    small_group(group_order, group_id),
                ))
            end
            A = [paper_perm(G, p) for p in case.A]
            B = [paper_perm(G, p) for p in case.B]
            H_A, G_A = local_codes[case.ca]
            H_B, G_B = local_codes[case.cb]
            c = QuantumTannerViaLeftRightActions(
                G,
                A,
                B,
                H_A,
                G_A,
                H_B,
                G_B;
                p1 = collect(case.p1),
                p2 = collect(case.p2),
            )
            hx, hz = parity_matrix_xz(c)
            @test code_n(c) == case.n
            @test code_k(c) == case.k
            @test maximum(vec(sum(hx, dims=2))) == case.wx
            @test maximum(vec(sum(hz, dims=2))) == case.wz
            dx_qdist, dz_qdist = compute_distance(hx, hz; num=10_000)
            println()
            println(case.name)
            println("  QDistRnd 50K : (dx, dz) = ($dx_qdist, $dz_qdist)")
            println("  sQetch paper : (dx, dz) = (≤$(case.dx), ≤$(case.dz))")
        end
    end
end
