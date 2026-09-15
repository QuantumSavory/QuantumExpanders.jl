@testitem "Lifted QT manuscript permutation regression" begin
    using Test
    using Oscar
    using QuantumExpanders

    G = codomain(isomorphism(PermGroup, small_group(12, 1)))
    x = cperm(G, [5, 6, 7])
    y = cperm(G, [1, 4, 3, 2], [6, 7])
    A = [one(G), one(G), x, x, y, y * x^2, y * x, inv(y)]

    H = [
        1 0 0 0 0 1 1 1
        0 1 0 0 1 0 1 1
        0 0 1 0 1 1 0 1
        0 0 0 1 1 1 1 0
    ]

    # This is the permutation recorded for the [[768, 24]] code in the
    # corrected manuscript appendix and in the Table 1 construction data.
    corrected_p2 = [1, 2, 3, 4, 5, 6, 8, 7]
    corrected_code = QuantumTannerViaLeftRightActions(
        G,
        A,
        A,
        H,
        H;
        p2 = corrected_p2,
    )

    @test code_n(corrected_code) == 768
    @test code_k(corrected_code) == 24
    hx, hz = parity_matrix_xz(corrected_code)
    @test maximum(vec(sum(hx, dims=2))) == 16
    @test maximum(vec(sum(hz, dims=2))) == 16
    @test iszero(mod.(hx * hz', 2))

    # The old manuscript permutation produced k = 16 and must not be used for
    # the [[768, 24]] row.  Keep this assertion to document the discrepancy.
    old_manuscript_p2 = [1, 2, 3, 4, 5, 7, 8, 6]
    old_manuscript_code = QuantumTannerViaLeftRightActions(
        G,
        A,
        A,
        H,
        H;
        p2 = old_manuscript_p2,
    )
    @test code_k(old_manuscript_code) == 16

    # Distances are randomized estimator metadata, so they are deliberately
    # not asserted by this deterministic regression test.
end
