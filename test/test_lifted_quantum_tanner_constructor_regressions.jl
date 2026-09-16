@testitem "Lifted QT constructors accept Oscar matrix types" begin
    using Test
    using Oscar
    using QuantumExpanders
    using QuantumClifford.ECC
    using QECCore

    G = cyclic_group(1)
    A = fill(one(G), 2)

    # Regression for the Oscar FqFieldElem conversion failure reported in
    # review.  The [1 1] repetition code is self-dual, so the same matrix can
    # be supplied as both its parity-check and generator matrix.
    H_fq = matrix(GF(2), [1 1])
    c_fq = QuantumTannerViaLeftRightActions(
        G,
        A,
        A,
        H_fq,
        H_fq,
        H_fq,
        H_fq,
    )

    # Regression for the parity-only constructor rejecting the zzModMatrix
    # returned by dual_code.
    H_zz = dual_code([1 1])
    c_zz = QuantumTannerViaLeftRightActions(G, A, A, H_zz, H_zz)

    @test code_n(c_fq) == 4
    @test code_k(c_fq) == 2
    @test code_n(c_zz) == 4
    @test code_k(c_zz) == 2

    hx_fq, hz_fq = parity_matrix_xz(c_fq)
    hx_zz, hz_zz = parity_matrix_xz(c_zz)
    @test hx_fq == hx_zz
    @test hz_fq == hz_zz
    @test iszero(mod.(hx_fq * hz_fq', 2))
    @test iszero(mod.(hx_zz * hz_zz', 2))
end
