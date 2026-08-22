function _lambda(a, idx)
    n = length(idx)
    rows = Int[]
    cols = Int[]
    for (g, i) in idx
        push!(rows, i)
        push!(cols, idx[a*g])
    end
    sparse(rows, cols, ones(Int, n), n, n)
end

function _rho(b, idx)
    n = length(idx); b_inv = inv(b)
    rows = Int[]
    cols = Int[]
    for (g, i) in idx
        push!(rows, i)
        push!(cols, idx[g*b_inv])
    end
    sparse(rows, cols, ones(Int, n), n, n)
end

function build_LA_RB(group, A, B)
    elems = collect(group)
    n_G = length(elems)
    idx   = Dict(g => k for (k,g) in enumerate(elems))
    blocks_LA = [_lambda(A[i], idx) for i in 1:length(A) for j in 1:length(B)]
    blocks_RB = [_rho(B[j],   idx) for i in 1:length(A) for j in 1:length(B)]
    blockdiag(blocks_LA...), blockdiag(blocks_RB...)
end

function parity_matrices(LA, RB, n_G, H0, H1, G0, G1, H0p, H1p, G0p, G1p)
    I_G = sparse(1:n_G, 1:n_G, ones(Int, n_G))
    HX  = mod.(vcat(
        Matrix(kron(kron(sparse(H0),  sparse(G0p)), I_G)),
        Matrix(kron(kron(sparse(H1),  sparse(G1p)), I_G)*LA*RB)), 2)
    HZ  = mod.(vcat(
        Matrix(kron(kron(sparse(G0),  sparse(H1p)), I_G)*RB),
        Matrix(kron(kron(sparse(G1),  sparse(H0p)), I_G)*LA)), 2)
    HX, HZ
end

struct QuantumTannerViaLeftRightActions{GT, ET}
    group::GT
    A::Vector{ET}
    B::Vector{ET}
    H0_A::Matrix{Int}
    H1_A::Matrix{Int}
    G0_A::Matrix{Int}
    G1_A::Matrix{Int}
    H0_B::Matrix{Int}
    H1_B::Matrix{Int}
    G0_B::Matrix{Int}
    G1_B::Matrix{Int}
    hx::Matrix{Int}
    hz::Matrix{Int}

    function QuantumTannerViaLeftRightActions(group, A::Vector, B::Vector,
            H0_A, H1_A, G0_A, G1_A, H0_B, H1_B, G0_B, G1_B)
        nA = size(H0_A, 2)
        nB = size(H0_B, 2)
        length(A) == nA || throw(ArgumentError("|A| = $(length(A)) ≠ n_A = $nA"))
        length(B) == nB || throw(ArgumentError("|B| = $(length(B)) ≠ n_B = $nB"))
        _check_duality(H0_A, G0_A, "A-side (H0,G0)")
        _check_duality(H1_A, G1_A, "A-side (H1,G1)")
        _check_duality(H0_B, G0_B, "B-side (H0',G0')")
        _check_duality(H1_B, G1_B, "B-side (H1',G1')")
        n_G    = length(collect(group))
        LA, RB = build_LA_RB(group, A, B)
        hx, hz = parity_matrices(LA, RB, n_G, H0_A, H1_A, G0_A, G1_A, H0_B, H1_B, G0_B, G1_B)
        iszero(mod.(hx * hz', 2)) || throw(ArgumentError("internal error: HX·HZᵀ ≠ 0 (CSS orthogonality violated)"))
        new{typeof(group), eltype(A)}(group, A, B,
            Matrix{Int}(H0_A), Matrix{Int}(H1_A),
            Matrix{Int}(G0_A), Matrix{Int}(G1_A),
            Matrix{Int}(H0_B), Matrix{Int}(H1_B),
            Matrix{Int}(G0_B), Matrix{Int}(G1_B), hx, hz)
    end
end

function _check_duality(H, G, label)
    F2 = GF(2)
    Hm = matrix(F2, Matrix{Int}(H))
    Gm = matrix(F2, Matrix{Int}(G))
    iszero(Hm * transpose(Gm)) || throw(ArgumentError("$label: H·Gᵀ ≠ 0 over F_2"))
    rank(Gm) == size(H,2) - rank(Hm) || throw(ArgumentError("$label: rows of G do not span ker H " * "(rank G = $(rank(Gm)), n − rank H = $(size(H,2)-rank(Hm)))"))
end

function QuantumTannerViaLeftRightActions(group, A::Vector, B::Vector, H_A, G_A, H_B, G_B; p1 = collect(1:size(H_A,2)), p2 = collect(1:size(H_B,2)))
    QuantumTannerViaLeftRightActions(group, A, B, H_A, H_A[:,p1], G_A, G_A[:,p1], H_B, H_B[:,p2], G_B, G_B[:,p2])
end

function QuantumTannerViaLeftRightActions(group, A::Vector, B::Vector, H_A::Matrix, H_B::Matrix; p1 = collect(1:size(H_A,2)), p2 = collect(1:size(H_B,2)))
    QuantumTannerViaLeftRightActions(group, A, B, H_A, dual_code(H_A), H_B, dual_code(H_B); p1 = p1, p2 = p2)
end
