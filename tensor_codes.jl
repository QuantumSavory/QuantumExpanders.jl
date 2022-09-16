using Nemo
using Oscar
using LinearAlgebra
using Random

function uniformly_random_CᴬCᴮ(ρ, Δ)
    @assert 0<ρ<1
    ρᴬ = ρ
    ρᴮ = 1-ρ
    rᴬ = Int(floor(ρᴬ*Δ))
    rᴮ = Δ-rᴬ
    spaceᴬ = MatrixSpace(ResidueRing(ZZ,2), rᴬ, Δ)
    spaceᴮ = MatrixSpace(ResidueRing(ZZ,2), rᴮ, Δ)
    local Hᴬ, Hᴮ
    while true
        Hᴬ = rand(spaceᴬ)
        rᴬ == rank(Hᴬ) && break
    end
    while true
        Hᴮ = rand(spaceᴮ)
        rᴮ == rank(Hᴮ) && break
    end
    Hᴬ, Hᴮ
end

function uniformly_random_code(ρ, Δ)
    @assert 0<ρ<1
    r = Int(floor(ρ*Δ))
    space = MatrixSpace(ResidueRing(ZZ,2), r, Δ)
    local H
    while true
        H = rand(space)
        r == rank(H) && break
    end
    H
end

function dual_code(H)
    #r, Δ = size(H)
    #H = MatrixSpace(ResidueRing(ZZ,2), r, Δ)(H)
    null = nullspace(H)[2]
    @assert all(H*null .== 0)
    @assert size(H,1) + size(null,2) == size(H,2)
    transpose(null)
end

using LinearAlgebra
function LinearAlgebra.kron(l::nmod_mat,r::nmod_mat)
    l1,l2 = size(l)
    r1,r2 = size(r)
    s1 = l1*r1
    s2 = l2*r2
    MatrixSpace(l.base_ring,s1,s2)(kron(Matrix(l),Matrix(r)))
end
