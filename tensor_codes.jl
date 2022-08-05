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

function dual_code(H)
    null = nullspace(H)[2]
    @assert all(a*null .== 0)
    @assert size(a,1) + size(null,2) == size(a,2)
    transpose(null)
end
