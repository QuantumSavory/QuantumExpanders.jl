using Nemo
using Oscar
using LinearAlgebra
using Random

using QuantumClifford

function uniformly_random_code_checkmatrix(ρ, Δ)
    @assert 0<ρ<1
    r = Int(floor(ρ*Δ))
    uniformly_random_code_checkmatrix(r, Δ)
end

function uniformly_random_code_checkmatrix(r::Integer, Δ)
    space = MatrixSpace(ResidueRing(ZZ,2), r, Δ)
    local H
    while true
        H = rand(space)
        r == rank(H) && break
    end
    H
end

function dual_code(H::nmod_mat)
    null = nullspace(H)[2]
    @assert all(H*null .== 0)
    #@assert size(H,1) + size(null,2) == size(H,2)
    transpose(null)
end

function dual_code(H)
    r, Δ = size(H)
    ring = ResidueRing(ZZ,2)
    H = MatrixSpace(ring, r, Δ)(ring.(H)) # TODO there must be a cleaner way to write this
    dual_code(H)
end

using LinearAlgebra
function LinearAlgebra.kron(l::nmod_mat,r::nmod_mat)
    l1,l2 = size(l)
    r1,r2 = size(r)
    s1 = l1*r1
    s2 = l2*r2
    MatrixSpace(l.base_ring,s1,s2)(kron(Matrix(l),Matrix(r)))
end

"""Check that two binary parity check matrices X and Z result in a good CSS code
(i.e., commutation constraints are fulfilled)"""
function good_css(X::nmod_mat,Z::nmod_mat)
    x = (x->Bool(x.data)).(X)
    z = (x->Bool(x.data)).(Z)
    good_css(x,z)
end

function good_css(X::nmod_mat,Z)
    x = (x->Bool(x.data)).(X)
    good_css(x,Z)
end

function good_css(X,Z::nmod_mat)
    z = (x->Bool(x.data)).(Z)
    good_css(X,z)
end

function good_css(X,Z)
    stab = css(X,Z)
    good_css(stab)
end

function good_css(stab)
    for row in stab
        !all(==(0), QuantumClifford.comm(row,stab)) && return false
    end
    return true
end

"""Create a CSS code from two binary parity check matrices, X and Z"""
css(X,Z) = vcat(Stabilizer(zero(X),X),Stabilizer(Z,zero(Z)))
