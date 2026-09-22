function _check_allrowscommute(Hx, Hz)
    for rowx in eachrow(Hx)
        for rowz in eachrow(Hz)
            comm = sum(rowx .& rowz)
            isodd(comm) && return false
        end
    end
    return true
end

"""
    random_quantum_Tanner_code(ρ, group, A, B; bipartite=true,
                               use_same_local_code=false,
                               stabilizer_convention=:leverrier_zemor,
                               rng=GLOBAL_RNG)

Generate random local codes and return `(H_X, H_Z)` for the resulting quantum
Tanner code.

By default, the X/Z assignment follows [leverrier2022quantum](@cite):
`C₀ = C_A ⊗ C_B` on `V₀` gives Z checks and
`C₁ = C_A^⊥ ⊗ C_B^⊥` on `V₁` gives X checks. This agrees with
[`QuantumTannerCode`](@ref). Set `stabilizer_convention=:gu` to use the
opposite assignment from [gu2022efficient](@cite).
"""
function random_quantum_Tanner_code(
    ρ::Real,
    group::Group,
    A::Vector{<:GroupElem},
    B::Vector{<:GroupElem};
    bipartite=true,
    use_same_local_code=false,
    stabilizer_convention::Symbol=:leverrier_zemor,
    rng::AbstractRNG=GLOBAL_RNG,
)
    stabilizer_convention in (:leverrier_zemor, :gu) ||
        throw(ArgumentError("stabilizer_convention must be :leverrier_zemor or :gu"))
    @show length(group), length(A), length(B)
    @show length(group)*length(A)*length(B)
    if bipartite
        @assert is_nonconjugate(group, A, B)
        @assert is_symmetric_gen(A)
        @assert is_symmetric_gen(B)
    end
    Δᴬ = length(A)
    Δᴮ = length(B)
    𝒢₀□, 𝒢₁□, edge₀_q_idx, edge₁_q_idx, edge₀_ab_idx, edge₁_ab_idx = bipartite ? cayley_complex_square_graphs(group, A, B) : cayley_complex_square_graphs_quadripartite(group, A, B)
    # "Let C_A, C_B ⊆ 𝔽₂^Δ be classical codes of rates ρ and (1-ρ) respectively" [gu2022efficient](@cite).
    # C_A has rate ρ, C_B has rate 1-ρ. For a code C with rate r = k/n, its dual C^⊥ has rate 1-r.
    # We consider H_A (parity check for C_A) defines C_A^⊥ which has rate 1-ρ.
    # We consider H_B (parity check for C_B) defines C_B^⊥ which has rate ρ.
    rate_Hᴬ = 1-ρ
    rate_Hᴮ = ρ
    Hᴬ = uniformly_random_code_checkmatrix(rate_Hᴬ, Δᴬ; rng=rng)
    Hᴮ = uniformly_random_code_checkmatrix(rate_Hᴮ, Δᴮ; rng=rng)
    @show Hᴬ
    @show Hᴮ
    if use_same_local_code
        Hᴮ = Hᴬ
    end
    # "The dual code of a code C is defined as C^⊥ = {x ∈ 𝔽₂ⁿ: ⟨x,y⟩=0 ∀ y ∈ C}" [gu2022efficient](@cite).
    Cᴬ = dual_code(Hᴬ)
    Cᴮ = dual_code(Hᴮ)
    @show Cᴬ
    @show Cᴮ
    C₀ = kronecker_product(Cᴬ, Cᴮ)
    C₁ = kronecker_product(Hᴬ, Hᴮ)
    @show size(C₀)
    @show size(C₁)
    @assert good_css(Hᴬ, Cᴬ)
    @assert good_css(Hᴮ, Cᴮ)
    @assert good_css(C₀, C₁)
    if bipartite
        checks_v₀ = tanner_code(𝒢₀□, edge₀_q_idx, edge₀_ab_idx, C₀)
        checks_v₁ = tanner_code(𝒢₁□, edge₁_q_idx, edge₁_ab_idx, C₁)
    else
        checks_v₀ = tanner_code_quadripartite(𝒢₀□, edge₀_q_idx, edge₀_ab_idx, C₀)
        checks_v₁ = tanner_code_quadripartite(𝒢₁□, edge₁_q_idx, edge₁_ab_idx, C₁)
    end
    if stabilizer_convention === :leverrier_zemor
        𝒞ˣ, 𝒞ᶻ = checks_v₁, checks_v₀
    else
        𝒞ˣ, 𝒞ᶻ = checks_v₀, checks_v₁
    end
    @show r1 = rank(𝒞ˣ)
    @show r2 = rank(𝒞ᶻ)
    @assert good_css(𝒞ˣ, 𝒞ᶻ)
    @assert _check_allrowscommute(𝒞ˣ, 𝒞ᶻ)
    return 𝒞ˣ, 𝒞ᶻ
end

"""
Generate a good Quantum Tanner code meeting minimum weight requirements.

### Arguments
- `ρ`: Rate parameter for local codes
- `group`, `A`, `B`: Expander graph components
- `minweight_x`: Minimum weight for X-stabilizers
- `minweight_z`: Minimum weight for Z-stabilizers
- `max_iterations`: Maximum attempts to find a good code
- `stabilizer_convention`: X/Z naming convention passed to
  [`random_quantum_Tanner_code`](@ref)
- `rng`: Random number generator used to construct candidate codes
"""
function gen_good_code(ρ::Real, group::Group, A::Vector{<:GroupElem}, B::Vector{<:GroupElem};
                       minweight_x=1, minweight_z=1, bipartite=true,
                       use_same_local_code=false, max_iterations=100,
                       stabilizer_convention::Symbol=:leverrier_zemor,
                       rng::AbstractRNG=GLOBAL_RNG)
    for i in 1:max_iterations
        𝒞ˣ, 𝒞ᶻ = random_quantum_Tanner_code(
            ρ,
            group,
            A,
            B;
            bipartite,
            use_same_local_code,
            stabilizer_convention,
            rng,
        )
        x_weight = minimum(unique(sum(𝒞ˣ, dims=1)))
        z_weight = minimum(unique(sum(𝒞ᶻ, dims=1)))
        if x_weight >= minweight_x && z_weight >= minweight_z
            @assert _check_allrowscommute(𝒞ˣ, 𝒞ᶻ)
            return 𝒞ˣ, 𝒞ᶻ
        end
    end
    throw(ArgumentError("Failed to generate code meeting weight constraints after $max_iterations iterations"))
end
