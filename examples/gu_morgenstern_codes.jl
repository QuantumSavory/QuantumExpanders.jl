# This file contains a "scratch pad" of various tests and demos related to generating
# Tanner codes based as described in Gu et al, but using the Morgenstern graphs.
#
# We do not have |A|=|B| so we are not exactly following the original prescription.
# In particular, we have taken a lot of liberties with the generation of the lower tensor codes.
#
# Not particularly careful with what is getting imported where, would clean up when
# turning it into a package.

using QuantumExpanders

using CairoMakie # Plotting libraries that are extremely slow to load due to
                 # lattency limitations in the Julia compiler. Otherwise great and fast.

using DataFrames
using AlgebraOfGraphics
using CSV

df = DataFrame(rank_a=Int[], rank_b=Int[], min_weight_a=Int[], min_weight_b=Int[], size_z=Tuple{Int, Int}[], size_x=Tuple{Int, Int}[], rank_z=Int[], rank_x=Int[], rank_overall=Int[], size_overall=Tuple{Int, Int}[], bipartite=Bool[], generators=String[])
for _ in 1:30
    for bipartite in [true, false]
        for generators in ["A,B", "A", "B"]    
            for ra in 1:3 # TODO - with quadripartite, the ranks and min weights need to be adjusted
                for rb in 1:2
                    for min_wt_a in 1:2
                        for min_wt_b in 1:2
                            res=gen_good_code(ra,rb, min_wt_a, min_wt_b, bipartite, generators)
                            if res === nothing
                                continue
                            end
                            𝒞ˣ,𝒞ᶻ=res;
                            size_z, size_x = size(𝒞ᶻ), size(𝒞ˣ)
                            rank_z, rank_x = rank(𝒞ᶻ), rank(𝒞ˣ)
                            (sum(𝒞ᶻ, dims=1) |> unique), (sum(𝒞ˣ, dims=1) |> unique)
                            code = css(𝒞ˣ,𝒞ᶻ)
                            size_overall = size(code)
                            rank_overall = rank(MixedDestabilizer(code))

                            push!(df, (ra, rb, min_wt_a, min_wt_b, size_z, size_x, rank_z, rank_x, rank_overall, size_overall, bipartite, generators))
                        end
                    end
                end
            end
        end
    end
end

CSV.write("codes_quad_parity_checks.csv", df)
df_test = AlgebraOfGraphics.data(df) * mapping(:bipartite, :rank_overall) * mapping(color=:generators) #* mapping(marker=:bipartite) 
CairoMakie.save("codes_quad.png", draw(df_test))

##

using QuantumClifford
#using QuantumCliffordPlots

stab = css(𝒞ˣ,𝒞ᶻ)
@assert good_css(stab)
for row in stab
    @assert all(==(0), QuantumClifford.comm(stab[1],stab))
end

QuantumClifford.stab_looks_good(stab) # internal function used for sanity checks
#stabilizerplot(stab)
