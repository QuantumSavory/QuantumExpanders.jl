# QuantumExpanders.jl

- `main_gu_morgenstern_codes.jl` - Create a CSS code out of two Tanner codes 𝒞ᶻ and 𝒞ˣ, each constructed out of two related graphs 𝒢₀□, 𝒢₁□, the graphs build out of a Cayley complex, which itself was based on a group G=SL₂qⁱ and two generator sets A and B. For consistency of indexing of "qubits" ≈ "graph edges" ≈ "squares" and "local bits" ≈ "generator pairs" we have the `edge_*_idx` maps.

- `main_morgenstern_graphs.jl` - Generate Cayley graphs based on G=SL₂qⁱ and two generator sets as given by Morgenstern and follow-up works. Study the expander properties of these graphs.