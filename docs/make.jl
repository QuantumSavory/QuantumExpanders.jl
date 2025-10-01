push!(LOAD_PATH,"../src/")

using Documenter
using DocumenterCitations
using QuantumExpanders

#bib = CitationBibliography(joinpath(@__DIR__,"src/references.bib"))

ENV["LINES"] = 80    # for forcing `displaysize(io)` to be big enough
ENV["COLUMNS"] = 80

bib = CitationBibliography(joinpath(@__DIR__,"src/references.bib"),style=:authoryear)

makedocs(
plugins = [bib],
doctest = false,
clean = true,
warnonly = :missing_docs,
sitename = "QuantumExpanders.jl",
format = Documenter.HTML(),
authors = "Feroz Ahmed Mian, Stefan Krastanov, Vaishnavi Addala, QuantumSavory community members",
pages = ["API" => "API.md"],
linkcheck = true,
)

deploydocs(
    repo = "github.com/QuantumSavory/QuantumExpanders.jl.git"
)
