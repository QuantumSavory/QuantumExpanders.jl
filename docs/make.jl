using Revise
push!(LOAD_PATH,"../src/")

using Documenter
using DocumenterCitations
using QuantumExpanders
using Oscar
using JuMP
using HiGHS
using QECCore
using QuantumClifford
using QuantumClifford.ECC

DocMeta.setdocmeta!(QuantumExpanders, :DocTestSetup, :(using QuantumClifford, QuantumExpanders); recursive=true)

ENV["LINES"] = 80
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
linkcheck = true
)

deploydocs(
    repo = "github.com/QuantumSavory/QuantumExpanders.jl.git"
)
