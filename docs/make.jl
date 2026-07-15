using Revise
push!(LOAD_PATH,"../src/")

using Documenter
using DocumenterCitations, DocumenterMermaid
using QuantumExpanders
using Oscar
using JuMP
using HiGHS
using QECCore
using QuantumClifford
using QuantumClifford.ECC

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
pages = [
    "QuantumExpanders.jl" => "index.md",
    "Morgenstern Ramanujan Graphs" => "morgenstern.md",
    "Lubotzky–Phillips–Sarnak Ramanujan Graphs" => "lps.md",
    "API" => "API.md",
],
linkcheck = true
)

deploydocs(
    repo = "github.com/QuantumSavory/QuantumExpanders.jl.git",
    devbranch = "main",
    deploy_config = Documenter.GitHubActions(),
    push_preview = true
)
