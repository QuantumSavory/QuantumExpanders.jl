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

ENV["LINES"] = 80
ENV["COLUMNS"] = 80

bib = CitationBibliography(
    joinpath(@__DIR__, "src/references.bib"),
    style = :authoryear,
)

makedocs(
    plugins = [bib],
    doctest = false,
    clean = true,
    warnonly = :missing_docs,
    sitename = "QuantumExpanders.jl",
    format = Documenter.HTML(),
    authors = "Feroz Ahmed Mian, Stefan Krastanov, Vaishnavi Addala, QuantumSavory community members",
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Quantum Tanner codes" => [
            "LRCC construction" => "quantum_tanner.md",
            "Lifted left-right actions" => "quantum_tanner_left_right_actions.md",
            "Reproducing the manuscript" => "paper_instances.md",
        ],
        "Ramanujan graphs" => [
            "Lubotzky–Phillips–Sarnak" => "lps.md",
            "Morgenstern" => "morgenstern.md",
        ],
        "API reference" => "API.md",
    ],
    linkcheck = true,
)

deploydocs(
    repo = "github.com/QuantumSavory/QuantumExpanders.jl.git",
    devbranch = "main",
    deploy_config = Documenter.GitHubActions(),
    push_preview = true,
)
