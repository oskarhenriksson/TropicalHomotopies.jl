using Documenter
using TropicalHomotopies
using Oscar

DocMeta.setdocmeta!(
    TropicalHomotopies,
    :DocTestSetup,
    quote
        using VerticalRootCounts
        using Oscar
    end;
    recursive = true,
)

makedocs(
    sitename = "TropicalHomotopies.jl",
    modules = [TropicalHomotopies],
)

deploydocs(
    repo = "github.com/oskarhenriksson/TropicalHomotopies.jl.git",
)