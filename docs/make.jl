using PlugFlowReactor
using Documenter

DocMeta.setdocmeta!(PlugFlowReactor, :DocTestSetup, :(using PlugFlowReactor); recursive=true)

makedocs(;
    modules=[PlugFlowReactor],
    authors="Vinod Janardhanan",
    repo="https://github.com/vinodjanardhanan/PlugFlowReactor.jl/blob/{commit}{path}#{line}",
    sitename="PlugFlowReactor.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://vinodjanardhanan.github.io/PlugFlowReactor.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/vinodjanardhanan/PlugFlowReactor.jl",
    devbranch="main",
)
