using ChanVese
using Documenter

DocMeta.setdocmeta!(ChanVese, :DocTestSetup, :(using ChanVese); recursive=true)

makedocs(;
    modules=[ChanVese],
    authors="Dale-Black <djblack@uci.edu> and contributors",
    repo="https://github.com/Dale-Black/ChanVese.jl/blob/{commit}{path}#{line}",
    sitename="ChanVese.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Dale-Black.github.io/ChanVese.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Dale-Black/ChanVese.jl",
)
