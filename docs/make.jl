using Dysts
using Documenter

DocMeta.setdocmeta!(Dysts, :DocTestSetup, :(using Dysts); recursive=true)

makedocs(;
    modules=[Dysts],
    authors="Nathanael Bosch <nathanael.bosch@uni-tuebingen.de> and contributors",
    repo="https://github.com/nathanaelbosch/Dysts.jl/blob/{commit}{path}#{line}",
    sitename="Dysts.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nathanaelbosch.github.io/Dysts.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/nathanaelbosch/Dysts.jl",
    devbranch="main",
)
