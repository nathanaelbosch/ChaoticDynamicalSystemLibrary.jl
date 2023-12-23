using ChaoticDynamicalSystemLibrary
using Documenter

DocMeta.setdocmeta!(ChaoticDynamicalSystemLibrary,
    :DocTestSetup,
    :(using ChaoticDynamicalSystemLibrary);
    recursive = true)

makedocs(;
    modules = [ChaoticDynamicalSystemLibrary],
    authors = "Nathanael Bosch <nathanael.bosch@uni-tuebingen.de> and contributors",
    repo = "https://github.com/nathanaelbosch/ChaoticDynamicalSystemLibrary.jl/blob/{commit}{path}#{line}",
    sitename = "ChaoticDynamicalSystemLibrary.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://nathanaelbosch.github.io/ChaoticDynamicalSystemLibrary.jl",
        edit_link = "main",
        assets = String[],
        size_threshold_ignore = ["index.md"],),
    pages = [
        "Home" => "index.md",
    ],)

deploydocs(;
    repo = "github.com/nathanaelbosch/ChaoticDynamicalSystemLibrary.jl",
    devbranch = "main",)
