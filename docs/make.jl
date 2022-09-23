using immunediscover
using Documenter

DocMeta.setdocmeta!(immunediscover, :DocTestSetup, :(using immunediscover); recursive=true)

makedocs(;
    modules=[immunediscover],
    authors="Mateusz Kaduk <mateusz.kaduk@gmail.com> and contributors",
    repo="https://gitlab.com/mateusz-kaduk/immunediscover.jl/blob/{commit}{path}#{line}",
    sitename="immunediscover.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mateusz-kaduk.gitlab.io/immunediscover.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
