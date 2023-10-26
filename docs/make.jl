using immunediscover
using Documenter

DocMeta.setdocmeta!(immunediscover, :DocTestSetup, :(using immunediscover); recursive=true)

makedocs(;
    modules=[immunediscover],
    authors="Mateusz Kaduk <mateusz.kaduk@gmail.com> and contributors",
    sitename="immunediscover.jl",
    checkdocs = :exports,
    repo="https://gitlab.com/mateusz-kaduk/immunediscover.jl/blob/{commit}{path}#{line}",
    format=Documenter.HTML(repolink="https://gitlab.com/mateusz-kaduk/immunediscover.jl/"),
    warnonly = true,
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(
    repo = "git@gitlab.com:mateusz-kaduk/immunediscover.jl.git"
)
