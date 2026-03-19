using immunediscover
using Documenter

DocMeta.setdocmeta!(immunediscover, :DocTestSetup, :(using immunediscover); recursive=true)

# Explicit branch/edit_link so Documenter does not need to call `git remote`
const REPO = "https://gitlab.com/mateusz-kaduk/immunediscover.jl"
const EDIT_BRANCH = "main"

makedocs(;
    modules=[immunediscover],
    authors="Mateusz Kaduk <mateusz.kaduk@gmail.com> and contributors",
    sitename="immunediscover.jl",
    checkdocs=:exports,
    repo="$REPO/blob/{commit}{path}#{line}",
    format=Documenter.HTML(
        repolink="$REPO/",
        edit_link=EDIT_BRANCH,
    ),
    warnonly=true,
    pages=[
        "Home" => "index.md",
        "Commands" => "commands.md",
        "User Guide" => [
            "Workflows" => "workflows.md",
            "Parameters" => "parameters.md",
        ],
        "Reference" => [
            "API" => "api.md",
            "Building a binary" => "building.md",
            "Output Columns" => "columns.md",
            "Troubleshooting" => "troubleshooting.md",
        ],
    ],
)

deploydocs(;
    repo="git@gitlab.com:mateusz-kaduk/immunediscover.jl.git",
    devbranch=EDIT_BRANCH,
)
