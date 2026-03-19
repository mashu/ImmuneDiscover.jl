#!/usr/bin/env julia
# Build a standalone immunediscover binary for the current platform.
# Uses PackageCompiler + optional precompile workload for fast startup.
#
# Usage (from repo root):
#   julia --project=build scripts/build_binary.jl [output_dir]
#
# Or from build/ with PackageCompiler already in the project:
#   julia -e 'using Pkg; Pkg.activate("build"); include("../scripts/build_binary.jl")'
#
# Output: output_dir/bin/immunediscover (or build/immunediscover_app/bin/immunediscover if no arg)
# The binary is relocatable; copy the whole output_dir to another machine with the same OS/arch.

using Pkg
using PackageCompiler

const SCRIPT_DIR = @__DIR__
const REPO_ROOT = abspath(joinpath(SCRIPT_DIR, ".."))
const DEFAULT_APP_DIR = joinpath(REPO_ROOT, "build", "immunediscover_app")
const PRECOMPILE_SCRIPT = joinpath(SCRIPT_DIR, "precompile_workload.jl")

function main()
    app_dir = length(ARGS) >= 1 ? abspath(ARGS[1]) : DEFAULT_APP_DIR
    if !isdir(REPO_ROOT) || !isfile(joinpath(REPO_ROOT, "Project.toml"))
        error("Repo root not found or not a Julia project: $REPO_ROOT")
    end
    precompile_file = isfile(PRECOMPILE_SCRIPT) ? PRECOMPILE_SCRIPT : String[]
    @info "Building immunediscover app" package_dir=REPO_ROOT app_dir=app_dir precompile=precompile_file
    create_app(
        REPO_ROOT,
        app_dir;
        executables = ["immunediscover" => "julia_main"],
        precompile_execution_file = precompile_file,
        incremental = false,
        force = true,
    )
    exe = joinpath(app_dir, "bin", "immunediscover")
    @info "Build complete" executable=exe
    println("Run: ", exe, " --help")
    return nothing
end

main()
