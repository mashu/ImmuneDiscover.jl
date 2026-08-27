#!/usr/bin/env julia
# Custom sysimage: native code for immunediscover and every dependency (CSV, DataFrames, …).
# Local use only (CPU-native). Rebuild after changing source or Manifest.toml.
#
# Usage (from repo root):
#   ./scripts/build_sysimage.sh
#   julia --project=build scripts/build_sysimage.jl [output.so]

using PackageCompiler

const SCRIPT_DIR = @__DIR__
const REPO_ROOT = abspath(joinpath(SCRIPT_DIR, ".."))
const DEFAULT_SYSIMAGE = joinpath(REPO_ROOT, "build", "immunediscover.so")
const PRECOMPILE_SCRIPT = joinpath(SCRIPT_DIR, "precompile_workload.jl")

function main()
    sysimage_path = length(ARGS) >= 1 ? abspath(ARGS[1]) : DEFAULT_SYSIMAGE
    if !isfile(joinpath(REPO_ROOT, "Project.toml"))
        error("Repo root not found or not a Julia project: $REPO_ROOT")
    end
    mkpath(dirname(sysimage_path))
    precompile_file = isfile(PRECOMPILE_SCRIPT) ? PRECOMPILE_SCRIPT : String[]
    @info "Building immunediscover sysimage" project=REPO_ROOT sysimage_path=sysimage_path
    create_sysimage(
        [:immunediscover];
        sysimage_path=sysimage_path,
        project=REPO_ROOT,
        precompile_execution_file=precompile_file,
        incremental=true,
        include_transitive_dependencies=true,
    )
    @info "Sysimage ready" sysimage_path=sysimage_path
    println("scripts/run.sh will use this automatically.")
    return nothing
end

main()
