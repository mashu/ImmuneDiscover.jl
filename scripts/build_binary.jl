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
        filter_stdlibs = true,
        include_transitive_dependencies = false,
        force = true,
    )
    exe = joinpath(app_dir, "bin", "immunediscover")
    install_fast_cli_wrapper!(app_dir)
    @info "Build complete" executable=exe
    println("Run: ", exe, " --help")
    return nothing
end

function read_project_version(toml_path)
    for line in eachline(toml_path)
        m = match(r"^version\s*=\s*\"(.+)\"\s*$", line)
        m !== nothing && return String(m.captures[1])
    end
    return "unknown"
end

function read_git_hash(repo)
    Sys.which("git") === nothing && return "unknown"
    proc = run(pipeline(`git -C $repo rev-parse HEAD`, stderr=devnull); wait=false)
    wait(proc)
    success(proc) || return "unknown"
    return strip(read(`git -C $repo rev-parse HEAD`, String))
end

"Wrap the PackageCompiler trampoline so --help/--version never load the sysimage."
function install_fast_cli_wrapper!(app_dir)
    Sys.iswindows() && return nothing
    help_src = joinpath(REPO_ROOT, "build", "help")
    if !isfile(joinpath(help_src, "root.txt"))
        julia = joinpath(Sys.BINDIR, Base.julia_exename())
        mkpath(help_src)
        expr = "using immunediscover; immunediscover.Cli.write_cli_help_pages!(ARGS[1])"
        run(`$julia --startup-file=no --quiet --project=$REPO_ROOT -e $expr -- $help_src`)
    end
    isdir(help_src) || return nothing
    trampoline = joinpath(app_dir, "bin", "immunediscover")
    isfile(trampoline) || return nothing

    share = joinpath(app_dir, "share", "immunediscover")
    help_dst = joinpath(share, "help")
    mkpath(help_dst)
    for file in readdir(help_src; join=true)
        endswith(file, ".txt") || continue
        cp(file, joinpath(help_dst, basename(file)); force=true)
    end
    cp(joinpath(SCRIPT_DIR, "cli_fastpath.sh"), joinpath(share, "cli_fastpath.sh"); force=true)
    write(joinpath(share, "VERSION"),
          "$(read_project_version(joinpath(REPO_ROOT, "Project.toml"))) (git $(read_git_hash(REPO_ROOT)))\n")

    realbin = joinpath(app_dir, "bin", "immunediscover.bin")
    mv(trampoline, realbin; force=true)
    write(trampoline, """
#!/bin/sh
set -eu
ROOT="\$(CDPATH= cd -- "\$(dirname "\$0")/.." && pwd)"
HELP_DIR="\$ROOT/share/immunediscover/help"
VERSION_FILE="\$ROOT/share/immunediscover/VERSION"
print_cli_version() { cat "\$VERSION_FILE"; }
. "\$ROOT/share/immunediscover/cli_fastpath.sh"
if try_fast_help_or_version "\$@"; then
    exit 0
fi
exec "\$ROOT/bin/immunediscover.bin" "\$@"
""")
    chmod(trampoline, 0o755)
    return nothing
end

main()
