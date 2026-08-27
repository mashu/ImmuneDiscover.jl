# ArgParse help text for the command tree, plus a gitignored on-disk cache so
# `scripts/run.sh` can print --help without starting Julia.

const CLI_HELP_DIR = joinpath(dirname(dirname(@__DIR__)), "build", "help")

help_invocation_args() = [page.args for page in cli_help_pages()]

"ArgParse `--help` text for `args` (includes generated --help/--version lines)."
function help_text(settings::ArgParseSettings, args::Vector{String})
    mktemp() do path, io
        redirect_stdout(io) do
            parse_args(args, settings)
        end
        flush(io)
        read(path, String)
    end
end

"""
    write_cli_help_pages!(dir=CLI_HELP_DIR)

Write one `.txt` file per argparse node (`root`, `<group>`, `<group>-<sub>`).
Uses a throwaway schema so the process-wide parse cache is not mutated.
"""
function write_cli_help_pages!(dir::AbstractString=CLI_HELP_DIR)
    settings = apply_parse_options!(build_argparse_settings(); exit_after_help=false)
    mkpath(dir)
    for page in cli_help_pages()
        write(joinpath(dir, page.name * ".txt"), help_text(settings, page.args))
    end
    return dir
end

function help_cache_sources()
    cmd_dir = @__DIR__
    sources = String[joinpath(cmd_dir, name) for name in readdir(cmd_dir) if endswith(name, ".jl")]
    push!(sources, project_toml)
    return sources
end

function help_cache_stale(dir::AbstractString=CLI_HELP_DIR)
    stamp = joinpath(dir, "root.txt")
    isfile(stamp) || return true
    t = mtime(stamp)
    return any(path -> isfile(path) && mtime(path) > t, help_cache_sources())
end

function ensure_help_pages!(dir::AbstractString=CLI_HELP_DIR)
    help_cache_stale(dir) && write_cli_help_pages!(dir)
    return dir
end

"Parse --help/--version on a throwaway schema so the package image does not cache ArgParseSettings."
function trace_cli_help_parse!()
    s = apply_parse_options!(build_argparse_settings(); exit_after_help=false)
    parse_args(String["--help"], s)
    parse_args(String["--version"], s)
    return nothing
end
