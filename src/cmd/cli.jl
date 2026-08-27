module Cli
    using ArgParse
    using ArgParse: @add_arg_table!
    using Dates
    using ..Option: Absent, Present, absent, optional

    export parse_commandline, apply_blast_presets!, show_blast_presets
    export BLAST_PRESETS, BLAST_DEFAULTS, blast_default, BLAST_PARAM_GROUPS

    import ArgParse.parse_item

    const GENES = ["V", "D", "J"]
    const CHAINS = ["IGKV", "IGLV", "IGHV"]

    include("cli_version.jl")
    include("cli_blast_presets.jl")
    include("cli_command.jl")

    # Per-group argument tables live alongside this file in src/cmd/ (one builder per
    # command group). To add a new command: register it in add_command_groups! and add
    # its `@add_arg_table! s[<group>][<name>]` block to the matching group file.
    include("preprocess.jl")
    include("discover.jl")
    include("search.jl")
    include("analyze.jl")
    include("table.jl")
    include("fasta.jl")

    """
        always_gz(file_path)

    Return path that ends with .gz, appending the extension if needed.
    """
    function always_gz(file_path)
        endswith(file_path, ".gz") ? file_path : file_path * ".gz"
    end

    "Append the invocation to immunediscover.log; never aborts the run."
    function log_invocation(args)
        logpath = "immunediscover.log"
        isdir(dirname(abspath(logpath))) || return nothing
        open(logpath, "a") do io
            println(io, software_version(), " ", Dates.now(),
                    " - Parsing command line arguments: ", args)
        end
        return nothing
    end

    "Build the ArgParse schema once. Version strings are stamped at parse time (no git at const-init)."
    function build_argparse_settings()
        s = ArgParseSettings("Tool for processing immune NGS data",
                            prog = "immunediscover",
                            commands_are_required = true,
                            version = "",
                            add_version = true,
                            usage = "usage: immunediscover <command> [-h|--help]",
                            epilog = "",
                            exit_after_help = true)
        add_command_groups!(s)
        add_preprocess_args!(s)
        add_discover_args!(s)
        add_search_args!(s)
        add_analyze_args!(s)
        add_table_args!(s)
        add_fasta_args!(s)
        return s
    end

    const CLI_SETTINGS = Ref{Any}(nothing)

    function foreach_command_settings!(f, s)
        f(s)
        seen = Set{String}()
        for cmd in COMMANDS
            group, sub = cli_path(cmd)
            if !(group in seen)
                f(s[group])
                push!(seen, group)
            end
            f(s[group][sub])
        end
        return nothing
    end

    function argparse_settings!(; exit_after_help::Bool, git_label::Bool=false)
        s = CLI_SETTINGS[]
        if s === nothing
            s = build_argparse_settings()
            CLI_SETTINGS[] = s
        end
        foreach_command_settings!(node -> (node.exit_after_help = exit_after_help), s)
        return stamp_cli_identity!(s, Val(git_label))
    end

    const CLI_HELP_DIR = joinpath(dirname(dirname(@__DIR__)), "build", "help")

    "Help-page stem: `root`, `search`, or `search-exact`."
    function help_page_key(args::Vector{String})
        parts = String[]
        for a in args
            (a == "--help" || a == "-h" || a == "--version" || a == "-V") && continue
            push!(parts, a)
        end
        return isempty(parts) ? "root" : join(parts, "-")
    end

    function help_page_args()
        pages = Vector{Vector{String}}()
        push!(pages, String["--help"])
        seen = Set{String}()
        for cmd in COMMANDS
            group, sub = cli_path(cmd)
            if !(group in seen)
                push!(pages, String[group, "--help"])
                push!(seen, group)
            end
            push!(pages, String[group, sub, "--help"])
        end
        return pages
    end

    function help_settings_for(args::Vector{String})
        s = argparse_settings!(; exit_after_help=false, git_label=false)
        node = s
        for a in args
            (a == "--help" || a == "-h" || a == "--version" || a == "-V") && continue
            node = node[a]
        end
        return node
    end

    function capture_cli_help(args::Vector{String})
        s = help_settings_for(args)
        sprint() do io
            ArgParse.show_help(io, s; exit_when_done=false)
        end
    end

    """
        write_cli_help_pages!(dir=CLI_HELP_DIR)

    Write ArgParse `--help` text for every command. `scripts/run.sh` regenerates this
    cache when `src/cmd/` is newer; the files are gitignored under `build/help/`.
    """
    function write_cli_help_pages!(dir::AbstractString=CLI_HELP_DIR)
        mkpath(dir)
        for args in help_page_args()
            write(joinpath(dir, help_page_key(args) * ".txt"), capture_cli_help(args))
        end
        return dir
    end

    "Drop the cached ArgParse schema so a precompile image does not serialize it."
    function reset_cli_settings!()
        CLI_SETTINGS[] = nothing
        return nothing
    end

    cli_has_flag(args, flag::AbstractString) = any(==(flag), args)
    cli_wants_version(args) = cli_has_flag(args, "--version") || cli_has_flag(args, "-V")
    cli_is_help_or_version(args) =
        cli_wants_version(args) || cli_has_flag(args, "--help") || cli_has_flag(args, "-h")

    stamp_cli_identity!(s, ::Val{false}) = stamp_cli_identity_label!(s, software_version())
    stamp_cli_identity!(s, ::Val{true}) =
        stamp_cli_identity_label!(s, "$(software_version()) (git $(software_git_hash()))")
    function stamp_cli_identity_label!(s, label::AbstractString)
        s.version = label
        s.epilog = "GKHLab, $label"
        return s
    end

    log_cli_invocation(args::Vector{String}) = log_cli_invocation(Val(cli_is_help_or_version(args)), args)
    log_cli_invocation(::Val{true}, _) = nothing
    log_cli_invocation(::Val{false}, args) = log_invocation(args)

    "Register the top-level command groups on the settings object."
    function add_command_groups!(s)
        @add_arg_table! s begin
            "discover"
                help = "De novo allele discovery (blast, hsmm)"
                action = :command
            "search"
                help = "Search against known references (exact, heptamer, bwa)"
                action = :command
            "analyze"
                help = "Downstream analysis (cooccurrence, haplotype)"
                action = :command
            "preprocess"
                help = "Data preparation (demultiplex)"
                action = :command
            "table"
                help = "Table utilities (outerjoin, leftjoin, transform, aggregate, unique, sort, filter, select, fasta, collect, exclude)"
                action = :command
            "fasta"
                help = "FASTA utilities (merge, diff, hash)"
                action = :command
        end
        return s
    end

    """
        parse_commandline(args)

    Handle command line. The ArgParse schema is built once per process and reused.
    """
    parse_commandline(args::AbstractVector{<:AbstractString}; kwargs...) =
        parse_commandline(String[a for a in args]; kwargs...)

    function parse_commandline(args::Vector{String}; exit_after_help::Bool=!isinteractive())
        s = argparse_settings!(; exit_after_help=exit_after_help, git_label=cli_wants_version(args))
        log_cli_invocation(args)

        # CLI-boundary catch: ArgParse throws ArgParseError by design for invalid user input.
        try
            return parse_args(args, s)
        catch e
            handle_parse_error(e, s)
        end
    end

    handle_parse_error(e::ArgParseError, s) = (println(e); ArgParse.show_help(s); nothing)
    handle_parse_error(e, _) = rethrow(e)
end
