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

    "Build the ArgParse schema. Version strings are stamped at parse time."
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
        cli_help_pages()

    Root, each group, and each subcommand as `(name, argv)` in tree order.
    Cache files are `name.txt`; `argv` is the ArgParse `--help` invocation.
    """
    function cli_help_pages()
        pages = @NamedTuple{name::String, args::Vector{String}}[]
        push!(pages, (name="root", args=String["--help"]))
        seen = Set{String}()
        for cmd in COMMANDS
            group, sub = cli_path(cmd)
            if !(group in seen)
                push!(pages, (name=group, args=String[group, "--help"]))
                push!(seen, group)
            end
            push!(pages, (name="$group-$sub", args=String[group, sub, "--help"]))
        end
        return pages
    end

    function argparse_node(settings::ArgParseSettings, args)
        node = settings
        for a in args
            (a == "--help" || a == "-h") && continue
            node = node[a]
        end
        return node
    end

    """
        argparse_named_nodes(settings)

    Root, each group, and each subcommand as `(cache_name => node)` in tree order.
    Nested ArgParse tables copy `exit_after_help` at construction; callers that
    change the flag must apply it to every node.
    """
    function argparse_named_nodes(settings::ArgParseSettings)
        return [page.name => argparse_node(settings, page.args) for page in cli_help_pages()]
    end

    cli_has_flag(args, flag::AbstractString) = any(==(flag), args)
    cli_wants_version(args) = cli_has_flag(args, "--version") || cli_has_flag(args, "-V")
    cli_is_help_or_version(args) =
        cli_wants_version(args) || cli_has_flag(args, "--help") || cli_has_flag(args, "-h")

    stamp_cli_identity!(s, ::Val{false}) = stamp_cli_identity_label!(s, software_version())
    stamp_cli_identity!(s, ::Val{true}) =
        stamp_cli_identity_label!(s, software_version_label())
    function stamp_cli_identity_label!(s, label::AbstractString)
        s.version = label
        s.epilog = "GKHLab, $label"
        return s
    end

    function apply_parse_options!(settings::ArgParseSettings; exit_after_help::Bool, git_label::Bool=false)
        for (_, node) in argparse_named_nodes(settings)
            node.exit_after_help = exit_after_help
        end
        return stamp_cli_identity!(settings, Val(git_label))
    end

    const CLI_SETTINGS = Ref{ArgParseSettings}()

    function argparse_settings!(; exit_after_help::Bool, git_label::Bool=false)
        if !isassigned(CLI_SETTINGS)
            CLI_SETTINGS[] = build_argparse_settings()
        end
        return apply_parse_options!(CLI_SETTINGS[]; exit_after_help, git_label)
    end

    log_cli_invocation(args::Vector{String}) = log_cli_invocation(Val(cli_is_help_or_version(args)), args)
    log_cli_invocation(::Val{true}, _) = nothing
    log_cli_invocation(::Val{false}, args) = log_invocation(args)

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

    include("cli_help.jl")
end
