module Cli
    using ArgParse
    using ArgParse: @add_arg_table!
    using Logging
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
            with_logger(ConsoleLogger(io)) do
                @info "$(software_version()) $(Dates.now()) - Parsing command line arguments: $args"
            end
        end
        return nothing
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
        parse_commandline(args)

    Handle command line
    """
    function parse_commandline(args; exit_after_help::Bool=!isinteractive())
        s = ArgParseSettings("Tool for processing immune NGS data",
                            commands_are_required = true,
                            version = "$(software_version()) (git $(software_git_hash()))",
                            add_version = true,
                            usage = "usage: immunediscover <command> [-h|--help]",
                            epilog = "GKHLab, $(software_version()) (git $(software_git_hash()))",
                            exit_after_help = exit_after_help)
        add_command_groups!(s)

        add_preprocess_args!(s)
        add_discover_args!(s)
        add_search_args!(s)
        add_analyze_args!(s)
        add_table_args!(s)
        add_fasta_args!(s)

        log_invocation(args)

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
