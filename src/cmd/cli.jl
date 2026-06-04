module Cli
    using ArgParse
    export parse_commandline, apply_blast_presets!, show_blast_presets, show_blast_params
    export BLAST_PRESETS, BLAST_CLI_DEFAULTS
    import ArgParse.parse_item
    using ArgParse: @add_arg_table!
    using Logging
    using Dates

    # Allowed values for the gene / chain selectors (shared across arg tables).
    const GENES = ["V", "D", "J"]
    const CHAINS = ["IGKV", "IGLV", "IGHV"]

    # Per-group argument tables live alongside this file in src/cmd/ (one builder per
    # command group). To add a new command: register it in add_command_groups! and add
    # its `@add_arg_table! s[<group>][<name>]` block to the matching group file.
    include("preprocess.jl")
    include("discover.jl")
    include("search.jl")
    include("analyze.jl")
    include("table.jl")
    include("fasta.jl")

    # CLI defaults for blast command - single source of truth (based on v0.0.66)
    const BLAST_CLI_DEFAULTS = Dict(
        "forward" => 20,
        "reverse" => 20,
        "minfullratio" => 0.1,
        "length" => 290,
        "maxdist" => 20,
        "minfullcount" => 5,
        "edge" => 0,
        "subjectcov" => 0.1,
        "minquality" => 0.75,
        "min-corecov" => 0.6,
        "args" => "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50",
        "work-dir" => ".immunediscover",
    )

    # V preset tuned for genomic IGHV novel recovery (see docs + tuning/README): relaxed trim/core coverage
    # and allelic ratio, longer affix extensions; BLAST capped at 5 alignments per query (see V_PRESET.md).
    const BLAST_PRESETS = Dict(
        "V" => Dict(
            "forward" => 20,
            "reverse" => 20,
            "minfullratio" => 0.035,
            "length" => 283,
            "maxdist" => 14,
            "minfullcount" => 5,
            "minquality" => 0.62,
            "min-corecov" => 0.50,
            "args" => "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50"
        ),
        "D" => Dict(
            "forward" => 40,
            "reverse" => 40,
            "minfullratio" => 0.2,
            "length" => 5,
            "maxdist" => 20,
            "minfullcount" => 10,
            "edge" => 10,
            "subjectcov" => 0.25,
            "minquality" => 0.5,
            "args" => "-task blastn -word_size 7 -xdrop_ungap 40 -xdrop_gap 40 -subject_besthit -num_alignments 10 -qcov_hsp_perc 5"
        ),
        "J" => Dict(
            "forward" => 12,
            "reverse" => 12,
            "minfullratio" => 0.1,
            "length" => 10,
            "maxdist" => 10,
            "minfullcount" => 10,
            "args" => "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 10"
        )
    )

    function show_blast_presets()
        println("BLAST Presets:")
        for (gene, settings) in BLAST_PRESETS
            println("\n$gene gene settings:")
            for (param, value) in settings
                println("  --$param = $value")
            end
        end
    end

    function get_blast_block(args)
        cmd = get(args, "%COMMAND%", "")
        cmd == "blast" && return args["blast"]
        cmd == "discover" && get(args["discover"], "%COMMAND%", "") == "blast" && return args["discover"]["blast"]
        return nothing
    end

    function show_blast_params(args)
        block = get_blast_block(args)
        block === nothing && return
        gene = block["gene"]
        if haskey(BLAST_PRESETS, gene)
            preset = BLAST_PRESETS[gene]
            println("BLAST Preset for $gene gene:")
            for (param, value) in preset
                println("  --$param = $(block[param])")
            end
        end
    end

    function apply_blast_presets!(parsed_args)
        block = get_blast_block(parsed_args)
        block === nothing && return parsed_args
        gene = block["gene"]
        if haskey(BLAST_PRESETS, gene)
            preset = BLAST_PRESETS[gene]

            # Only apply preset if the current value equals the CLI default
            # This means the user didn't explicitly override it
            for (key, preset_value) in preset
                if haskey(block, key) && haskey(BLAST_CLI_DEFAULTS, key)
                    current_value = block[key]
                    default_value = BLAST_CLI_DEFAULTS[key]
                    
                    # Apply preset only if current value equals the default
                    if current_value == default_value
                        @info "Applying $gene preset $key: $current_value → $preset_value"
                        block[key] = preset_value
                    else
                        @info "Keeping user override for $key: $current_value (not default $default_value)"
                    end
                else
                    # Apply preset for keys not in CLI defaults
                    if haskey(block, key)
                        @info "Applying $gene preset $key: $(get(block, key, nothing)) → $preset_value"
                        block[key] = preset_value
                    end
                end
            end
        end
        return parsed_args
    end

    # Lazy version detection — avoids running external `git` at const-initialization
    # time, which triggers a Julia 1.12 compiler inference bug.
    const version_cache = Ref{String}("")
    const hash_cache = Ref{String}("")

    function run_git_or_unknown(args::Cmd)
        Sys.which("git") === nothing && return "unknown"
        proc = run(pipeline(args, stderr=devnull), wait=false)
        wait(proc)
        success(proc) || return "unknown"
        return strip(read(args, String))
    end

    function software_version()
        if isempty(version_cache[])
            version_cache[] = run_git_or_unknown(`git -C $(@__DIR__) describe --tags --abbrev=0`)
        end
        return version_cache[]
    end

    function software_git_hash()
        if isempty(hash_cache[])
            hash_cache[] = run_git_or_unknown(`git -C $(@__DIR__) rev-parse HEAD`)
        end
        return hash_cache[]
    end

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
    function parse_commandline(args)
        s = ArgParseSettings("Tool for processing immune NGS data",
                            commands_are_required = true,
                            version = "$(software_version()) (git $(software_git_hash()))",
                            add_version = true,
                            usage = "usage: immunediscover <command> [-h|--help]",
                            epilog = "GKHLab, $(software_version()) (git $(software_git_hash()))")
        add_command_groups!(s)

        add_preprocess_args!(s)
        add_discover_args!(s)
        add_search_args!(s)
        add_analyze_args!(s)
        add_table_args!(s)
        add_fasta_args!(s)

        log_invocation(args)

        # CLI-boundary catch: ArgParse throws ArgParseError by design for invalid user input.
        # This is the standard pattern and acceptable at the CLI entry point.
        try
            return parse_args(args, s)
        catch e
            if e isa ArgParseError
                println(e)
                ArgParse.show_help(s)
            else
                rethrow()
            end
        end
    end
end
