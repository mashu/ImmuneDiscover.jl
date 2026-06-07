module Cli
    using ArgParse
    export parse_commandline, apply_blast_presets!, show_blast_presets
    export BLAST_PRESETS, BLAST_DEFAULTS, blast_default, BLAST_PARAM_GROUPS

    # Logical grouping for displaying `discover blast` parameters (order matters; any key not
    # listed falls under "other" so nothing is hidden).
    const BLAST_PARAM_GROUPS = [
        "Inputs / outputs"          => ["input", "fasta", "pseudo", "output", "full-output", "work-dir"],
        "Gene preset"               => ["gene", "show-presets"],
        "Extension & trimming"      => ["forward", "reverse", "minquality", "min-corecov"],
        "BLAST search"              => ["args", "maxdist", "edge", "subjectcov", "min-read-length"],
        "Cluster & output filters"  => ["minfullcount", "minfullratio", "min-reads-total", "length", "isin", "keep-failed"],
        "Quality-metric filters"    => ["min-recurrence", "max-homopolymer"],
        "Run control"               => ["overwrite", "verbose"],
    ]
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

    # ─── `discover blast` tunable parameters: ONE source of truth ────────────────────────
    # `BLAST_DEFAULTS` holds every tunable's global default; `discover.jl` reads them via
    # `blast_default(key)` for its ArgParse `default=`, so the arg table and this table can
    # never drift. Each gene preset lists ONLY the keys it changes from the default — a value
    # equal to the default is omitted (that is why the log no longer shows "20 → 20" no-ops).
    const BLAST_DEFAULTS = Dict{String,Any}(
        "forward"         => 20,
        "reverse"         => 20,
        "minquality"      => 0.75,
        "min-corecov"     => 0.6,
        "args"            => "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 50",
        "maxdist"         => 20,
        "edge"            => 0,
        "subjectcov"      => 0.1,
        "min-read-length" => 0,
        "minfullcount"    => 5,
        "minfullratio"    => 0.1,
        "min-reads-total" => 0,
        "length"          => 290,
        "min-recurrence"  => 0,
        "max-homopolymer" => 0,
        "work-dir"        => ".immunediscover",
    )

    "Global default for a `discover blast` parameter (single source for the ArgParse table)."
    blast_default(key::AbstractString) = BLAST_DEFAULTS[key]

    # Gene presets tuned on KI IGH self-tests (recovery of known-novel alleles; see selftest).
    #   V: `minfullratio 0.08` is the key false-positive cut — a germline allele is a major
    #      allele (peak per-donor allelic ratio ≥ 0.085) in at least one carrier, while PCR /
    #      sequencing artifacts never are. 0.08 keeps every truth-novel allele (recall 1.0)
    #      while removing ~60% of false novel calls vs the old 0.035.
    const BLAST_PRESETS = Dict(
        "V" => Dict{String,Any}(
            "minfullratio" => 0.08,
            "length"       => 283,
            "maxdist"      => 14,
            "minquality"   => 0.62,
            "min-corecov"  => 0.50,
        ),
        "D" => Dict{String,Any}(
            "forward"      => 40,
            "reverse"      => 40,
            "minfullratio" => 0.2,
            "length"       => 5,
            "minfullcount" => 10,
            "edge"         => 10,
            "subjectcov"   => 0.25,
            "minquality"   => 0.5,
            "args"         => "-task blastn -word_size 7 -xdrop_ungap 40 -xdrop_gap 40 -subject_besthit -num_alignments 10 -qcov_hsp_perc 5",
        ),
        "J" => Dict{String,Any}(
            "forward"      => 12,
            "reverse"      => 12,
            "length"       => 10,
            "maxdist"      => 10,
            "minfullcount" => 10,
            "args"         => "-task megablast -subject_besthit -num_alignments 5 -qcov_hsp_perc 10",
        ),
    )

    "One parameter the gene preset wants to set; `applied` is false when a user override is kept."
    struct PresetChange
        key::String
        current::Any
        preset::Any
        applied::Bool
    end

    """
        preset_changes(block, gene) -> Vector{PresetChange}

    Resolve what the `gene` preset would do to the parsed `block`. A preset value is applied
    only when the current value still equals the global default (i.e. the user did not pass
    it explicitly); otherwise the user override is kept. Sorted by key for stable logging.
    """
    function preset_changes(block::AbstractDict, gene::AbstractString)
        preset = BLAST_PRESETS[gene]
        changes = PresetChange[]
        for key in sort!(collect(keys(preset)))
            haskey(block, key) || continue
            current = block[key]
            applied = current == blast_default(key)
            push!(changes, PresetChange(key, current, preset[key], applied))
        end
        return changes
    end

    "Single aligned block describing the preset outcome (replaces the per-key @info spam)."
    function log_preset_changes(gene::AbstractString, changes::AbstractVector{PresetChange})
        printstyled("━━ $gene gene preset ", "━"^40, "\n"; color=:cyan, bold=true)
        applied = filter(c -> c.applied, changes)
        kept    = filter(c -> !c.applied, changes)
        kew = isempty(changes) ? 0 : maximum(length(c.key) for c in changes)
        if !isempty(applied)
            println("  applied (param was default → preset value)")
            for c in applied
                println("    ", rpad(c.key, kew), "  ", c.current, " → ", c.preset)
            end
        end
        if !isempty(kept)
            println("  kept (explicit user override; preset skipped)")
            for c in kept
                println("    ", rpad(c.key, kew), "  ", c.current, "  (preset ", c.preset, ")")
            end
        end
        return nothing
    end

    "Print every gene preset as the delta from the defaults (wired to `--show-presets`)."
    function show_blast_presets()
        printstyled("BLAST gene presets (only keys that differ from the global default)\n";
                    color=:cyan, bold=true)
        for gene in sort!(collect(keys(BLAST_PRESETS)))
            println("\n$gene gene:")
            preset = BLAST_PRESETS[gene]
            for key in sort!(collect(keys(preset)))
                println("  --$key  $(blast_default(key)) → $(preset[key])")
            end
        end
        return nothing
    end

    function get_blast_block(args)
        cmd = get(args, "%COMMAND%", "")
        cmd == "blast" && return args["blast"]
        cmd == "discover" && get(args["discover"], "%COMMAND%", "") == "blast" && return args["discover"]["blast"]
        return nothing
    end

    """
        apply_blast_presets!(parsed_args) -> parsed_args

    Apply the selected gene preset in place: each preset key takes the preset value unless the
    user passed it explicitly (detected as "current value ≠ default"). Logs one tidy block.
    """
    function apply_blast_presets!(parsed_args)
        block = get_blast_block(parsed_args)
        block === nothing && return parsed_args
        gene = block["gene"]
        haskey(BLAST_PRESETS, gene) || return parsed_args
        changes = preset_changes(block, gene)
        log_preset_changes(gene, changes)
        for c in changes
            c.applied && (block[c.key] = c.preset)
        end
        return parsed_args
    end

    # Lazy version detection — avoids running external `git` at const-initialization
    # time, which triggers a Julia 1.12 compiler inference bug.
    const version_cache = Ref{String}("")
    const hash_cache = Ref{String}("")
    const project_toml = abspath(joinpath(@__DIR__, "..", "..", "Project.toml"))

    function run_git_or_unknown(args::Cmd)
        Sys.which("git") === nothing && return "unknown"
        proc = run(pipeline(args, stderr=devnull), wait=false)
        wait(proc)
        success(proc) || return "unknown"
        return strip(read(args, String))
    end

    "Read `version` from the package Project.toml; empty when unavailable."
    function read_project_version()
        isfile(project_toml) || return ""
        for line in eachline(project_toml)
            m = match(r"""^version\s*=\s*"(.+)"\s*$""", line)
            m !== nothing && return m.captures[1]
        end
        return ""
    end

    function software_version()
        if isempty(version_cache[])
            pkg = read_project_version()
            version_cache[] = if !isempty(pkg)
                pkg
            else
                run_git_or_unknown(`git -C $(dirname(project_toml)) describe --tags --abbrev=0`)
            end
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

    # ===================== Command registry (type-stable dispatch) =====================
    #
    # One singleton type per subcommand. Routing is by multiple dispatch on the concrete
    # type (no Dict{String,Function}, no boxed closures), so each `run_command` method is
    # its own specialization. The only run-time step is mapping the parsed (group, name)
    # strings to the matching singleton once per invocation (command_for); everything after
    # is statically dispatched. `run_command` is declared here and given methods by the top
    # module, where the handlers (in the various submodules) are in scope.

    abstract type Command end

    struct PreprocessDemultiplex <: Command end
    struct DiscoverBlast        <: Command end
    struct DiscoverHsmm         <: Command end
    struct DiscoverSelftest     <: Command end
    struct SearchExact          <: Command end
    struct SearchHeptamer       <: Command end
    struct SearchBwa            <: Command end
    struct AnalyzeCooccurrence  <: Command end
    struct AnalyzeHaplotype     <: Command end
    struct TableOuterjoin       <: Command end
    struct TableLeftjoin        <: Command end
    struct TableTransform       <: Command end
    struct TableAggregate       <: Command end
    struct TableUnique          <: Command end
    struct TableSort            <: Command end
    struct TableFilter          <: Command end
    struct TableSelect          <: Command end
    struct TableFasta           <: Command end
    struct TableCollect         <: Command end
    struct TableExclude         <: Command end
    struct FastaMerge           <: Command end
    struct FastaDiff            <: Command end
    struct FastaHash            <: Command end

    "(group, subcommand) path a Command is reached by on the command line."
    cli_path(::PreprocessDemultiplex) = ("preprocess", "demultiplex")
    cli_path(::DiscoverBlast)        = ("discover", "blast")
    cli_path(::DiscoverHsmm)         = ("discover", "hsmm")
    cli_path(::DiscoverSelftest)     = ("discover", "selftest")
    cli_path(::SearchExact)          = ("search", "exact")
    cli_path(::SearchHeptamer)       = ("search", "heptamer")
    cli_path(::SearchBwa)            = ("search", "bwa")
    cli_path(::AnalyzeCooccurrence)  = ("analyze", "cooccurrence")
    cli_path(::AnalyzeHaplotype)     = ("analyze", "haplotype")
    cli_path(::TableOuterjoin)       = ("table", "outerjoin")
    cli_path(::TableLeftjoin)        = ("table", "leftjoin")
    cli_path(::TableTransform)       = ("table", "transform")
    cli_path(::TableAggregate)       = ("table", "aggregate")
    cli_path(::TableUnique)          = ("table", "unique")
    cli_path(::TableSort)            = ("table", "sort")
    cli_path(::TableFilter)          = ("table", "filter")
    cli_path(::TableSelect)          = ("table", "select")
    cli_path(::TableFasta)           = ("table", "fasta")
    cli_path(::TableCollect)         = ("table", "collect")
    cli_path(::TableExclude)         = ("table", "exclude")
    cli_path(::FastaMerge)           = ("fasta", "merge")
    cli_path(::FastaDiff)            = ("fasta", "diff")
    cli_path(::FastaHash)            = ("fasta", "hash")

    const COMMANDS = (
        PreprocessDemultiplex(),
        DiscoverBlast(), DiscoverHsmm(), DiscoverSelftest(),
        SearchExact(), SearchHeptamer(), SearchBwa(),
        AnalyzeCooccurrence(), AnalyzeHaplotype(),
        TableOuterjoin(), TableLeftjoin(), TableTransform(), TableAggregate(),
        TableUnique(), TableSort(), TableFilter(), TableSelect(), TableFasta(),
        TableCollect(), TableExclude(),
        FastaMerge(), FastaDiff(), FastaHash(),
    )

    """
        command_for(parsed_args) -> Command or nothing

    Resolve the parsed top-level group and its `%COMMAND%` subcommand to the matching
    Command singleton (the single run-time mapping; dispatch is static thereafter).
    """
    function command_for(parsed_args)
        group = get(parsed_args, "%COMMAND%", "")
        group == "" && return nothing
        sub = get(get(parsed_args, group, Dict{String,Any}()), "%COMMAND%", "")
        for c in COMMANDS
            cli_path(c) == (group, sub) && return c
        end
        return nothing
    end

    "Run a resolved command. Methods are defined by the top module (handlers live there)."
    function run_command end

    export Command, command_for, run_command
end
