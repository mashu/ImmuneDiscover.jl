module immunediscover
    include("cli.jl")
    include("data.jl")
    include("demultiplex.jl")
    include("simulate.jl")
    include("profile.jl")
    include("exact.jl")
    include("heptamer.jl")
    include("hsmm.jl")
    include("bwa.jl")
    include("blast.jl")
    include("keyedsets.jl")
    include("cooccurrence.jl")
    include("haplotype.jl")
    include("fasta.jl")
    include("merge.jl")
    include("table.jl")

    using .cli
    using .data
    using .demultiplex
    using .simulate
    using .profile
    using .exact
    using .heptamer
    using .HSMM
    using .bwa
    using .blast
    using .keyedsets
    using .cooccurrence
    using .haplotype
    using .fasta
    using .merge
    using .table

    using CSV
    using DataFrames
    using UnicodePlots
    using Glob
    using Statistics
    using DataStructures
    using FASTX
    using PrecompileTools

    export load_fasta, blast_discover

    # Tests import immunediscover.association — alias required until tests are migrated
    const association = cooccurrence

    # Delegated utilities (callers reference immunediscover.concatenate_columns etc.)
    const concatenate_columns = data.concatenate_columns
    const validate_types = data.validate_types
    const get_ratio_threshold = data.get_ratio_threshold

    # --- Command dispatch tables ---

    const SEARCH_HANDLERS = Dict{String, Function}(
        "heptamer" => (pa) -> heptamer.handle_heptamer(pa, immunediscover, cli.always_gz),
        "exact"    => (pa) -> exact.handle_exact(pa, immunediscover, cli.always_gz),
        "hsmm"     => (pa) -> HSMM.handle_hsmm(pa),
        "blast"    => (pa) -> blast.handle_blast(pa, immunediscover, cli.always_gz),
    )

    const ANALYZE_HANDLERS = Dict{String, Function}(
        "cooccurrence" => (pa) -> cooccurrence.handle_cooccurrence(pa, cli.always_gz),
        "haplotype"    => (pa) -> haplotype.handle_haplotype(pa),
        "bwa"          => (pa) -> bwa.handle_bwa(pa, immunediscover, cli.always_gz),
    )

    const FASTA_HANDLERS = Dict{String, Function}(
        "merge" => (pa) -> merge.handle_merge(pa),
        "diff"  => (pa) -> fasta.handle_fasta_diff(pa, immunediscover),
        "hash"  => (pa) -> fasta.handle_fasta_hash(pa, immunediscover),
    )

    const TOPLEVEL_HANDLERS = Dict{String, Function}(
        "demultiplex" => (pa) -> demultiplex.handle_demultiplex(pa, cli.always_gz),
        "table"       => (pa) -> table.handle_table(pa, immunediscover, cli.always_gz),
    )

    # Grouped commands: the Dict key doubles as the parsed_args group key
    const GROUP_HANDLERS = Dict{String, Dict{String, Function}}(
        "search"  => SEARCH_HANDLERS,
        "analyze" => ANALYZE_HANDLERS,
        "fasta"   => FASTA_HANDLERS,
    )

    """
        dispatch_subcommand(parsed_args, group_key, handlers)

    Look up and execute the subcommand handler from a dispatch table.
    """
    function dispatch_subcommand(parsed_args, group_key, handlers)
        subcmd = get(parsed_args[group_key], "%COMMAND%", "")
        handler = get(handlers, subcmd, nothing)
        if handler !== nothing
            handler(parsed_args)
        else
            @warn "Unknown $group_key subcommand: $subcmd"
        end
    end

    """
        real_main(args=[])

    Main entry point — routes top-level and grouped subcommands via dispatch tables.
    """
    function real_main(args=[])
        parsed_args = parse_commandline(args)
        if parsed_args === nothing
            return
        end

        cmd = get(parsed_args, "%COMMAND%", "")

        # Top-level commands
        toplevel = get(TOPLEVEL_HANDLERS, cmd, nothing)
        if toplevel !== nothing
            toplevel(parsed_args)
            return
        end

        # Grouped commands
        handlers = get(GROUP_HANDLERS, cmd, nothing)
        if handlers !== nothing
            dispatch_subcommand(parsed_args, cmd, handlers)
        end
    end

    """
        julia_main()::Cint

    Entry point for PackageCompiler. Top-level catch is the standard CLI boundary pattern.
    """
    function julia_main()::Cint
        try
            real_main(ARGS)
        catch
            Base.invokelatest(Base.display_error, Base.catch_stack())
            return 1
        end
        return 0
    end

    @compile_workload begin
        parse_commandline(["--help"])
        io = IOBuffer()
        write(io, "well\tcase\tname\tgenomic_sequence\n1\tD1\tread1\tATCG\n")
        seekstart(io)
        CSV.File(io, delim='\t') |> DataFrame
    end
end
