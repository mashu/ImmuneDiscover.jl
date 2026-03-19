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

    using .Cli
    using .Data
    using .Demultiplex
    using .Simulate
    using .Profile
    using .Exact
    using .Heptamer
    using .HSMM
    using .Bwa
    using .Blast
    using .KeyedSets
    using .Cooccurrence
    using .Haplotype
    using .Fasta
    using .Merge
    using .Table

    using CSV
    using DataFrames
    using Glob
    using Statistics
    using DataStructures
    using FASTX
    using PrecompileTools

    export load_fasta, blast_discover

    # Tests import immunediscover.association — alias required until tests are migrated
    const association = Cooccurrence

    # Delegated utilities (callers reference immunediscover.concatenate_columns etc.)
    const concatenate_columns = Data.concatenate_columns
    const validate_types = Data.validate_types
    const get_ratio_threshold = Data.get_ratio_threshold

    # --- Command dispatch tables ---

    const SEARCH_HANDLERS = Dict{String, Function}(
        "heptamer" => (pa) -> Heptamer.handle_heptamer(pa, immunediscover, Cli.always_gz),
        "exact"    => (pa) -> Exact.handle_exact(pa, immunediscover, Cli.always_gz),
        "hsmm"     => (pa) -> HSMM.handle_hsmm(pa),
        "blast"    => (pa) -> Blast.handle_blast(pa, immunediscover, Cli.always_gz),
    )

    const ANALYZE_HANDLERS = Dict{String, Function}(
        "cooccurrence" => (pa) -> Cooccurrence.handle_cooccurrence(pa, Cli.always_gz),
        "haplotype"    => (pa) -> Haplotype.handle_haplotype(pa),
        "bwa"          => (pa) -> Bwa.handle_bwa(pa, immunediscover, Cli.always_gz),
    )

    const FASTA_HANDLERS = Dict{String, Function}(
        "merge" => (pa) -> Merge.handle_merge(pa),
        "diff"  => (pa) -> Fasta.handle_fasta_diff(pa, immunediscover),
        "hash"  => (pa) -> Fasta.handle_fasta_hash(pa, immunediscover),
    )

    const TOPLEVEL_HANDLERS = Dict{String, Function}(
        "demultiplex" => (pa) -> Demultiplex.handle_demultiplex(pa, Cli.always_gz),
        "table"       => (pa) -> Table.handle_table(pa, immunediscover, Cli.always_gz),
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
        io = IOBuffer()
        write(io, "well\tcase\tname\tgenomic_sequence\n1\tD1\tread1\tATCG\n")
        seekstart(io)
        CSV.File(io, delim='\t') |> DataFrame
    end
end
