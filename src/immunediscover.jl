module immunediscover
    # --- CLI scaffolding ---
    include("cli.jl")

    # --- Utils (no inter-module deps) ---
    include("utils/data.jl")
    include("utils/filters.jl")
    include("utils/keyedsets.jl")
    include("discover/profile.jl")

    # --- Search / discover (may depend on utils) ---
    include("search/exact.jl")
    include("preprocess/demultiplex.jl")
    include("preprocess/simulate.jl")
    include("search/heptamer.jl")
    include("discover/hsmm.jl")
    include("search/bwa.jl")
    include("discover/blast.jl")

    # --- Analyze (may depend on search) ---
    include("analyze/cooccurrence.jl")
    include("analyze/haplotype.jl")

    # --- Table / FASTA utilities ---
    include("utils/fasta.jl")
    include("utils/merge.jl")
    include("utils/table.jl")

    using .Cli
    using .Data
    using .Filters
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

        toplevel = get(TOPLEVEL_HANDLERS, cmd, nothing)
        if toplevel !== nothing
            toplevel(parsed_args)
            return
        end

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
