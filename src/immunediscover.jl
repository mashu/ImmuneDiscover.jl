module immunediscover
    # --- CLI scaffolding ---
    include("cmd/cli.jl")

    # --- Utils (no inter-module deps) ---
    include("utils/data.jl")
    include("utils/filters.jl")
    include("utils/seqstats.jl")
    include("utils/align.jl")
    include("utils/mosaic.jl")
    include("utils/report.jl")
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
    include("discover/selftest.jl")

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
    using .SeqStats
    using .Mosaic
    using .Report
    using .Demultiplex
    using .Simulate
    using .Profile
    using .Exact
    using .Heptamer
    using .HSMM
    using .Bwa
    using .Blast
    using .Selftest
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

    # --- Command dispatch: one method per subcommand, dispatched on the Cli.Command
    #     singleton (defined in cmd/cli.jl). Handlers live in the submodules. ---

    Cli.run_command(::Cli.PreprocessDemultiplex, pa) = Demultiplex.handle_demultiplex(pa, Cli.always_gz)
    Cli.run_command(::Cli.DiscoverBlast, pa)        = Blast.handle_blast(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.DiscoverHsmm, pa)         = HSMM.handle_hsmm(pa)
    Cli.run_command(::Cli.DiscoverSelftest, pa)     = Selftest.handle_selftest(pa)
    Cli.run_command(::Cli.SearchExact, pa)          = Exact.handle_exact(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.SearchHeptamer, pa)       = Heptamer.handle_heptamer(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.SearchBwa, pa)            = Bwa.handle_bwa(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.AnalyzeCooccurrence, pa)  = Cooccurrence.handle_cooccurrence(pa)
    Cli.run_command(::Cli.AnalyzeHaplotype, pa)     = Haplotype.handle_haplotype(pa)
    Cli.run_command(::Cli.TableOuterjoin, pa)       = Table.handle_outerjoin(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.TableLeftjoin, pa)        = Table.handle_leftjoin(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.TableTransform, pa)       = Table.handle_transform(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.TableAggregate, pa)       = Table.handle_aggregate(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.TableUnique, pa)          = Table.handle_unique(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.TableSort, pa)            = Table.handle_sort(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.TableFilter, pa)          = Table.handle_filter(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.TableSelect, pa)          = Table.handle_select(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.TableFasta, pa)           = Table.handle_fasta_export(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.TableCollect, pa)         = Table.handle_collect(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.TableExclude, pa)         = Table.handle_exclude(pa, immunediscover, Cli.always_gz)
    Cli.run_command(::Cli.FastaMerge, pa)           = Merge.handle_merge(pa)
    Cli.run_command(::Cli.FastaDiff, pa)            = Fasta.handle_fasta_diff(pa, immunediscover)
    Cli.run_command(::Cli.FastaHash, pa)            = Fasta.handle_fasta_hash(pa, immunediscover)

    """
        real_main(args=[])

    Main entry point — parse the command line, resolve the Command, and run it.
    """
    function real_main(args=[])
        parsed_args = parse_commandline(args)
        parsed_args === nothing && return
        cmd = Cli.command_for(parsed_args)
        if cmd === nothing
            @warn "Unknown or missing command: $(get(parsed_args, "%COMMAND%", ""))"
            return
        end
        Cli.run_command(cmd, parsed_args)
        return
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
        df = CSV.File(io, delim='\t') |> DataFrame

        for args in [
            ["discover", "blast", "i.tsv", "d.fa", "o.tsv", "-g", "V"],
            ["discover", "hsmm", "i.tsv", "d.fa", "o.tsv.gz"],
            ["search", "exact", "i.tsv.gz", "d.fa", "o.tsv.gz"],
            ["search", "heptamer", "i.tsv.gz", "d.fa", "o.tsv.gz", "s.tsv"],
            ["search", "bwa", "i.tsv", "o.tsv", "g.fa"],
            ["analyze", "cooccurrence", "i.tsv"],
            ["analyze", "haplotype", "i.tsv", "o.tsv"],
            ["preprocess", "demultiplex", "i.fq", "idx.tsv", "o.tsv"],
            ["table", "outerjoin", "l.tsv", "r.tsv", "o.tsv", "-k", "key"],
            ["table", "leftjoin", "l.tsv", "r.tsv", "o.tsv", "-k", "key"],
            ["table", "transform", "i.tsv", "o.tsv", "-c", "col", "-p", "(.*)", "-r", "\\1"],
            ["table", "aggregate", "i.tsv", "o.tsv", "-g", "col"],
            ["table", "unique", "i.tsv", "o.tsv", "-c", "col"],
            ["table", "sort", "i.tsv", "o.tsv", "-c", "col"],
            ["table", "filter", "i.tsv", "o.tsv", "-c", "col", "--pattern", "x"],
            ["table", "select", "i.tsv", "o.tsv", "-c", "col"],
            ["table", "fasta", "i.tsv", "o.fa"],
            ["table", "collect", "*.tsv", "o.tsv"],
            ["table", "exclude", "i.tsv", "o.tsv", "r.fa"],
            ["fasta", "merge", "o.fa", "a.fa", "b.fa"],
            ["fasta", "diff", "a.fa", "b.fa"],
            ["fasta", "hash", "i.fa"],
        ]
            parse_commandline(args)
        end
    end
end
