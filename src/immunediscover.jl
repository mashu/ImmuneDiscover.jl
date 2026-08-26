module immunediscover
    # --- Utils needed by CLI (Option) then CLI scaffolding ---
    include("utils/option.jl")
    include("cmd/cli.jl")

    # --- Utils (no inter-module deps) ---
    include("utils/gene.jl")
    include("utils/spans.jl")
    include("utils/dna.jl")
    include("utils/data.jl")
    include("utils/ratio_columns.jl")
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
    include("precompile_workload.jl")

    using .Cli
    using .Option
    using .Gene
    using .Spans
    using .DNA
    using .Data
    using .RatioColumns
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
    using PrecompileTools: @setup_workload, @compile_workload

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
        real_main(args=String[])

    Main entry point — parse the command line, resolve the Command, and run it.
    """
    real_main() = real_main(String[])
    function real_main(args::Vector{String})
        parsed_args = parse_commandline(args)
        run_parsed(optional(parsed_args))
        return
    end
    real_main(args::AbstractVector{<:AbstractString}) = real_main(String[a for a in args])

    run_parsed(::Absent) = nothing
    function run_parsed(pa::Present)
        run_resolved(Cli.command_for(pa.value), pa.value)
        return
    end

    run_resolved(::Absent, parsed_args) =
        (@warn "Unknown or missing command: $(get(parsed_args, "%COMMAND%", ""))"; nothing)
    run_resolved(cmd::Present, parsed_args) = Cli.run_command(cmd.value, parsed_args)

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

    # Top-level --help/--version only. Tracing every subcommand or running handlers
    # here produces a Julia 1.12 package image that fails to load.
    @setup_workload begin
        @compile_workload begin
            redirect_stderr(devnull) do
                redirect_stdout(devnull) do
                    parse_commandline(String["--help"]; exit_after_help=false)
                    parse_commandline(String["--version"]; exit_after_help=false)
                end
            end
        end
        Cli.reset_cli_settings!()
    end
end
