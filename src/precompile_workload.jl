# CLI and kernel paths exercised during package precompilation and PackageCompiler builds.
# Uses tiny temp fixtures so handlers compile without touching the repo.
#
# Commands that can run from fixtures call `real_main` with no swallowed errors — a silent
# try/catch here used to hide broken traces (the workload compiled nothing). Commands that
# need external binaries (blastn, BWA index) are parse-only plus `precompile(run_command, …)`.
# Help/version for every subcommand is compiled with `exit_after_help=false`.

const DEMUX_HEADER = "well\tcase\tname\tgenomic_sequence\n"

function write_line(path::AbstractString, content::AbstractString)
    open(path, "w") do io
        print(io, content)
    end
    return path
end

function write_tsv(path::AbstractString, header::AbstractString, rows::AbstractVector{<:AbstractString})
    open(path, "w") do io
        print(io, header)
        for row in rows
            println(io, row)
        end
    end
    return path
end

function write_fasta(path::AbstractString)
    return write_line(path, ">IGHV1-1*01\nATCGATCGATCGATCG\n")
end

"""
    precompile_fixture_dir() -> String

Create a temporary directory of minimal inputs for tracing command handlers.
"""
function precompile_fixture_dir()
    root = mktempdir(; cleanup=true)
    p(name) = joinpath(root, name)

    write_tsv(p("demux.tsv"), DEMUX_HEADER, [
        "1\tD1\tr1\tAAAAAAAAAAAAAAACACAGTGCCCCCCCCCC",
        "1\tD1\tr2\tAAAAAAAAAAAAAAACACAGTGCCCCCCCCCC",
    ])
    write_tsv(p("i.tsv"), "col\n", ["1"])
    write_tsv(p("l.tsv"), "key\tval\n", ["1\ta", "2\tb"])
    write_tsv(p("r.tsv"), "key\tother\n", ["1\tx", "3\ty"])
    write_tsv(p("part1.tsv"), "col\n", ["1"])
    write_tsv(p("part2.tsv"), "col\n", ["2"])
    write_tsv(p("idx.tsv"), "forward_index\treverse_index\tcase\n", ["AAAAAAAAAA\tCCCCCCCCCC\tD1"])
    write_tsv(p("seqs.tsv"), "allele_name\tseq\tcase\n",
              ["IGHV1-1*01\tATCGATCG\tD1", "IGHV1-1*02\tGGGGAAAA\tD2"])
    write_tsv(p("exact_like.tsv"),
              "case\tdb_name\tgene\tcount\tsequence\n",
              ["D1\tIGHV1-1*01\tIGHV1-1\t10\tATCGATCG",
               "D1\tIGHV1-1*02\tIGHV1-1\t8\tGGGGAAAA",
               "D2\tIGHV1-1*01\tIGHV1-1\t12\tATCGATCG"])
    write_tsv(p("filt.tsv"), "name\tval\n", ["alpha\t10", "beta\t3"])
    write_line(p("i.fq"), "@r1\nAAAAAAAAAACACAGTGCCCCCCCCC\n+\nIIIIIIIIIIIIIIIIIIIIIIIIII\n")
    write_line(p("d.fa"), ">IGHV1-1*01\nAAAAAAAAAAAAAAA\n")
    write_fasta(p("base.fa"))
    write_line(p("truth.fa"), ">IGHV1-1*01\nATCGATCGATCGATCG\n>IGHV1-1*99\nTTTTAAAACCCCGGGG\n")
    write_line(p("a.fa"), ">IGHV1-1*01\nATCGATCGATCGATCG\n")
    write_line(p("b.fa"), ">IGHV1-2*01\nGGGGAAAACCCCTTTT\n")
    write_fasta(p("i.fa"))
    write_fasta(p("r.fa"))
    write_fasta(p("g.fa"))
    write_line(p("heptamers.json"), "{\"IGHV\":[\"CACAGTG\"]}\n")
    return root
end

"Subcommand --help paths for every registered command (ArgParse, no handler)."
function cli_help_args()
    args = [String["--help"], String["--version"]]
    seen = Set{String}()
    for cmd in Cli.COMMANDS
        group, sub = Cli.cli_path(cmd)
        if !(group in seen)
            push!(args, String[group, "--help"])
            push!(seen, group)
        end
        push!(args, String[group, sub, "--help"])
    end
    return args
end

"Handlers that run to completion on the fixture dir (no external binaries)."
function cli_runnable_args(root::AbstractString)
    p(name) = joinpath(root, name)
    return [
        ["search", "exact", p("demux.tsv"), p("d.fa"), p("out-exact.tsv"),
         "--noplot", "--min-fullcount", "1"],
        ["search", "exact", p("demux.tsv"), p("d.fa"), p("out-exact-diag.tsv"),
         "--noplot", "--min-fullcount", "1", "--diagnostic"],
        ["search", "heptamer", p("demux.tsv"), p("d.fa"), p("out-hept.tsv"),
         p("out-hept-summary.tsv"), "-j", p("heptamers.json")],
        ["analyze", "cooccurrence", p("exact_like.tsv")],
        ["analyze", "haplotype", p("exact_like.tsv"), p("out-haplo.tsv"), "-c", "1"],
        ["preprocess", "demultiplex", p("i.fq"), p("idx.tsv"), p("out-demux.tsv"), "-l", "1"],
        ["table", "outerjoin", p("l.tsv"), p("r.tsv"), p("out-join.tsv"), "-k", "key"],
        ["table", "leftjoin", p("l.tsv"), p("r.tsv"), p("out-ljoin.tsv"), "-k", "key"],
        ["table", "transform", p("seqs.tsv"), p("out-xform.tsv"), "-c", "allele_name", "-p", "(.*)", "-r", "\\1"],
        ["table", "aggregate", p("seqs.tsv"), p("out-agg.tsv"), "-g", "case"],
        ["table", "unique", p("i.tsv"), p("out-uniq.tsv"), "-c", "col"],
        ["table", "sort", p("i.tsv"), p("out-sort.tsv"), "-c", "col"],
        ["table", "filter", p("filt.tsv"), p("out-filter.tsv"), "-c", "name", "--pattern", "a"],
        ["table", "filter", p("filt.tsv"), p("out-filter-num.tsv"), "-c", "val", "--operator", ">=", "--threshold", "5"],
        ["table", "select", p("l.tsv"), p("out-select.tsv"), "-c", "key"],
        ["table", "fasta", p("seqs.tsv"), p("out-export.fa")],
        ["table", "collect", joinpath(root, "part*.tsv"), p("out-collect.tsv")],
        ["table", "exclude", p("seqs.tsv"), p("out-excl.tsv"), p("r.fa")],
        ["fasta", "merge", p("out-merge.fa"), p("a.fa"), p("b.fa")],
        ["fasta", "diff", p("a.fa"), p("b.fa")],
        ["fasta", "hash", p("i.fa")],
    ]
end

"Parse-only: needs blastn, a BWA index, or a discovery full-table that the fixture is not."
function cli_parse_only_args(root::AbstractString)
    p(name) = joinpath(root, name)
    return [
        ["discover", "blast", p("demux.tsv"), p("d.fa"), p("out-blast.tsv"), "-g", "V"],
        ["discover", "blast", "-G", p("demux.tsv"), p("d.fa"), p("out-blast-presets.tsv")],
        ["discover", "hsmm", p("demux.tsv"), p("d.fa"), p("out-hsmm.tsv.gz")],
        ["discover", "selftest", p("exact_like.tsv"), p("base.fa"), p("truth.fa"),
         p("out-selftest.tsv"), "--seq-col", "sequence"],
        ["search", "bwa", p("demux.tsv"), p("out-bwa.tsv"), p("g.fa")],
    ]
end

function with_cli_args(f, args)
    saved = copy(ARGS)
    empty!(ARGS)
    append!(ARGS, args)
    f()
    empty!(ARGS)
    append!(ARGS, saved)
    return nothing
end

precompile_parse!(args::Vector{String}) = (parse_commandline(args); nothing)
precompile_help!(args::Vector{String}) = (parse_commandline(args; exit_after_help=false); nothing)
precompile_run!(args::Vector{String}) = (real_main(args); nothing)

function precompile_library_kernels!()
    io = IOBuffer()
    write(io, DEMUX_HEADER, "1\tD1\tread1\tATCGATCGATCGATCG\n")
    seekstart(io)
    CSV.File(io, delim='\t') |> DataFrame

    tbl = DataFrame(well = [1, 1], case = ["D1", "D1"], name = ["r1", "r2"],
                    genomic_sequence = ["AAAAAAAAAAAAAAACACAGTGCCCCCCCCCC",
                                        "AAAAAAAAAAAAAAACACAGTGCCCCCCCCCC"])
    db = [("IGHV1-1*01", "AAAAAAAAAAAAAAA")]
    Exact.exact_search(tbl, db, "V"; N=1)
    Exact.select_exact_output_columns(Exact.exact_search(tbl, db, "V"; N=1),
                                      Gene.VGene(), absent; diagnostic=false)
    Gene.gene_type_from_name("IGHV1-1*01")
    Gene.gene_type_from_name("GAPDH")
    Blast.consensus_prefix(["AAAC", "AAAT"]; min_fraction=0.5)
    Blast.consensus_suffix(["CGT", "CGA"]; min_fraction=0.5)
    Data.round_floats!(DataFrame(x=[0.123456], n=[1]))
    Data.barplot_if_available(["a", "b"], [3, 1]; title="precompile")
    Data.boxplot_if_available(["a"], [[1, 2, 3]])
    return nothing
end

function precompile_entry_signatures!()
    precompile(real_main, (Vector{String},))
    precompile(julia_main, ())
    precompile(parse_commandline, (Vector{String},))
    precompile(Cli.command_for, (Dict{String,Any},))
    precompile(run_parsed, (Present{Dict{String,Any}},))
    for cmd in Cli.COMMANDS
        precompile(Cli.run_command, (typeof(cmd), Dict{String,Any}))
    end
    return nothing
end

"""
    precompile_cli_workload!(root=precompile_fixture_dir())

Exercise CSV/DataFrame kernels, ArgParse help for every command, and every handler that
can run on fixtures so native code for `real_main` / `julia_main` is in the pkgimage.
"""
function precompile_cli_workload!(root::AbstractString=precompile_fixture_dir())
    redirect_stderr(devnull) do
        redirect_stdout(devnull) do
            Cli.build_argparse_settings()
            precompile_library_kernels!()

            for args in cli_help_args()
                precompile_help!(args)
            end
            for args in cli_parse_only_args(root)
                precompile_parse!(String[args...])
            end
            for args in cli_runnable_args(root)
                precompile_run!(String[args...])
            end

            with_cli_args(["fasta", "hash", joinpath(root, "i.fa")]) do
                julia_main()
            end
            with_cli_args(["--help"]) do
                parse_commandline(ARGS; exit_after_help=false)
            end

            precompile_entry_signatures!()
        end
    end
    return nothing
end
