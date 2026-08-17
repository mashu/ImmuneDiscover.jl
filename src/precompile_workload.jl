# CLI and kernel paths exercised during package precompilation and PackageCompiler builds.
# Uses tiny temp fixtures so handlers compile without touching the repo or printing errors.
# The compile workload must call `real_main(::Vector{String})` (the live CLI signature)
# without NullLogger, so Logging / ArgParse / handler native code is in the pkgimage.

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

    write_tsv(p("demux.tsv"), DEMUX_HEADER, ["1\tD1\tr1\tATCGATCGATCGATCGATCG"])
    write_tsv(p("i.tsv"), "col\n", ["1"])
    write_tsv(p("l.tsv"), "key\tval\n", ["1\ta", "2\tb"])
    write_tsv(p("r.tsv"), "key\tother\n", ["1\tx", "3\ty"])
    write_tsv(p("part1.tsv"), "col\n", ["1"])
    write_tsv(p("part2.tsv"), "col\n", ["2"])
    write_tsv(p("idx.tsv"), "forward_index\treverse_index\tcase\n", ["ATCGATCGAT\tGCTAGCTAGC\tD1"])
    write_tsv(p("seqs.tsv"), "allele_name\tseq\tcase\n",
              ["IGHV1-1*01\tATCGATCG\tD1", "IGHV1-1*02\tGGGGAAAA\tD2"])
    write_tsv(p("exact_like.tsv"),
              "case\tdb_name\tgene\tcount\tsequence\n",
              ["D1\tIGHV1-1*01\tIGHV1-1\t10\tATCGATCG",
               "D1\tIGHV1-1*02\tIGHV1-1\t8\tGGGGAAAA",
               "D2\tIGHV1-1*01\tIGHV1-1\t12\tATCGATCG"])
    write_tsv(p("filt.tsv"), "name\tval\n", ["alpha\t10", "beta\t3"])
    write_line(p("i.fq"), "@r1\nATCGATCGATCGATCGATCG\n+\nIIIIIIIIIIIIIIIIIIII\n")
    write_fasta(p("d.fa"))
    write_fasta(p("base.fa"))
    write_line(p("truth.fa"), ">IGHV1-1*01\nATCGATCGATCGATCG\n>IGHV1-1*99\nTTTTAAAACCCCGGGG\n")
    write_line(p("a.fa"), ">IGHV1-1*01\nATCGATCGATCGATCG\n")
    write_line(p("b.fa"), ">IGHV1-2*01\nGGGGAAAACCCCTTTT\n")
    write_fasta(p("i.fa"))
    write_fasta(p("r.fa"))
    write_fasta(p("g.fa"))
    return root
end

function cli_command_args(root::AbstractString)
    p(name) = joinpath(root, name)
    return [
        ["discover", "blast", p("demux.tsv"), p("d.fa"), p("out-blast.tsv"), "-g", "V"],
        ["discover", "blast", "-G", p("demux.tsv"), p("d.fa"), p("out-blast-presets.tsv")],
        ["discover", "hsmm", p("demux.tsv"), p("d.fa"), p("out-hsmm.tsv.gz")],
        ["discover", "selftest", p("exact_like.tsv"), p("base.fa"), p("truth.fa"), p("out-selftest.tsv"), "--seq-col", "sequence"],
        ["search", "exact", p("demux.tsv"), p("d.fa"), p("out-exact.tsv"), "--noplot"],
        ["search", "heptamer", p("demux.tsv"), p("d.fa"), p("out-hept.tsv"), p("out-hept-summary.tsv")],
        ["search", "bwa", p("demux.tsv"), p("out-bwa.tsv"), p("g.fa")],
        ["analyze", "cooccurrence", p("exact_like.tsv")],
        ["analyze", "haplotype", p("exact_like.tsv"), p("out-haplo.tsv"), "-c", "1"],
        ["preprocess", "demultiplex", p("i.fq"), p("idx.tsv"), p("out-demux.tsv")],
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

function with_cli_args(f, args)
    saved = copy(ARGS)
    empty!(ARGS)
    append!(ARGS, args)
    f()
    empty!(ARGS)
    append!(ARGS, saved)
    return nothing
end

"Compile a command through the live CLI entry (`real_main(::Vector{String})`)."
function precompile_run!(args::Vector{String})
    try
        real_main(args)
    catch
    end
    return nothing
end

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
    Gene.gene_type_from_name("IGHV1-1*01")
    Gene.gene_type_from_name("GAPDH")
    Blast.consensus_prefix(["AAAC", "AAAT"]; min_fraction=0.5)
    Blast.consensus_suffix(["CGT", "CGA"]; min_fraction=0.5)
    Data.round_floats!(DataFrame(x=[0.123456], n=[1]))
    return nothing
end

function precompile_entry_signatures!()
    precompile(real_main, (Vector{String},))
    precompile(julia_main, ())
    precompile(parse_commandline, (Vector{String},))
    precompile(Cli.command_for, (Dict{String,Any},))
    precompile(run_parsed, (Present{Dict{String,Any}},))
    precompile(Cli.run_command, (Cli.FastaHash, Dict{String,Any}))
    precompile(Cli.run_command, (Cli.FastaDiff, Dict{String,Any}))
    precompile(Cli.run_command, (Cli.FastaMerge, Dict{String,Any}))
    precompile(Cli.run_command, (Cli.TableUnique, Dict{String,Any}))
    precompile(Cli.run_command, (Cli.TableSort, Dict{String,Any}))
    precompile(Cli.run_command, (Cli.TableFilter, Dict{String,Any}))
    precompile(Cli.run_command, (Cli.TableSelect, Dict{String,Any}))
    precompile(Cli.run_command, (Cli.TableOuterjoin, Dict{String,Any}))
    precompile(Cli.run_command, (Cli.TableLeftjoin, Dict{String,Any}))
    precompile(Cli.run_command, (Cli.SearchExact, Dict{String,Any}))
    precompile(Cli.run_command, (Cli.AnalyzeHaplotype, Dict{String,Any}))
    precompile(Cli.run_command, (Cli.AnalyzeCooccurrence, Dict{String,Any}))
    return nothing
end

"""
    precompile_cli_workload!(root=precompile_fixture_dir())

Exercise CSV/DataFrame kernels, ArgParse, and every command handler so native code for
`real_main` / `julia_main` is available at runtime (no first-invocation JIT in the binary).
"""
function precompile_cli_workload!(root::AbstractString=precompile_fixture_dir())
    redirect_stderr(devnull) do
        redirect_stdout(devnull) do
            Cli.build_argparse_settings()
            precompile_library_kernels!()

            parse_commandline(String["--help"]; exit_after_help=false)
            parse_commandline(String["--version"]; exit_after_help=false)
            parse_commandline(String["fasta", "hash", "x"])
            parse_commandline(String["table", "unique", "a.tsv", "b.tsv", "-c", "col"])
            parse_commandline(String["search", "exact", "a.tsv", "b.fa", "c.tsv"])
            parse_commandline(String["discover", "blast", "a.tsv", "b.fa", "c.tsv", "-g", "V"])

            for args in cli_command_args(root)
                parse_commandline(String[args...])
                precompile_run!(String[args...])
            end

            with_cli_args(["fasta", "hash", joinpath(root, "i.fa")]) do
                julia_main()
            end

            precompile_entry_signatures!()
        end
    end
    return nothing
end
