# CLI and kernel paths exercised during package precompilation and PackageCompiler builds.
# Uses tiny temp fixtures so handlers compile without touching the repo or printing errors.

using Logging: NullLogger, with_logger

const _DEMUX_HEADER = "well\tcase\tname\tgenomic_sequence\n"

function _write_line(path::AbstractString, content::AbstractString)
    open(path, "w") do io
        print(io, content)
    end
    return path
end

function _write_tsv(path::AbstractString, header::AbstractString, rows::AbstractVector{<:AbstractString})
    open(path, "w") do io
        print(io, header)
        for row in rows
            println(io, row)
        end
    end
    return path
end

function _write_fasta(path::AbstractString)
    return _write_line(path, ">ref\nATCGATCG\n")
end

"""
    precompile_fixture_dir() -> String

Create a temporary directory of minimal inputs for tracing command handlers.
"""
function precompile_fixture_dir()
    root = mktempdir(; cleanup=true)
    p(name) = joinpath(root, name)

    _write_tsv(p("demux.tsv"), _DEMUX_HEADER, ["1\tD1\tr1\tATCGATCGATCG"])
    _write_tsv(p("i.tsv"), "col\n", ["1"])
    _write_tsv(p("l.tsv"), "key\n", ["1"])
    _write_tsv(p("r.tsv"), "key\n", ["1"])
    _write_tsv(p("part1.tsv"), "col\n", ["1"])
    _write_tsv(p("part2.tsv"), "col\n", ["2"])
    _write_tsv(p("idx.tsv"), "forward_index\treverse_index\tcase\n", ["ATCGATCGAT\tGCTAGCTAGC\tD1"])
    _write_line(p("i.fq"), "@r1\nATCGATCGATCGATCGATCG\n+\nIIIIIIIIIIIIIIIIIIII\n")
    for fa in ("d.fa", "base.fa", "truth.fa", "a.fa", "b.fa", "i.fa", "r.fa", "g.fa")
        _write_fasta(p(fa))
    end
    return root
end

function _cli_parse_args(root::AbstractString)
    p(name) = joinpath(root, name)
    return [
        ["discover", "blast", p("demux.tsv"), p("d.fa"), p("out-blast.tsv"), "-g", "V"],
        ["discover", "hsmm", p("demux.tsv"), p("d.fa"), p("out-hsmm.tsv.gz")],
        ["discover", "selftest", p("demux.tsv"), p("base.fa"), p("truth.fa"), p("out-selftest.tsv")],
        ["search", "exact", p("demux.tsv"), p("d.fa"), p("out-exact.tsv")],
        ["search", "heptamer", p("demux.tsv"), p("d.fa"), p("out-hept.tsv"), p("out-hept-summary.tsv")],
        ["search", "bwa", p("demux.tsv"), p("out-bwa.tsv"), p("g.fa")],
        ["analyze", "cooccurrence", p("i.tsv")],
        ["analyze", "haplotype", p("i.tsv"), p("out-haplo.tsv")],
        ["preprocess", "demultiplex", p("i.fq"), p("idx.tsv"), p("out-demux.tsv")],
        ["table", "outerjoin", p("l.tsv"), p("r.tsv"), p("out-join.tsv"), "-k", "key"],
        ["table", "leftjoin", p("l.tsv"), p("r.tsv"), p("out-ljoin.tsv"), "-k", "key"],
        ["table", "transform", p("i.tsv"), p("out-xform.tsv"), "-c", "col", "-p", "(.*)", "-r", "\\1"],
        ["table", "aggregate", p("i.tsv"), p("out-agg.tsv"), "-g", "col"],
        ["table", "unique", p("i.tsv"), p("out-uniq.tsv"), "-c", "col"],
        ["table", "sort", p("i.tsv"), p("out-sort.tsv"), "-c", "col"],
        ["table", "filter", p("i.tsv"), p("out-filter.tsv"), "-c", "col", "--pattern", "x"],
        ["table", "select", p("i.tsv"), p("out-select.tsv"), "-c", "col"],
        ["table", "fasta", p("i.tsv"), p("out-export.fa")],
        ["table", "collect", joinpath(root, "part*.tsv"), p("out-collect.tsv")],
        ["table", "exclude", p("i.tsv"), p("out-excl.tsv"), p("r.fa")],
        ["fasta", "merge", p("out-merge.fa"), p("a.fa"), p("b.fa")],
        ["fasta", "diff", p("a.fa"), p("b.fa")],
        ["fasta", "hash", p("i.fa")],
    ]
end

function _cli_dispatch_args(root::AbstractString)
    p(name) = joinpath(root, name)
    return [
        ["discover", "blast", "-G", p("demux.tsv"), p("d.fa"), p("out-blast.tsv")],
        ["discover", "hsmm", p("demux.tsv"), p("d.fa"), p("out-hsmm.tsv.gz")],
        ["discover", "selftest", p("demux.tsv"), p("base.fa"), p("truth.fa"), p("out-selftest.tsv")],
        ["search", "exact", p("demux.tsv"), p("d.fa"), p("out-exact.tsv")],
        ["search", "heptamer", p("demux.tsv"), p("d.fa"), p("out-hept.tsv"), p("out-hept-summary.tsv")],
        ["search", "bwa", p("demux.tsv"), p("out-bwa.tsv"), p("g.fa")],
        ["analyze", "cooccurrence", p("i.tsv")],
        ["analyze", "haplotype", p("i.tsv"), p("out-haplo.tsv")],
        ["preprocess", "demultiplex", p("i.fq"), p("idx.tsv"), p("out-demux.tsv")],
        ["table", "outerjoin", p("l.tsv"), p("r.tsv"), p("out-join.tsv"), "-k", "key"],
        ["table", "leftjoin", p("l.tsv"), p("r.tsv"), p("out-ljoin.tsv"), "-k", "key"],
        ["table", "transform", p("i.tsv"), p("out-xform.tsv"), "-c", "col", "-p", "(.*)", "-r", "\\1"],
        ["table", "aggregate", p("i.tsv"), p("out-agg.tsv"), "-g", "col"],
        ["table", "unique", p("i.tsv"), p("out-uniq.tsv"), "-c", "col"],
        ["table", "sort", p("i.tsv"), p("out-sort.tsv"), "-c", "col"],
        ["table", "filter", p("i.tsv"), p("out-filter.tsv"), "-c", "col", "--pattern", "x"],
        ["table", "select", p("i.tsv"), p("out-select.tsv"), "-c", "col"],
        ["table", "fasta", p("i.tsv"), p("out-export.fa")],
        ["table", "collect", joinpath(root, "part*.tsv"), p("out-collect.tsv")],
        ["table", "exclude", p("i.tsv"), p("out-excl.tsv"), p("r.fa")],
        ["fasta", "merge", p("out-merge.fa"), p("a.fa"), p("b.fa")],
        ["fasta", "diff", p("a.fa"), p("b.fa")],
        ["fasta", "hash", p("i.fa")],
    ]
end

function _with_cli_args(f, args)
    saved = copy(ARGS)
    empty!(ARGS)
    append!(ARGS, args)
    f()
    empty!(ARGS)
    append!(ARGS, saved)
    return nothing
end

# Build-time only: route through julia_main (catches errors) with output discarded.
function _precompile_run!(args)
    redirect_stderr(devnull) do
        redirect_stdout(devnull) do
            _with_cli_args(() -> julia_main(), args)
        end
    end
    return nothing
end

"""
    precompile_cli_workload!()

Exercise CSV/DataFrame kernels, ArgParse, and every command handler so native code for
`real_main` / `julia_main` is available at runtime (no first-invocation JIT in the binary).
"""
function precompile_cli_workload!()
    io = IOBuffer()
    write(io, _DEMUX_HEADER, "1\tD1\tread1\tATCG\n")
    seekstart(io)
    CSV.File(io, delim='\t') |> DataFrame

    tbl = DataFrame(well = [1, 1], case = ["D1", "D1"], name = ["r1", "r2"],
                    genomic_sequence = ["AAAAAAAAAAAAAAACACAGTGCCCCCCCCCC",
                                        "AAAAAAAAAAAAAAACACAGTGCCCCCCCCCC"])
    db = [("IGHV1-1*01", "AAAAAAAAAAAAAAA")]
    Exact.exact_search(tbl, db, "V"; N=1)

    root = precompile_fixture_dir()

    with_logger(NullLogger()) do
        for args in _cli_parse_args(root)
            parse_commandline(args)
        end
        redirect_stdout(devnull) do
            parse_commandline(["--help"]; exit_after_help=false)
            parse_commandline(["--version"]; exit_after_help=false)
        end
        for args in _cli_dispatch_args(root)
            _precompile_run!(args)
        end
        saved = copy(ARGS)
        try
            empty!(ARGS)
            append!(ARGS, ["discover", "blast", "-G",
                           joinpath(root, "demux.tsv"), joinpath(root, "d.fa"),
                           joinpath(root, "out-blast.tsv")])
            redirect_stdout(devnull) do
                julia_main()
            end
        finally
            empty!(ARGS)
            append!(ARGS, saved)
        end
    end
    return nothing
end
