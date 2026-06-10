# CLI and kernel paths exercised during package precompilation and PackageCompiler builds.
# Handlers are run against dummy paths; missing-file errors are ignored via `julia_main`.

const _CLI_PARSE_ARGS = [
    ["discover", "blast", "i.tsv", "d.fa", "o.tsv", "-g", "V"],
    ["discover", "hsmm", "i.tsv", "d.fa", "o.tsv.gz"],
    ["discover", "selftest", "d.full.tsv.gz", "base.fa", "truth.fa", "out.tsv"],
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

const _CLI_DISPATCH_ARGS = [
    ["discover", "blast", "-G", "i.tsv", "d.fa", "o.tsv"],
    ["discover", "hsmm", "i.tsv", "d.fa", "o.tsv.gz"],
    ["discover", "selftest", "d.full.tsv.gz", "base.fa", "truth.fa", "out.tsv"],
    ["search", "exact", "i.tsv", "d.fa", "o.tsv"],
    ["search", "heptamer", "i.tsv", "d.fa", "o.tsv", "s.tsv"],
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

function _with_cli_args(f, args)
    saved = copy(ARGS)
    empty!(ARGS)
    append!(ARGS, args)
    f()
    empty!(ARGS)
    append!(ARGS, saved)
    return nothing
end

"""
    precompile_cli_workload!()

Exercise CSV/DataFrame kernels, ArgParse, and every command handler so native code for
`real_main` / `julia_main` is available at runtime (no first-invocation JIT in the binary).
"""
function precompile_cli_workload!()
    io = IOBuffer()
    write(io, "well\tcase\tname\tgenomic_sequence\n1\tD1\tread1\tATCG\n")
    seekstart(io)
    CSV.File(io, delim='\t') |> DataFrame

    tbl = DataFrame(well = [1, 1], case = ["D1", "D1"], name = ["r1", "r2"],
                    genomic_sequence = ["AAAAAAAAAAAAAAACACAGTGCCCCCCCCCC",
                                        "AAAAAAAAAAAAAAACACAGTGCCCCCCCCCC"])
    db = [("IGHV1-1*01", "AAAAAAAAAAAAAAA")]
    Exact.exact_search(tbl, db, "V"; N=1)

    for args in _CLI_PARSE_ARGS
        parse_commandline(args)
    end
    parse_commandline(["--help"]; exit_after_help=false)
    parse_commandline(["--version"]; exit_after_help=false)

    for args in _CLI_DISPATCH_ARGS
        _with_cli_args(() -> julia_main(), args)
    end
    return nothing
end
