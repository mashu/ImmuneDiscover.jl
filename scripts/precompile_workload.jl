# Precompile workload for PackageCompiler: exercises CLI parsing and the analysis
# hot paths so the compiled binary starts fast. Run during create_app via
# precompile_execution_file.
#
# NOTE: do NOT use real_main(["--help"]) here — ArgParse's help action calls exit(0)
# (exit_after_help defaults to true), which would terminate this script before any
# later statement is traced. We trace CLI parsing with parse_commandline (no exit,
# no file IO) and run the real kernels directly instead.

using immunediscover
using CSV, DataFrames

# DataFrame/CSV path (consistent with @compile_workload in the main module).
io = IOBuffer()
write(io, "well\tcase\tname\tgenomic_sequence\n1\tD1\tread1\tATCG\n")
seekstart(io)
CSV.File(io, delim='\t') |> DataFrame

# Analysis hot path: a tiny end-to-end exact search so its specializations land in
# the sysimage. Guarded — a build must never fail because of the workload.
try
    tbl = DataFrame(well = [1, 1], case = ["D1", "D1"], name = ["r1", "r2"],
                    genomic_sequence = ["AAAAAAAAAAAAAAACACAGTGCCCCCCCCCC",
                                        "AAAAAAAAAAAAAAACACAGTGCCCCCCCCCC"])
    db = [("IGHV1-1*01", "AAAAAAAAAAAAAAA")]
    immunediscover.Exact.exact_search(tbl, db, "V"; N=1)
catch err
    @warn "precompile workload: exact_search exercise skipped" exception=err
end

# Trace ArgParse / dispatch-table parsing across representative commands (no exit).
for args in (
        ["search", "exact", "i.tsv", "d.fa", "o.tsv"],
        ["discover", "blast", "i.tsv", "d.fa", "o.tsv", "-g", "V"],
        ["discover", "hsmm", "i.tsv", "d.fa", "o.tsv.gz"],
        ["analyze", "cooccurrence", "i.tsv"],
        ["analyze", "haplotype", "i.tsv", "o.tsv"],
        ["preprocess", "demultiplex", "i.fq", "idx.tsv", "o.tsv"],
        ["table", "sort", "i.tsv", "o.tsv", "-c", "col"],
        ["fasta", "merge", "o.fa", "a.fa", "b.fa"],
    )
    immunediscover.parse_commandline(args)
end

nothing
