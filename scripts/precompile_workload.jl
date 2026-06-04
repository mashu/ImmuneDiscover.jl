# Precompile workload for PackageCompiler: exercises CLI and hot paths
# so the compiled binary starts fast. Run during create_app via precompile_execution_file.

using immunediscover
using CSV, DataFrames

# Exercise CLI parsing (ArgParse, dispatch tables)
immunediscover.real_main(["--help"])
immunediscover.real_main(["search", "--help"])
immunediscover.real_main(["analyze", "--help"])
immunediscover.real_main(["table", "--help"])

# Touch DataFrame/CSV path (consistent with @compile_workload in main module)
io = IOBuffer()
write(io, "well\tcase\tname\tgenomic_sequence\n1\tD1\tread1\tATCG\n")
seekstart(io)
CSV.File(io, delim='\t') |> DataFrame

# Exercise the analysis hot paths so their specializations land in the sysimage,
# making the compiled binary fast on first real run. Guarded: a build must never
# fail because of the workload.
try
    tbl = DataFrame(well = [1, 1], case = ["D1", "D1"], name = ["r1", "r2"],
                    genomic_sequence = ["AAAAAAAAAAAAAAACACAGTGCCCCCCCCCC",
                                        "AAAAAAAAAAAAAAACACAGTGCCCCCCCCCC"])
    db = [("IGHV1-1*01", "AAAAAAAAAAAAAAA")]
    immunediscover.Exact.exact_search(tbl, db, "V"; mincount=1, minratio=0.0, N=1)
catch err
    @warn "precompile workload: exact_search exercise skipped" exception=err
end

nothing
