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

nothing
