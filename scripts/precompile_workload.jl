# Precompile workload for PackageCompiler: exercises CLI parsing and every command handler
# so the compiled binary starts without first-invocation JIT. Run during create_app via
# precompile_execution_file.

using immunediscover

immunediscover.precompile_cli_workload!()

nothing
