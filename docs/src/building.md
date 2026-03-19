# Building a standalone binary

Standalone binaries use **PackageCompiler** and **PrecompileTools** for fast startup. Build on each target platform (the binary is not cross-compiled: build on Linux for Linux, macOS for macOS, etc.).

## Quick build

From the repository root:

```bash
./scripts/build_binary.sh
```

Output: `build/immunediscover_app/bin/immunediscover`. Copy the entire `build/immunediscover_app` directory to relocate the app.

To choose a different output directory:

```bash
./scripts/build_binary.sh /path/to/output_dir
```

## Requirements

- Julia 1.9 or later
- The `build/` project with PackageCompiler (one-time setup):

  ```bash
  julia --project=build -e 'using Pkg; Pkg.instantiate()'
  ```

## How it works

The script `scripts/build_binary.jl` calls `PackageCompiler.create_app()` on the immunediscover package. The optional file `scripts/precompile_workload.jl` is used as a precompile execution script: it runs `--help` and a few CLI paths so the compiled sysimage already has those code paths compiled, reducing startup time.

You can also run the build step directly:

```bash
julia --project=build scripts/build_binary.jl [output_dir]
```
