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

From source (already precompiled via PrecompileTools; skip `~/.julia/config/startup.jl`):

```bash
./scripts/run.sh --help
./scripts/run.sh --version
```

## Requirements

- Julia 1.9 or later
- The `build/` project with PackageCompiler (one-time setup):

  ```bash
  julia --project=build -e 'using Pkg; Pkg.instantiate()'
  ```

## How it works

1. **PrecompileTools `@compile_workload`** (in `src/immunediscover.jl`) traces top-level `--help` and `--version` during `Pkg.precompile` so those paths are native code in the package image. A larger fixture CLI is not used here: on Julia 1.12 it yields a cache that fails to load.
2. **PackageCompiler `create_app`** with `scripts/precompile_workload.jl` runs the full fixture CLI (every subcommand `--help`, plus handlers that do not need blastn/BWA) into the standalone binary.

You can also run the build step directly:

```bash
julia --project=build scripts/build_binary.jl [output_dir]
```

**Note:** UnicodePlots is used for terminal bar plots in `demultiplex` and `search exact`. The standalone binary is built without it so that PackageCompiler can succeed; plotting is skipped if UnicodePlots is not available. `--noplot` still disables plots when the package is present.

**CLI snappiness:** `--help` / `--version` do not write `immunediscover.log`. Ordinary commands read the version from `Project.toml` and do not shell out to `git`; `--version` appends the git hash when git is available.
