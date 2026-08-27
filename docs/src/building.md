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

`--help` / `--version` (including `search exact --help`) are printed from ArgParse text
cached under `build/help/` (gitignored). `scripts/run.sh` regenerates that cache when
`src/cmd/` is newer. The first `--help` after a clone or CLI change starts Julia once;
later help does not.

## Faster local commands (sysimage)

This compiles immunediscover **and its dependencies** (CSV, DataFrames, FASTX, …) into one
native image for this machine:

```bash
./scripts/build_sysimage.sh
```

Output: `build/immunediscover.so`. `scripts/run.sh` uses it automatically. Rebuild after
changing source or `Manifest.toml`. This is CPU-native and not for distributing binaries.

## Standalone binary (PackageCompiler)

```bash
./scripts/build_binary.sh
```

## Requirements

- Julia 1.9 or later
- The `build/` project with PackageCompiler (one-time setup):

  ```bash
  julia --project=build -e 'using Pkg; Pkg.instantiate()'
  ```

## How it works

1. **PrecompileTools `@compile_workload`** (in `src/immunediscover.jl`) traces top-level `--help` and `--version` during `Pkg.precompile` so those paths are native code in the package image. A larger fixture CLI is not used here: on Julia 1.12 it yields a cache that fails to load.
2. **PackageCompiler `create_sysimage`** (`scripts/build_sysimage.sh`) bakes the package and every dependency into `build/immunediscover.so` for fast `scripts/run.sh` on this machine.
3. **PackageCompiler `create_app`** with `scripts/precompile_workload.jl` bakes handlers into a relocatable sysimage. After the build, a shell wrapper around the trampoline serves `--help` / `--version` from ArgParse pages generated into `build/help/` so those paths do not load the sysimage.

### Juliac.jl

[Juliac](https://github.com/JuliaLang/Juliac.jl) is a CLI in front of PackageCompiler, plus optional `--trim` on Julia 1.12+ to drop unreachable code. **`--trim` does not work for this package.** Trim needs the compiler to see every call from `@main`. DataFrames, CSV, ArgParse, and UnicodePlots are deliberately type-unstable (dynamic column schemas, argparse tables). Juliac then fails with verifier errors such as unresolved `DataFrames` calls. Without `--trim`, Juliac produces the same class of large bundled executable as `create_app` — which we already build. The package defines `@main` so an untrimmed Juliac frontend can be tried later; do not expect a small static binary until those dependencies themselves trim.

You can also run the app build step directly:

```bash
julia --project=build scripts/build_binary.jl [output_dir]
```

**Note:** UnicodePlots is used for terminal bar plots in `demultiplex` and `search exact`. The standalone binary is built without it so that PackageCompiler can succeed; plotting is skipped if UnicodePlots is not available. `--noplot` still disables plots when the package is present.

**CLI snappiness:** `--version` never starts Julia. `--help` regenerates ArgParse pages into
`build/help/` when `src/cmd/` is newer, then prints them without Julia. Ordinary commands
read the version from `Project.toml` and do not shell out to `git`; `--version` appends the
git hash when git is available. `--help` / `--version` also do not write `immunediscover.log`.
