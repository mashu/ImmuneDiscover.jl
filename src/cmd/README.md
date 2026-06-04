# Command layer (`src/cmd/`)

This directory is the **CLI/command layer**. It is deliberately kept separate from the
**reusable library** under `src/{search,discover,analyze,preprocess,utils}/`, which
contains the algorithms and has *no* dependency on ArgParse.

```
src/
├── immunediscover.jl     # top module: includes everything, routes subcommands to handlers
├── cmd/                  # ← you are here: the command/CLI layer
│   ├── cli.jl            #    Cli module: ArgParseSettings skeleton, presets, parse_commandline
│   ├── preprocess.jl     #    arg tables for the `preprocess` group
│   ├── discover.jl       #    arg tables for the `discover` group (blast, hsmm)
│   ├── search.jl         #    arg tables for the `search` group (exact, heptamer, bwa)
│   ├── analyze.jl        #    arg tables for the `analyze` group (cooccurrence, haplotype)
│   ├── table.jl          #    arg tables for the `table` group
│   └── fasta.jl          #    arg tables for the `fasta` group
└── search/exact.jl …     # reusable library module: exposes exact_search(...) + handle_exact(...)
```

A **command** has two parts that follow the same skeleton everywhere:

1. **Argument table** — lives in the matching `src/cmd/<group>.jl`, inside
   `add_<group>_args!(s)`:

   ```julia
   @add_arg_table! s["<group>"]["<name>"] begin
       "input"
           help = "…"
           required = true
       "-x", "--option"
           help = "…"
           default = …
           arg_type = …
           range_tester = (x -> …)
   end
   ```

2. **Handler** — `handle_<name>(parsed_args, immunediscover_module, always_gz)` in the
   library module that owns the algorithm (e.g. `Exact.handle_exact`). The handler reads
   `parsed_args[<group>][<name>]`, calls the library function, and writes output.

## Adding a new subcommand

1. Declare it under its group in `cli.jl :: add_command_groups!` (if it's a new group) and
   add a `… action = :command` entry to the group list inside `add_<group>_args!(s)`.
2. Add its `@add_arg_table! s["<group>"]["<name>"] begin … end` block to `src/cmd/<group>.jl`.
3. Write `handle_<name>(...)` in the relevant library module (keep the science in the library;
   the handler is only glue: parse args → call library → write output).
4. Register the handler in the dispatch table in `src/immunediscover.jl`
   (`<GROUP>_HANDLERS`).

## Direction

The next refactor (tracked separately, CI-gated) folds steps 1–4 into a single
self-contained command file per subcommand via a `Subcommand` registry
`(group, name, configure!, handler)`, so the parser and the dispatch are both derived
from one declaration and adding a command touches exactly one file.
