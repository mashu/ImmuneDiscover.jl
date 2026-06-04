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

A **command** has three small parts that follow the same skeleton everywhere:

1. **Argument table** — in the matching `src/cmd/<group>.jl`, inside `add_<group>_args!(s)`:

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

2. **Command identity** — a singleton type in `cli.jl`, dispatched on rather than looked up:

   ```julia
   struct SearchExact <: Command end
   cli_path(::SearchExact) = ("search", "exact")
   # …and add SearchExact() to the COMMANDS tuple
   ```

3. **Handler binding** — one `run_command` method in `src/immunediscover.jl` that routes the
   command to the library handler (the science stays in the library module):

   ```julia
   Cli.run_command(::Cli.SearchExact, pa) = Exact.handle_exact(pa, immunediscover, Cli.always_gz)
   ```

Routing is multiple dispatch on the concrete `Command` type — no `Dict{String,Function}`,
no boxed closures. The only run-time step is `Cli.command_for`, which maps the parsed
`(group, subcommand)` strings to the matching singleton once per invocation; every call
after that is statically dispatched and precompilable.

## Adding a new subcommand

1. If it's a new group, declare the group in `cli.jl :: add_command_groups!`; add a
   `"<name>" … action = :command` entry to the group list in `add_<group>_args!(s)` and the
   subcommand's `@add_arg_table! s["<group>"]["<name>"] begin … end` block there.
2. In `cli.jl`: add a `struct <Name> <: Command end`, a `cli_path(::<Name>)` method, and the
   instance to `COMMANDS`.
3. In `src/immunediscover.jl`: add `Cli.run_command(::Cli.<Name>, pa) = <Module>.handle_<name>(…)`,
   and write `handle_<name>` in the library module (glue only: parse args → call library →
   write output).
