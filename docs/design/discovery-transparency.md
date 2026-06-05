# Design: Discovery transparency, discriminative metrics, and self-test

Status: proposed (implement as a separate PR series, Phase 1 first).

Goal: make `discover blast` explain *why* every candidate is kept or rejected, emit a full
annotated table next to the filtered one, print fancy staged diagnostics, add metrics that
better separate true novel alleles from artifacts, and provide a `selftest` command to
measure (and tune) recovery of known-vs-novel alleles.

The design generalizes to `search exact` and `discover hsmm`, but Phase 1 targets
`discover blast`.

---

## Phase 1 — Transparency: reject reasons, dual output, staged reports

### 1.1 Filters: an annotate path beside the drop path

In `src/utils/filters.jl`, add (keeping `GermlineFilter` for the existing drop path):

```julia
"""
    annotate_rejections!(df, criteria; reason_col=:reject_reason, stage_col=:reject_stage)

For each row, record the label of the FIRST criterion it fails (and a short stage tag),
or "" if it passes all of them. Rows are NOT removed. Returns df.
"""
function annotate_rejections!(df, criteria::Vector{<:FilterCriterion};
                              reason_col=:reject_reason, stage_col=:reject_stage)
```

- Reuses the existing `passes(row, criterion)` dispatch.
- Unit test (mirror the `filters.jl` testset): first-failing label is recorded; passing rows
  get `""`.

### 1.2 `blast_discover` / `handle_blast`: annotate, then partition

In `src/discover/blast.jl`, change the staged `filter!`s
(edge → scov → pseudo → mismatch → core-cov → `full_count`/`full_ratio`/`length`/`aln_mismatch`)
so each stage **annotates** `reject_reason`/`reject_stage` on rows not already rejected,
instead of dropping. Carry all candidates to the end.

- A small helper `mark!(df, mask, reason, stage)` keeps the non-GermlineFilter stages uniform
  (it only marks rows that are still un-rejected).
- The final `GermlineFilter` block becomes `annotate_rejections!` with the same criteria.

### 1.3 Two outputs

- `<output>` (existing `always_gz` path) = rows with `reject_reason == ""`; identical schema
  to today.
- `<output>` with a `.full` infix (e.g. `KI_IMD_IGHV.full.tsv.gz`) or an explicit
  `--full-output PATH` = **all** candidates plus `reject_reason` and `reject_stage`, so the
  user can inspect every assignment and the reason it was rejected.
- Add `--full-output` (default: derive from `<output>`) to the blast arg table in
  `src/cmd/discover.jl`.

### 1.4 Staged reporting (fancy)

New `src/utils/report.jl` (`Report` module) reusing the colored style in
`Filters.report_filter_step`:

- `stage_report(name, kept, before; values=nothing)` — colored `▸ <stage>  kept K/B (−D)`;
  when the UnicodePlots extension is active (`Data.barplot_fn[]` set), also draw a histogram of
  `values` (mismatch, scov, corecov, full_count) via `Data.barplot_if_available`.
- Degrades to plain text when UnicodePlots is absent.
- **Binary note:** UnicodePlots is a `[weakdeps]` extension (kept out of the default load for
  speed). To get plots in the compiled binary, add `UnicodePlots` to `build/Project.toml` so
  `create_app` bundles it; the library stays optional.

### 1.5 Tests

- `annotate_rejections!` unit test.
- BLAST itself needs `blastn`; instead test the annotate/partition logic on a synthetic
  `clusters`-shaped DataFrame so the reason/stage columns and the full-vs-filtered split are
  verified without BLAST.

### Touchpoints

`src/utils/filters.jl`, `src/discover/blast.jl`, `src/cmd/discover.jl`,
`src/utils/report.jl` (new), `test/runtests.jl`, optionally `build/Project.toml`.

---

## Phase 2 — Discriminative metrics (true novel vs artifact)

Per-candidate cluster metrics, computed from the reads assigned to each candidate:

- `n_reads` (support),
- `consensus_support` — mean fraction of reads matching the cluster consensus,
- `mean_pairwise_hamming` within the cluster,
- **positional Shannon entropy** — averaged over aligned columns,
- `nearest_known_dist` — edit distance to the closest database allele,
- composition flags — GC content, homopolymer runs, N content.

Rationale: a *true* novel allele yields a tight cluster (high consensus support, low entropy)
at a meaningful edit distance from known alleles; artifacts are diffuse / low-support, or
trivially close to a known allele. Expose these as columns **and** optional `FilterCriterion`s,
so their effect is visible in `reject_reason` from Phase 1.

---

## Phase 3 — `discover selftest`

`discover selftest <reads.tsv> <base.fasta> <truth.fasta> [blast opts]`:

- Run discovery against **base** (so the novel alleles must be *discovered*).
- Define truth-novel = sequences in `truth` not in `base`.
- Match candidates to truth-novel (exact / small edit distance).
- Report **recall** (recovered novel / total novel) and **precision** (true novel / discovered
  novel), a per-allele table, and the FP / FN lists annotated with the Phase-2 metrics and the
  Phase-1 `reject_stage` — so you can see *which* filter killed a true allele and tune it.
- Optional UnicodePlots summary (recovered vs missed by edit-distance bin).
- Implemented as one more `Command` singleton in the registry
  (`struct DiscoverSelftest <: Command end`, `cli_path`, `run_command`, arg table in
  `src/cmd/discover.jl`).

This rides on Phases 1–2: the rejection reasons and metrics are exactly what make the
self-test actionable for tuning the blast thresholds.
