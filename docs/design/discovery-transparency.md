# Discovery transparency, metrics, and self-test

The discovery and search commands (`discover blast`, `search exact`, `discover hsmm`) explain
*why* every candidate is kept or rejected, expose metrics that separate genuine novel alleles
from artifacts, and let you measure recovery of known-vs-novel alleles. This page describes how
those pieces fit together; the richest set (quality metrics) lives in `discover blast`.

## Transparency: reject reasons and two output tables

`discover blast`, `search exact`, and `discover hsmm` all annotate their filters rather than
dropping rows, and each writes two tables: the filtered result, and a full annotated table
(`<output>.full.tsv.gz`) with `reject_reason` / `reject_stage` for every candidate.

Filtering annotates rather than drops. Each candidate carries `reject_reason` (the label of the
first filter it failed, empty if accepted) and `reject_stage` (the stage that rejected it). The
run therefore produces two tables:

- `<output>` — the filtered discoveries (candidates with an empty `reject_reason`).
- `<output>.full.tsv.gz` (or `--full-output PATH`) — every candidate plus `reject_reason` /
  `reject_stage`, so you can inspect each assignment and the reason it was dropped.

`search exact` writes a **slim** TSV by default (identifiers, counts, allelic ratios, flanks).
Pass `--diagnostic` to include intermediate statistics (cohort folds, `chimera_score`,
denominators). Filters still use those values when the columns are omitted. `--raw PATH` is a
separate uncollapsed match dump from before filtering.

Each stage prints a colored summary line (kept/before, removed), and the per-criterion output
filters are reported separately, so it is clear how many candidates each filter removes.
`subject coverage`, `BLAST mismatch`, and `core coverage` also draw a unicode histogram of their
metric.

Every discovery/search command closes with a findings report (shared `Report` helpers):

- `report_rejections` — a bar plot of candidates per outcome (accepted bucket + one bar per reject
  reason), so it's obvious which filter removed what.
- `recurrence_report` — cross-donor recurrence of the accepted candidates (single-donor share ≈
  artifacts) as a histogram.
- `filter_quality_report` — accepted-vs-rejected medians of metrics the filters do *not* key on
  (e.g. `n_donors`): higher-for-accepted is independent evidence the filters keep the real alleles.
- `rss_consistency` (`search exact`, RSS mode) — the heptamer **consensus** (e.g. `CACAGTG`), mean
  conservation, and a per-position **variation** bar plot, so a clean RSS motif and where it varies
  are obvious (replaces the unlabelled composition heatmap).

`search exact` additionally reports: accepted alleles per gene; **per-donor QC** — genes/alleles
per donor and read depth as box plots, plus the weakest donors by name (to spot failed donors);
**per-gene amplification** — the read-count distribution per gene as a box plot, genes sorted by
median (best-amplifying first); and which reference genes never matched a read versus matched but
were fully filtered out.

The annotate path lives in `Filters` (`annotate_rejections!`, `mark_rejected!`,
`init_rejection_columns!`, `accepted`) and is shared with the rest of the pipeline; the colored
diagnostics live in `Report`.

## Discriminative metrics (true novel vs artifact)

Every candidate gets quality columns (in both output tables). A genuine novel allele tends to
recur across donors with solid read support and benign composition; an artifact is usually seen
in one donor, has few reads, sits a single base from a much more abundant allele, or is
composition-extreme.

- `n_donors` — distinct donors sharing the exact trimmed core (cross-donor recurrence; the
  strongest single signal).
- `n_reads_total` — reads backing the core across the run.
- `gc_content`, `max_homopolymer` — composition of the trimmed core.
- `nn_dist` / `parent_ratio` — between-cluster separation: edit distance to the nearest
  *more-abundant* core of the same gene and the parent's read ratio. A small `nn_dist` with a
  large `parent_ratio` marks an error satellite of a dominant allele.

`SeqStats` provides the reusable building blocks (`gc_content`, `max_homopolymer`, `n_content`,
`shannon_entropy`, `consensus_fraction`); the run summary shows a donor-recurrence histogram and
a base-composition heatmap (with mean positional entropy) of the accepted candidates.

Two of these are also available as optional filters (off by default): `--min-recurrence` requires
`n_donors ≥ N`, and `--max-homopolymer` drops cores with a longer homopolymer run. When set, they
annotate `reject_reason` like any other filter.

## Self-test (`discover selftest`)

`discover selftest <discovery.full.tsv.gz> <base.fasta> <truth.fasta> <report.tsv>` scores how
well discovery recovers known-novel alleles. Run `discover blast` against the BASE reference
(known alleles only) to produce the full annotated table, then:

- truth-novel = sequences in TRUTH not in BASE — the alleles that must be discovered;
- each truth-novel allele is classified **recovered** (an accepted candidate matches it),
  **rejected** (a candidate matched it but was filtered — annotated with the `reject_stage` that
  dropped it), or **missed** (never seen as a candidate);
- it reports **recall** (recovered / truth-novel), **precision** (true / accepted-novel cores),
  an outcome bar plot, and a rejected-by-stage bar plot showing which filter to relax.

### Metric separation (threshold tuning)

The self-test also labels every *novel* candidate row (core absent from BASE) as **true** (its
core matches a truth-novel allele) or **false**, then, for each metric column in the full table
(`n_donors`, `n_reads_total`, `scov`, `corecov`, mismatch, `gc_content`, `max_homopolymer`,
`nn_dist`, `parent_ratio`, …), finds the single threshold — and keep-direction (`≥` or `≤`) —
that best separates true from false by Youden's J (TP-rate − FP-rate). The metrics are printed
sorted by separation, so the top row is the most discriminative filter and its suggested cut. This
turns the self-test from "did we recover allele X" into "which blast threshold to set, and to
what value", and is saved with `--metrics-output`.

Because it consumes the full table, the self-test ties the reject reasons and quality metrics
together: it shows exactly which filter killed a true allele, and which metric best separates true
from false candidates — which is what makes the blast thresholds tunable.
