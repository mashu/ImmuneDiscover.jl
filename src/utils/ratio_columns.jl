module RatioColumns
    """
    Shared ratio column names (match CLI flags: hyphens → underscores).

    | Column | Formula | CLI flag (0 = filter off) |
    |--------|---------|---------------------------|
    | `allelic_ratio` | `count` ÷ max(`count`) in donor+gene | `--min-allelic-ratio` |
    | `full_allelic_ratio` | `full_count` ÷ max(`full_count`) in donor+gene | `--min-full-allelic-ratio` |
    | `peak_allelic_ratio` | max `full_allelic_ratio` across donors | `--min-peak-allelic-ratio` |
    | `gene_fraction` | `count` ÷ sum(accepted `count`) in donor+gene | `--min-gene-fraction` |

    Count floors: `--min-count` (`count`), `--min-fullcount` (`full_count`). 0 = off.
    Per-gene TSV overrides: `--expect` (`allelic_ratio`), `--expect-full` (`full_allelic_ratio`).
    """

    const ALLELIC_RATIO = :allelic_ratio
    const FULL_ALLELIC_RATIO = :full_allelic_ratio
    const PEAK_ALLELIC_RATIO = :peak_allelic_ratio
    const GENE_FRACTION = :gene_fraction

    export ALLELIC_RATIO, FULL_ALLELIC_RATIO, PEAK_ALLELIC_RATIO, GENE_FRACTION
end
