# ========================== Counting & metrics (stages) ==========================
"""
    add_counts!(result_df, sequence_lookup) -> df

Add `full_count` (identical full rows), `count` (per well/case/db_name/sequence), `gene`,
optional `isin_db`, and IgDiscover-style `full_allelic_ratio` / `allelic_ratio` (÷max).
"""
function add_counts!(result_df::DataFrame, sequence_lookup)
    df = transform(groupby(result_df, names(result_df)), nrow => :full_count)
    transform!(groupby(df, [:well, :case, :db_name, :sequence]), nrow => :count)
    transform!(df, :db_name => ByRow(x -> first(split(x, '*'))) => :gene)
    mark_novel!(df, optional(sequence_lookup))
    add_group_ratio!(df, :full_count, [:well, :case, :gene], FULL_ALLELIC_RATIO)
    add_group_ratio!(df, :count, [:well, :case, :gene], ALLELIC_RATIO)
    return df
end

mark_novel!(_, ::Absent) = nothing
function mark_novel!(df, lookup::Present)
    @info "Adding isin_db column based on reference FASTA"
    df[!, :isin_db] = map(row -> get(lookup.value, row.sequence, false) ? "" : "Novel", eachrow(df))
    return nothing
end

"""
    add_quality_metrics!(udf) -> udf

Add per-candidate-core metrics over distinct rows: `n_donors` (cross-donor recurrence),
`n_reads_total` (read support), `peak_allelic_ratio` (best per-donor full allelic ratio). Must be
called on the de-duplicated table so `full_count` is summed once per distinct row.
"""
function add_quality_metrics!(udf::DataFrame)
    transform!(groupby(udf, :sequence), :case => (x -> length(unique(x))) => :n_donors)
    transform!(groupby(udf, :sequence), :full_count => sum => :n_reads_total)
    transform!(groupby(udf, :sequence), FULL_ALLELIC_RATIO => maximum => PEAK_ALLELIC_RATIO)
    return udf
end

"""
    in_analysis_locus(db_name, locus) -> Bool

Whether `db_name` participates in locus-scoped frequency denominators. When `locus` is empty,
every allele is included (TCR, housekeeping, immunoglobulin, controls). Otherwise only names
starting with `locus` are included — use e.g. `IGHV` or `IG` to keep spike-in controls out of
IG/TCR frequency totals without dropping those rows from the output table.
"""
function in_analysis_locus(db_name::AbstractString, locus::AbstractString)
    prefix = strip(String(locus))
    return isempty(prefix) || startswith(String(db_name), prefix)
end

"""
    locus_group_stat!(df, groupcols, srccol, destcol, locus, statfn; default=0)

Within each group, set `destcol` to `statfn` over `srccol` for accepted rows (`reject_reason == ""`)
that pass `in_analysis_locus`; `default` when none qualify.
"""
function locus_group_stat!(df::DataFrame, groupcols, srccol::Symbol, destcol::Symbol, locus, statfn; default=0)
    transform!(groupby(df, groupcols)) do g
        fg = filter(r -> in_analysis_locus(r.db_name, locus) && isempty(r.reject_reason), g)
        DataFrame(destcol => fill(isempty(fg) ? default : statfn(fg[!, srccol]), nrow(g)))
    end
    return df
end

"""
    add_frequency_columns!(df, locus) -> df

Add the locus-frequency columns used by the reference-frequency filters: `gene_count` /
`case_count` (per-well/case totals), `allele_cohort_median` / `gene_cohort_median` (across-donor
medians), `gene_case_freq` (gene-usage fraction), `allele_cohort_fold` / `gene_cohort_fold`
(fold-change vs cohort median), and `gene_fraction` (count÷sum of accepted counts in gene).
All aggregates are over accepted (count/ratio-passing) rows. Requires `reject_reason`
(count/ratio annotated first).

Pass an empty `locus` to include every allele in the frequency denominators. Pass a prefix such
as `IGHV` or `TRBV` to scope denominators to that locus and zero out locus-scoped stats for other
alleles (e.g. spike-in controls) while still retaining their rows in the table.
"""
function add_frequency_columns!(df::DataFrame, locus::AbstractString)
    # Per-well/case totals (over accepted, in-locus rows).
    locus_group_stat!(df, [:well, :case, :gene], :count, :gene_count, locus, sum)
    locus_group_stat!(df, [:well, :case], :count, :case_count, locus, sum)
    # Cohort (across-donor) medians — robust central tendency used for the fold-change below.
    locus_group_stat!(df, [:db_name], :count, :allele_cohort_median, locus, median; default=0.0)
    locus_group_stat!(df, [:gene], :gene_count, :gene_cohort_median, locus, median; default=0.0)

    # gene_case_freq = gene-usage fraction in the case (low ⇒ possible deletion; see --deletion).
    df[:, :gene_case_freq] = safe_ratio.(df.gene_count, df.case_count)
    # *_cohort_fold = this donor's count ÷ the allele's/gene's cohort-median (robust fold-change
    # vs typical); a tiny fold flags a sporadic low-support observation.
    df[:, :allele_cohort_fold] = safe_ratio.(df.count, df.allele_cohort_median)
    df[:, :gene_cohort_fold] = safe_ratio.(df.gene_count, df.gene_cohort_median)
    # gene_fraction = allele reads ÷ all accepted reads in gene (÷sum, not IgDiscover allelic_ratio).
    transform!(groupby(df, [:well, :case, :gene])) do g
        denom = sum((r.count for r in eachrow(g) if isempty(r.reject_reason)); init=0)
        DataFrame(GENE_FRACTION => safe_ratio.(g.count, denom))
    end
    return df
end

# ========================== Filter assembly & annotation (stages) ==========================

"""
    exact_filter_criteria(; min_fullcount, min_count, min_allelic_ratio, min_full_allelic_ratio,
                          expect_dict, expect_full_dict, min_recurrence, min_seqlen,
                          min_peak_allelic_ratio, min_reads_total)

Count and allelic-ratio criteria for exact candidates. Each flag filters one output column;
set a flag to 0 to disable that filter. Per-gene overrides: `expect_dict` / `expect_full_dict`.
"""
function exact_filter_criteria(; min_fullcount::Int=0, min_count::Int=0,
                               min_allelic_ratio::Float64, min_full_allelic_ratio::Float64,
                               expect_dict, expect_full_dict,
                               min_recurrence::Int=0, min_seqlen::Int=0,
                               min_peak_allelic_ratio::Float64=0.0, min_reads_total::Int=0)
    crit = FilterCriterion[]
    min_count > 0 && push!(crit,
        MinThreshold(:count, Float64(min_count), "min count (--min-count $min_count)"))
    min_fullcount > 0 && push!(crit,
        MinThreshold(:full_count, Float64(min_fullcount),
                     "min full count (--min-fullcount $min_fullcount)"))
    min_seqlen > 0 && push!(crit,
        MinStringLength(:sequence, min_seqlen, "min sequence length (--min-seqlen $min_seqlen)"))
    if min_allelic_ratio > 0
        push!(crit, CustomFilter(r -> r.allelic_ratio >= get_ratio_threshold(expect_dict, r,
                                 type="allelic_ratio", default=min_allelic_ratio),
                     "min allelic ratio (--min-allelic-ratio $min_allelic_ratio)"))
    end
    if min_full_allelic_ratio > 0
        push!(crit, CustomFilter(r -> r.full_allelic_ratio >= get_ratio_threshold(expect_full_dict, r,
                                 type="full_allelic_ratio", default=min_full_allelic_ratio),
                     "min full allelic ratio (--min-full-allelic-ratio $min_full_allelic_ratio)"))
    end
    min_recurrence > 0 && push!(crit,
        MinThreshold(:n_donors, Float64(min_recurrence), "min donor recurrence (--min-recurrence $min_recurrence)"))
    min_peak_allelic_ratio > 0 && push!(crit,
        MinThreshold(PEAK_ALLELIC_RATIO, min_peak_allelic_ratio,
                     "min peak allelic ratio (--min-peak-allelic-ratio $min_peak_allelic_ratio)"))
    min_reads_total > 0 && push!(crit,
        MinThreshold(:n_reads_total, Float64(min_reads_total),
                     "min reads total (--min-reads-total $min_reads_total)"))
    return crit
end

"""
    exact_frequency_criteria(mod, expect_dict, deletion_dict, min_allele_fold, min_gene_fold)

The reference-frequency criteria: allele/gene-case frequency floors (control-gene-aware via
`get_ratio_threshold`) and the cross-case median ratios.
"""
function exact_frequency_criteria(mod, deletion_dict, min_gene_fraction,
                                min_gene_case_freq, min_allele_fold, min_gene_fold)
    crit = FilterCriterion[
        CustomFilter(x -> x.gene_case_freq >= mod.get_ratio_threshold(deletion_dict, x,
                         type="gene_case_freq", default=min_gene_case_freq),
                     "min gene case freq (--min-gene-case-freq $min_gene_case_freq)"),
        MinThreshold(:allele_cohort_fold, min_allele_fold, "min allele cohort fold (--min-allele-cohort-fold)"),
        MinThreshold(:gene_cohort_fold, min_gene_fold, "min gene cohort fold (--min-gene-cohort-fold)"),
    ]
    min_gene_fraction > 0 && pushfirst!(crit,
        CustomFilter(x -> x.gene_fraction >= min_gene_fraction,
                     "min gene fraction (--min-gene-fraction $min_gene_fraction)"))
    return crit
end

"""
    annotate_stage!(df, criteria, stage)

Apply each criterion in turn, marking (not dropping) the first rejection per row and printing
a per-criterion kept/removed line. Shared shape with the discovery pipelines.
"""
function annotate_stage!(df::DataFrame, criteria, stage::AbstractString)
    for criterion in criteria
        before = count(isempty, df.reject_reason)
        fail = Bool[!passes(row, criterion) for row in eachrow(df)]
        mark_rejected!(df, fail, criterion.label, stage)
        stage_report(criterion.label, count(isempty, df.reject_reason), before)
    end
    return df
end
