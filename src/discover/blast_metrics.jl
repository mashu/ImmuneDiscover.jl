"Summarize cross-donor recurrence of accepted candidates (single-donor cores ≈ likely artifacts)."
function report_recurrence(kept)
    nrow(kept) == 0 && return nothing
    uniq = unique(select(kept, [:aln_qseq, :n_donors]))
    n = nrow(uniq)
    single = count(==(1), uniq.n_donors)
    println("  $n distinct accepted candidate sequences; ",
            single, " seen in a single donor (", round(100 * single / n; digits=1),
            "% — more likely artifacts), ", n - single, " in two or more donors.")
    println("  donors per candidate (x = number of donors, bar height = candidates):")
    histogram_if_available(uniq.n_donors; nbins=20)
    return nothing
end

"""
    add_blast_support_columns!(df) -> df

After trimming (`aln_qseq` set): `count` = sum of `full_count` per (donor, allele, core);
then per-donor allelic ratios on `count` and `full_count` separately.
"""
function add_blast_support_columns!(df::DataFrame)
    transform!(groupby(df, [:well, :case, :sseqid, :aln_qseq]), :full_count => sum => :count)
    add_group_ratio!(df, :count, [:well, :case, :gene], ALLELIC_RATIO)
    add_group_ratio!(df, :full_count, [:well, :case, :gene], FULL_ALLELIC_RATIO)
    return df
end

"Post-output diagnostic columns (no `discover blast` threshold flag; self-test still scans them)."
const BLAST_DIAGNOSTIC_METRICS = [
    "nn_dist", "parent_ratio", "satellite_score", "chimera_score", "gc_content",
]

"Cluster / trim stages before output filters (rows failing these never reach the full table)."
const BLAST_UPSTREAM_METRICS = ["scov", "corecov", "blast_mismatch"]

"""
    build_blast_output_criteria(b; keep_failed, include_inactive) -> Vector{FilterCriterion}

Same criterion list as `handle_blast` output filters. `include_inactive=true` lists every
optional threshold column (for self-test metric discovery), not only flags > 0 in `b`.
"""
function build_blast_output_criteria(b::AbstractDict;
                                     keep_failed::Bool=true, include_inactive::Bool=false)
    min_count = get(b, "min-count", 0)
    min_fullcount = get(b, "min-fullcount", 0)
    min_allelic = get(b, "min-allelic-ratio", 0.0)
    min_full_allelic = get(b, "min-full-allelic-ratio", 0.0)
    min_peak_allelic = get(b, "min-peak-allelic-ratio", 0.0)
    min_length = get(b, "length", 0)
    min_recurrence = get(b, "min-recurrence", 0)
    max_homop = get(b, "max-homopolymer", 0)
    min_reads_total = get(b, "min-reads-total", 0)
    max_aln_mismatch = get(b, "max-aln-mismatch", 0)

    criteria = FilterCriterion[
        MinStringLength(:aln_qseq, min_length, "min trimmed length (--length $min_length)"),
    ]
    if !keep_failed
        push!(criteria, NonNegative(:core_aln_mismatch, "trimming failed (core_aln_mismatch < 0)"))
    end
    push!(criteria, MaxThreshold(:core_aln_mismatch, Float64(max_aln_mismatch),
                                 "max trimmed-core distance (--max-aln-mismatch $max_aln_mismatch)"))
    (include_inactive || min_count > 0) && push!(criteria,
        MinThreshold(:count, Float64(min_count), "min count (--min-count $min_count)"))
    (include_inactive || min_fullcount > 0) && push!(criteria,
        MinThreshold(:full_count, Float64(min_fullcount),
                     "min full count (--min-fullcount $min_fullcount)"))
    (include_inactive || min_allelic > 0) && push!(criteria,
        MinThreshold(ALLELIC_RATIO, min_allelic,
                     "min allelic ratio (--min-allelic-ratio $min_allelic)"))
    (include_inactive || min_full_allelic > 0) && push!(criteria,
        MinThreshold(FULL_ALLELIC_RATIO, min_full_allelic,
                     "min full allelic ratio (--min-full-allelic-ratio $min_full_allelic)"))
    (include_inactive || min_peak_allelic > 0) && push!(criteria,
        MinThreshold(PEAK_ALLELIC_RATIO, min_peak_allelic,
                     "min peak allelic ratio (--min-peak-allelic-ratio $min_peak_allelic)"))
    (include_inactive || min_reads_total > 0) && push!(criteria,
        MinThreshold(:n_reads_total, Float64(min_reads_total),
                     "min total reads (--min-reads-total $min_reads_total)"))
    (include_inactive || min_recurrence > 0) && push!(criteria,
        MinThreshold(:n_donors, Float64(min_recurrence),
                     "min donor recurrence (--min-recurrence $min_recurrence)"))
    (include_inactive || max_homop > 0) && push!(criteria,
        MaxThreshold(:max_homopolymer, Float64(max_homop),
                     "max homopolymer (--max-homopolymer $max_homop)"))
    return criteria
end

"""
    blast_discoverable_metrics(df) -> Vector{String}
    blast_discoverable_metrics(df, blast_block) -> Vector{String}

Every tunable / diagnostic metric column that `discover blast` can filter on or report,
intersected with columns present in `df`. Single source for self-test metric scans.
"""
blast_discoverable_metrics(df::DataFrame) = blast_discoverable_metrics(df, absent)
blast_discoverable_metrics(df::DataFrame, ::Absent) = blast_discoverable_metrics(df, Dict{String,Any}())
blast_discoverable_metrics(df::DataFrame, b::Present) = blast_discoverable_metrics(df, b.value)
function blast_discoverable_metrics(df::DataFrame, blast_block::AbstractDict)
    cols = Set(string.(names(df)))
    out = String[]
    seen = Set{String}()
    for c in build_blast_output_criteria(blast_block; include_inactive=true)
        m = criterion_column(c)
        isempty(m) || m in seen || (push!(seen, m); push!(out, m))
    end
    for m in vcat(BLAST_UPSTREAM_METRICS, BLAST_DIAGNOSTIC_METRICS, ["max_homopolymer"])
        m in cols && m in seen && continue
        m in cols && (push!(seen, m); push!(out, m))
    end
    return out
end

"""
    blast_cli_suggestion(metric, direction, threshold) -> String

Map a discovery-table column to the matching `discover blast` CLI flag when one exists.
"""
function blast_cli_suggestion(metric::AbstractString, direction::AbstractString, threshold::Real)
    t4 = round(Float64(threshold); digits=4)
    ti = round(Int, threshold)
    metric == "peak_allelic_ratio" && direction == "keep ≥" && return "--min-peak-allelic-ratio $t4"
    metric == "full_allelic_ratio" && direction == "keep ≥" && return "--min-full-allelic-ratio $t4"
    metric == "allelic_ratio" && direction == "keep ≥" && return "--min-allelic-ratio $t4"
    metric == "count" && direction == "keep ≥" && return "--min-count $ti"
    metric == "full_count" && direction == "keep ≥" && return "--min-fullcount $ti"
    metric == "aln_qseq" && direction == "keep ≥" && return "--length $ti"
    metric == "corecov" && direction == "keep ≥" && return "--min-corecov $t4"
    metric == "scov" && direction == "keep ≥" && return "--subjectcov $t4"
    metric == "n_reads_total" && direction == "keep ≥" && return "--min-reads-total $ti"
    metric == "n_donors" && direction == "keep ≥" && return "--min-recurrence $ti"
    metric == "blast_mismatch" && direction == "keep ≤" && return "--max-blast-mismatch $ti"
    metric == "core_aln_mismatch" && direction == "keep ≤" && return "--max-aln-mismatch $ti"
    metric == "max_homopolymer" && direction == "keep ≤" && return "--max-homopolymer $ti"
    return "$metric $direction $t4"
end
