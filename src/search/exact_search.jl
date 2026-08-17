# ========================== exact_search (orchestrator) ==========================

"""
    exact_search(table, query, gene; kwargs...) -> DataFrame

Find exact occurrences of each `query` allele in the reads and return the UNFILTERED candidate
table (one row per distinct flank/sequence, capped at `N` flank records per allele) with counts,
allelic ratios and quality metrics. Callers annotate/filter as needed (`handle_exact` runs the
full transparency cascade; `hsmm` keeps `full_count ≥ select_min_fullcount`).
"""
function exact_search(table, query, gt::GeneType; affix=13, rss=["heptamer", "spacer", "nonamer"],
                      extension=absent, N=10, raw=absent, sequence_lookup=absent,
                      border::Int=0, adjust_per_gene_extension::Bool=false, adjust_percent::Float64=1.0)
    @assert all([name in names(table) for name in ["well","case","name","genomic_sequence"]]) "File must contain following columns: well, case, name, genomic_sequence"

    ext = optional(extension)
    per_gene_prefix, per_gene_suffix = calibrated_extensions(ext, table, query, gt, border,
                                                            adjust_per_gene_extension, adjust_percent)
    result_df, totals_all, accepted_all = collect_matches(table, query, gt, affix, rss, ext,
        border, adjust_per_gene_extension, per_gene_prefix, per_gene_suffix)
    summarize_border_stats!(totals_all, accepted_all, border_filter_active(ext, border))

    isempty(result_df) && return result_df
    write_raw_matches!(optional(raw), result_df)

    df = add_counts!(result_df, optional(sequence_lookup))
    sort!(df, [:full_count, :count], rev=[true, true])
    udf = sort(unique(df), [:well, :case, :gene, :db_name, :sequence])
    add_quality_metrics!(udf)

    priority_columns = ["well", "case", "gene", "db_name", "count", "full_count",
                        "allelic_ratio", "full_allelic_ratio"]
    remaining_columns = setdiff(names(udf), priority_columns)
    udf = udf[:, vcat(priority_columns, remaining_columns)]
    gdf = groupby(udf, [:well, :case, :gene, :db_name, :sequence])
    udf_indexed = transform(gdf, :well => (x -> 1:length(x)) => :flank_index)
    return filter(x -> x.flank_index <= N, udf_indexed)
end

exact_search(table, query, gene::AbstractString; kwargs...) =
    exact_search(table, query, parse_gene_type(gene); kwargs...)

calibrated_extensions(::Absent, _, _, _, _, _, _) = (Dict{String,Int}(), Dict{String,Int}())
function calibrated_extensions(e::Present, table, query, gt, border, adjust, pct)
    (border > 0 && adjust) || return (Dict{String,Int}(), Dict{String,Int}())
    return calibrate_extension(table, query, gt, e.value, border, pct)
end

border_filter_active(::Absent, _) = false
border_filter_active(::Present, border::Int) = border > 0

write_raw_matches!(::Absent, _) = nothing
write_raw_matches!(p::Present, df) = CSV.write(p.value * ".gz", df, delim='\t', compress=true)
