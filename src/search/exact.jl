module Exact
    using CSV
    using DataFrames
    using ProgressMeter
    using Folds
    using FASTX
    using Statistics
    using ..Gene: GeneType, VGene, DGene, JGene, parse_gene_type, gene_type_from_name,
                  gene_string, majority_gene
    using ..Spans: each_exact_span, flank_slice
    using ..RatioColumns: ALLELIC_RATIO, FULL_ALLELIC_RATIO, PEAK_ALLELIC_RATIO, GENE_FRACTION
    using ..Filters: FilterCriterion, MinThreshold, MinStringLength, CustomFilter, add_group_ratio!,
                     init_rejection_columns!, mark_rejected!, accepted, passes
    using ..Mosaic: refs_by_gene, add_chimera_scores!
    using ..Data: barplot_if_available, boxplot_if_available, round_floats!, load_fasta, get_ratio_threshold
    using ..Report: section, stage_report, report_rejections,
                    filter_quality_report, rss_consistency

    include("exact_flanks.jl")
    include("exact_match.jl")
    include("exact_metrics.jl")
    include("exact_report.jl")
    include("exact_search.jl")
    include("exact_handle.jl")

    export grouped_ratios, transform_counts, build_sequence_lookup, handle_exact, exact_search
    export in_analysis_locus
end
