module Selftest
    # Evaluate how well `discover blast` recovers known-novel alleles. Inputs:
    #   - a discovery FULL table (every candidate + reject_reason / reject_stage),
    #   - a BASE reference FASTA (the reference used for discovery),
    #   - a TRUTH FASTA (known + novel).

    using DataFrames
    using CSV
    using Statistics: median
    using ..Blast: blast_discoverable_metrics, blast_cli_suggestion, build_blast_output_criteria
    using ..Cli: BLAST_DEFAULTS, BLAST_PRESETS
    using ..Data: load_fasta, barplot_if_available
    using ..Filters: passes, MetricView, LengthView, NumericView, metric_view
    using ..Option: Absent, Present, absent, optional
    using ..Report: section, stage_report

    export handle_selftest, evaluate_recovery, classify_allele, is_novel, metric_separation,
           recall_safe_filters, blast_discoverable_metrics, marginal_filter_shadowing, selftest_blast_block

    include("selftest_recovery.jl")
    include("selftest_metrics.jl")
    include("selftest_handle.jl")
end
