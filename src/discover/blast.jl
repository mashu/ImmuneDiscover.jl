module Blast
    using CSV
    using CodecZlib
    using DataFrames
    using FASTX
    using DataStructures
    using BioAlignments
    using BioSequences
    using Folds
    using MD5
    using ProgressMeter: Progress, next!
    using Base.Threads: nthreads

    using ..Align
    using ..Mosaic: refs_by_gene, add_chimera_scores!
    using ..Data: load_fasta as data_load_fasta, unique_name, histogram_if_available
    using ..SeqStats: gc_content, max_homopolymer
    using ..Spans: each_exact_span
    using ..RatioColumns: ALLELIC_RATIO, FULL_ALLELIC_RATIO, PEAK_ALLELIC_RATIO
    using ..Filters: FilterCriterion, MinThreshold, MaxThreshold, MinStringLength, NonNegative,
                     add_group_ratio!, init_rejection_columns!, mark_rejected!, accepted, passes,
                     criterion_column
    using ..Report: stage_report, section, cluster_profile_heatmap, params_report, report_rejections

    export blast_discover, save_to_fasta, accumulate_affixes, save_extended, handle_blast
    export resolve_work_dir, blast_hits_gz_path, blast_cache_key
    export build_blast_output_criteria, blast_discoverable_metrics, blast_cli_suggestion
    export consensus_prefix, consensus_suffix, name_candidate

    include("blast_io.jl")
    include("blast_process.jl")
    include("blast_affix.jl")
    include("blast_neighbors.jl")
    include("blast_discover.jl")
    include("blast_output.jl")
end
