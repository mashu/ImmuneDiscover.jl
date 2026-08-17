# ========================== CLI handler (orchestrator) ==========================

function handle_exact(parsed_args, immunediscover_module, always_gz)
    @info "Exact search"
    ex = parsed_args["search"]["exact"]
    extension = ex["extension"]
    border = get(ex, "border", 0)
    adjust_per_gene_extension = get(ex, "adjust-per-gene-extension", false)
    adjust_percent = get(ex, "adjust-percent", 1.0)
    limit = ex["limit"]
    refgenes = ex["refgene"]
    length(refgenes) > 0 && @info "Using reference genes $refgenes"

    table = immunediscover_module.load_demultiplex(ex["tsv"])
    limit > 0 && (@info "Limiting reads to $limit"; table = table[1:limit, :])
    db = immunediscover_module.load_fasta(ex["fasta"])
    min_fullcount = ex["min-fullcount"]
    min_allelic_ratio = ex["min-allelic-ratio"]
    min_fullcount > 0 && min_fullcount < 5 && @warn "Decreasing --min-fullcount below 5 may lead to false positives"
    top = ex["top"]
    affix = ex["affix"]
    gene = ex["gene"]
    locus = ex["locus"]

    # Foot-gun guard: -g sets the RSS-extraction orientation (V/J extract one heptamer, D both
    # sides). Warn early if the reference looks like a different gene than -g.
        ref_types = GeneType[]
        for (n, _) in db
            t = gene_type_from_name(string(n))
            t === nothing || push!(ref_types, t)
        end
        if !isempty(ref_types)
            majority = gene_string(majority_gene(ref_types))
            majority != gene && @warn "You passed -g $gene but the reference FASTA looks like $majority genes — RSS is extracted in $gene orientation (only the 3' heptamer for V, 5' for J; both sides for D). Pass -g $majority to extract the correct RSS."
        end

    local rss
    if extension !== nothing
        @info "Using extension mode with length $extension"; rss = String[]
    else
        rss = split(ex["rss"], ',')
        immunediscover_module.validate_types(rss)
        @info "Extract RSS: $(join(rss,','))"
    end
    top != 1 && @info "Uncollapsed mode; at most $top full records returned."

    # `expect`/`deletion` control-gene threshold files serve two distinct, name-keyed roles:
    # expect_dict: per-gene/per-allele overrides for --min-allelic-ratio (IgDiscover allele_ratio).
    # deletion_dict: per-gene gene_case_freq overrides only.
    expect_dict = load_ratio_dict(ex["expect"])
    expect_full_dict = load_ratio_dict(ex["expect-full"])
    deletion_dict = load_ratio_dict(ex["deletion"])

    raw = ex["raw"]
    sequence_lookup = ex["ref-fasta"] !== nothing ? build_sequence_lookup(ex["ref-fasta"]) : nothing

    counts_df = exact_search(table, db, gene; affix=affix, rss=rss, extension=extension, N=top,
        raw=raw, sequence_lookup=sequence_lookup, border=border,
        adjust_per_gene_extension=adjust_per_gene_extension, adjust_percent=adjust_percent)
    if nrow(counts_df) == 0
        @warn "No exact matches"
        return
    end
    add_chimera_scores!(counts_df, refs_by_gene(db); seq_col=:sequence, gene_col=:gene)
    sort!(counts_df, [:case, :db_name])

    # Count/ratio filters — annotate (don't drop) so the full table records every candidate.
    section("Exact search — count and ratio filters")
    init_rejection_columns!(counts_df)
    annotate_stage!(counts_df,
        exact_filter_criteria(; min_fullcount=min_fullcount,
            min_count=ex["min-count"],
            min_allelic_ratio=min_allelic_ratio,
            min_full_allelic_ratio=ex["min-full-allelic-ratio"],
            expect_dict=expect_dict, expect_full_dict=expect_full_dict,
            min_recurrence=get(ex, "min-recurrence", 0),
            min_seqlen=get(ex, "min-seqlen", 0),
            min_peak_allelic_ratio=get(ex, "min-peak-allelic-ratio", 0.0),
            min_reads_total=get(ex, "min-reads-total", 0)),
        "count and ratio filter")

    if !ex["noplot"]
        plotdf = accepted(counts_df)
        n_donors_input = length(unique(table.case))
        nrow(plotdf) > 0 ? immunediscover_module.plotgenes(plotdf; n_donors_in_input=n_donors_input) :
            @warn "No exact matches to plot"
    end

    if isempty(strip(locus))
        @info "Computing frequency columns across all alleles (no --locus prefix filter)"
    else
        @info "Scoping frequency denominators to db_name prefixes starting with $locus"
    end
    add_frequency_columns!(counts_df, locus)

    section("Exact search — frequency filters")
    annotate_stage!(counts_df,
        exact_frequency_criteria(immunediscover_module, deletion_dict,
            ex["min-gene-fraction"], ex["min-gene-case-freq"],
            ex["min-allele-cohort-fold"], ex["min-gene-cohort-fold"]),
        "frequency filter")

    reason_cols = [:reject_reason, :reject_stage]
    kept = accepted(counts_df)
    if length(refgenes) > 0
        for refgene in refgenes
            kept = grouped_ratios(kept, refgene, count_col=:count)
            transform!(groupby(kept, [:well, :case, :gene]), :count => sum => :ref_gene_count)
            kept = grouped_ratios(kept, refgene, count_col=:ref_gene_count)
        end
    end
    output = always_gz(ex["output"])
    full_output = always_gz(replace(replace(output, r"\.gz$" => ""), r"\.tsv$" => "") * ".full.tsv")
    gt = parse_gene_type(gene)
    report_exact_findings(counts_df, kept, db, gt, extension, table)
    # Readable output: metrics left, long sequence/flank columns right (genomic order); round
    # float columns to 4 dp instead of full Float64 precision.
    kept = round_floats!(order_exact_columns(kept, gt, extension))
    counts_df = round_floats!(order_exact_columns(counts_df, gt, extension))
    CSV.write(output, select(kept, Not(reason_cols)), compress=true, delim='\t')
    printstyled("  ✓ "; color=:green, bold=true); println("filtered → $output  ($(nrow(kept)) rows)")
    CSV.write(full_output, counts_df, compress=true, delim='\t')
    printstyled("  ✓ "; color=:green, bold=true); println("full     → $full_output  ($(nrow(counts_df)) candidates + reject_reason)")
    return
end
