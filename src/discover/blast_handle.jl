# CLI orchestrator for `discover blast`. Extension vs no-extension is a type, not a flag.

struct NoExtension end
struct AffixExtension
    forward::Int
    reverse::Int
end

extension_plan(forward::Int, reverse::Int) =
    (forward == 0 && reverse == 0) ? NoExtension() : AffixExtension(forward, reverse)

function warn_short_extension(n::Int, label::AbstractString)
    (n < 7) && (n > 0) && @warn "$label extension $n is short and may lead to false positives"
    return nothing
end

function combined_reference(fasta_path, pseudo, work_dir, fasta_stem)
    db_p = Vector{Tuple{String, String}}()
    if !isempty(pseudo)
        for (name, seq) in data_load_fasta(pseudo)
            push!(db_p, ("P" * name, seq))
        end
    end
    for (name, seq) in data_load_fasta(fasta_path)
        push!(db_p, (name, seq))
    end
    combined_fasta_path = joinpath(work_dir, fasta_stem * "-combined.fasta")
    save_to_fasta(db_p, combined_fasta_path)
    return db_p, combined_fasta_path
end

function prepare_extended_fasta(::NoExtension, combined_fasta_path, affixes_path, _, _)
    @info "No sequence extension requested, using original sequences"
    empty_affixes = [(name="", prefix="", suffix="")]
    CSV.write(affixes_path, DataFrame(empty_affixes), delim='\t')
    return combined_fasta_path
end

function prepare_extended_fasta(ext::AffixExtension, combined_fasta_path, affixes_path, b, work_dir)
    fasta_stem = file_stem(b["fasta"])
    ext_fasta_path = joinpath(work_dir, fasta_stem * "-combined-extended.fasta")
    if isfile(ext_fasta_path) && !b["overwrite"]
        @info "Using existing extended sequences from $ext_fasta_path"
        return ext_fasta_path
    end
    demux = load_csv(b["input"])
    db_p = data_load_fasta(combined_fasta_path)
    @info "Extending gene sequences by $(ext.forward) forward and $(ext.reverse) reverse nucleotides" reads=nrow(demux) references=length(db_p) julia_threads=nthreads()
    extended = accumulate_affixes(db_p, demux,
        forward_extension=ext.forward, reverse_extension=ext.reverse)
    affixes = save_extended(extended, ext_fasta_path)
    CSV.write(affixes_path, DataFrame(affixes, [:name, :prefix, :suffix]), delim='\t')
    @info "Saved affixes in $affixes_path"
    return ext_fasta_path
end

function load_affix_dict(affixes_path)
    affix_dict = Dict{String, Tuple{String, String}}()
    for row in eachrow(CSV.File(affixes_path, delim='\t') |> DataFrame)
        name = strip(ismissing(row.name) ? "" : String(row.name))
        prefix = ismissing(row.prefix) ? "" : String(row.prefix)
        suffix = ismissing(row.suffix) ? "" : String(row.suffix)
        isempty(name) && continue
        affix_dict[name] = (prefix, suffix)
    end
    return affix_dict
end

assign_cores!(::NoExtension, blast_clusters, _, _, _, _) = begin
    blast_clusters[:, :aln_qseq] = blast_clusters[:, :qseq]
    blast_clusters[:, :core_aln_mismatch] = blast_clusters[:, :blast_mismatch]
    return blast_clusters
end

function assign_cores!(ext::AffixExtension, blast_clusters, affixes_path, DB, minquality, verbose)
    isfile(affixes_path) || return blast_clusters
    @info "Loading affixes from $affixes_path"
    affix_dict = load_affix_dict(affixes_path)
    stats = AlignmentStats()
    check_affix_quality_warning(ext.forward, minquality)
    check_affix_quality_warning(ext.reverse, minquality)

    db_keys = Set(strip(String(name)) for (name, seq) in DB)
    db_dict = Dict(strip(String(name)) => String(seq) for (name, seq) in DB)
    db_lengths = Dict(strip(String(name)) => length(seq) for (name, seq) in DB)

    results = Folds.map(eachrow(blast_clusters)) do row
        local_stats = AlignmentStats()
        sseqid = strip(String(row.sseqid))
        base_key = sseqid_to_db_key(sseqid, db_keys)
        qseq = String(row.qseq)
        prefix, suffix = get(affix_dict, base_key, ("", ""))
        if !haskey(db_dict, base_key)
            push!(local_stats.failed_genes[sseqid], "no unextended reference for $base_key")
            return ("", -1, local_stats)
        end
        reference = db_dict[base_key]
        trimmed, distance = trim_and_align_sequence(qseq, prefix, suffix, reference, local_stats,
            min_quality=minquality, sseqid=sseqid)
        return (trimmed, distance, local_stats)
    end
    blast_clusters[:, :aln_qseq] = [r[1] for r in results]
    blast_clusters[:, :core_aln_mismatch] = [r[2] for r in results]
    for r in results
        merge_stats!(stats, r[3])
    end

    blast_clusters[:, :db_length] = [get(db_lengths, sseqid_to_db_key(strip(String(row.sseqid)), db_keys), 0) for row in eachrow(blast_clusters)]
    blast_clusters[:, :corecov] = map(eachrow(blast_clusters)) do row
        row.db_length > 0 ? length(row.aln_qseq) / row.db_length : 0.0
    end
    @info "Alignment stats: $(stats.total_attempts) attempts, $(stats.prefix_failures) prefix, $(stats.suffix_failures) suffix failures"
    if verbose
        for (gene_name, messages) in stats.failed_genes
            for msg in messages
                @info "  Failed [$gene_name]: $msg"
            end
        end
    end
    return blast_clusters
end

function mark_corecov!(blast_clusters, min_corecov)
    hasproperty(blast_clusters, :corecov) || return blast_clusters
    before_corecov = nrow(blast_clusters)
    section("BLAST discovery — trimming & core coverage")
    mark_rejected!(blast_clusters, blast_clusters.corecov .< min_corecov,
                   "core coverage < $min_corecov (core÷DB len)", "core coverage")
    kept_corecov = count(isempty, blast_clusters.reject_reason)
    stage_report("core coverage ≥ $min_corecov (core÷DB len)", kept_corecov, before_corecov;
                 values=blast_clusters.corecov, histogram=true, hist_closed=:right)
    return blast_clusters
end

function add_quality_metrics!(blast_clusters)
    blast_clusters[:, :gc_content] = gc_content.(blast_clusters.aln_qseq)
    blast_clusters[:, :max_homopolymer] = max_homopolymer.(blast_clusters.aln_qseq)
    transform!(groupby(blast_clusters, :aln_qseq), :case => (x -> length(unique(x))) => :n_donors)
    transform!(groupby(blast_clusters, :aln_qseq), :full_count => sum => :n_reads_total)
    transform!(groupby(blast_clusters, :aln_qseq), FULL_ALLELIC_RATIO => maximum => PEAK_ALLELIC_RATIO)
    return blast_clusters
end

function apply_output_filters!(blast_clusters, criteria)
    section("BLAST discovery — output filters")
    init_rejection_columns!(blast_clusters)
    for criterion in criteria
        before = count(isempty, blast_clusters.reject_reason)
        fail = Bool[!passes(row, criterion) for row in eachrow(blast_clusters)]
        mark_rejected!(blast_clusters, fail, criterion.label, "output filter")
        stage_report(criterion.label, count(isempty, blast_clusters.reject_reason), before)
    end
    return blast_clusters
end

function name_and_score!(blast_clusters, DB, isin)
    db_seqs = [(strip(String(n)), String(s)) for (n, s) in DB]
    db_by_gene = refs_by_gene(db_seqs)
    exact_lookup = seq_to_name_lookup(db_seqs)
    n_accepted = count(isempty, blast_clusters.reject_reason)
    @info "Post-filter analysis: $n_accepted accepted / $(nrow(blast_clusters)) candidates"
    @info "Computing chimera scores (accepted candidates only)"
    add_chimera_scores!(blast_clusters, db_by_gene; seq_col=:aln_qseq, gene_col=:gene,
                        only_accepted=true)
    @info "Naming candidates"
    blast_clusters[:, :allele_name] = map(r -> name_candidate(String(r.aln_qseq), r.sseqid,
                                                              r.core_aln_mismatch, exact_lookup,
                                                              db_seqs, isin),
                                          eachrow(blast_clusters))
    @info "Computing neighbor statistics"
    add_neighbor_stats!(blast_clusters)
    return db_seqs
end

function blast_output_paths(b, always_gz)
    output = always_gz(b["output"])
    full = optional(get(b, "full-output", nothing))
    return output, full_table_path(full, output, always_gz)
end

full_table_path(::Absent, output, always_gz) =
    always_gz(replace(replace(output, r"\.gz$" => ""), r"\.tsv$" => "") * ".full.tsv")
full_table_path(p::Present, _, always_gz) = always_gz(p.value)

function write_blast_tables(blast_clusters, output, full_output, db_seqs)
    reason_cols = [:reject_reason, :reject_stage]
    section("BLAST discovery — summary")
    kept = accepted(blast_clusters)
    stage_report("accepted (passed all filters)", nrow(kept), nrow(blast_clusters))
    report_rejections(blast_clusters.reject_reason)
    report_recurrence(kept)
    n_close = count((kept.nn_dist .>= 0) .& (kept.nn_dist .<= 1))
    n_shadow = count(kept.likely_satellite)
    n_close > 0 && @info "$n_close accepted candidate(s) are within 1 bp of a more-abundant core; $n_shadow flagged likely_satellite (satellite_score ≥ 0.5). See nn_dist / parent_ratio / satellite_score / chimera_score — not a default filter."
    hm_kept = filter(r -> !r.likely_satellite, kept)
    if nrow(hm_kept) >= 1
        uc = combine(groupby(hm_kept, :aln_qseq),
                     :gene => first => :gene,
                     :allele_name => first => :allele_name,
                     :sseqid => first => :sseqid,
                     :core_aln_mismatch => first => :core_aln_mismatch,
                     :n_reads_total => first => :n_reads_total)
        cluster_profile_heatmap(String.(uc.gene), String.(uc.aln_qseq),
                                String.(uc.allele_name), String.(uc.sseqid);
                                reads=uc.n_reads_total,
                                core_aln_mismatch=uc.core_aln_mismatch,
                                db_seqs=db_seqs,
                                title="novel alleles")
    end
    CSV.write(output, select(kept, Not(reason_cols)), compress=true, delim='\t')
    printstyled("  ✓ "; color=:green, bold=true); println("filtered  → $output  ($(nrow(kept)) rows)")
    CSV.write(full_output, blast_clusters, compress=true, delim='\t')
    printstyled("  ✓ "; color=:green, bold=true); println("full      → $full_output  ($(nrow(blast_clusters)) candidates + reject_reason)")
end

function handle_blast(parsed_args, immunediscover_module, always_gz)
    using_cli = immunediscover_module.Cli
    if get(parsed_args["discover"]["blast"], "show-presets", false)
        using_cli.show_blast_presets()
        return
    end
    parsed_args = using_cli.apply_blast_presets!(parsed_args)
    gene = parsed_args["discover"]["blast"]["gene"]
    haskey(using_cli.BLAST_PRESETS, gene) && @info "Applied $gene gene preset (overrides logged above)"
    params_report(parsed_args["discover"]["blast"], using_cli.BLAST_PARAM_GROUPS;
                  title="discover blast — parameters")
    @info "Discovery with BLAST assignments"

    b = parsed_args["discover"]["blast"]
    fasta_path = b["fasta"]
    work_dir = resolve_work_dir(b["work-dir"])
    isdir(work_dir) || mkpath(work_dir)
    @info "Work directory (caches, intermediates): $work_dir"
    fasta_stem = file_stem(fasta_path)
    affixes_path = joinpath(work_dir, fasta_stem * ".affixes")
    DB = immunediscover_module.load_fasta(fasta_path)

    verbose = b["verbose"]
    minquality = b["minquality"]
    ext = extension_plan(b["forward"], b["reverse"])
    keep_failed = b["keep-failed"]
    min_corecov = b["min-corecov"]
    isin = b["isin"]
    warn_short_extension(b["forward"], "Forward")
    warn_short_extension(b["reverse"], "Reverse")

    _, combined_fasta_path = combined_reference(fasta_path, b["pseudo"], work_dir, fasta_stem)
    ext_fasta_path = prepare_extended_fasta(ext, combined_fasta_path, affixes_path, b, work_dir)

    blast_clusters = blast_discover(
        b["input"],
        ext_fasta_path;
        work_dir=work_dir,
        max_dist=b["max-blast-mismatch"],
        min_edge=b["edge"],
        min_scov=b["subjectcov"],
        args=b["args"],
        verbose=verbose,
        overwrite=b["overwrite"],
        min_read_length=get(b, "min-read-length", 0),
    )

    assign_cores!(ext, blast_clusters, affixes_path, DB, minquality, verbose)
    mark_corecov!(blast_clusters, min_corecov)
    add_blast_support_columns!(blast_clusters)
    add_quality_metrics!(blast_clusters)
    apply_output_filters!(blast_clusters, build_blast_output_criteria(b; keep_failed=keep_failed))
    db_seqs = name_and_score!(blast_clusters, DB, isin)
    output, full_output = blast_output_paths(b, always_gz)
    write_blast_tables(blast_clusters, output, full_output, db_seqs)
end
