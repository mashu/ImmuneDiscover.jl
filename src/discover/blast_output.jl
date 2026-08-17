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
    blast_discoverable_metrics(df; blast_block) -> Vector{String}

Every tunable / diagnostic metric column that `discover blast` can filter on or report,
intersected with columns present in `df`. Single source for self-test metric scans.
"""
blast_discoverable_metrics(df::DataFrame) = blast_discoverable_metrics(df, Dict{String,Any}())
blast_discoverable_metrics(df::DataFrame, ::Nothing) = blast_discoverable_metrics(df)
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

function seq_to_name_lookup(db_seqs)
    lookup = Dict{String, String}()
    for (name, seq) in db_seqs
        haskey(lookup, seq) || (lookup[seq] = name)
    end
    return lookup
end

"""
    name_candidate(core, sseqid, core_aln_mismatch, db_seqs, isin) -> String

Name a discovery candidate from its trimmed `core`. "Novel" means genuine sequence variation
in the (non-extended) gene relative to EVERY known allele, so the resolution order is:

  1. `core_aln_mismatch == 0`  → the best-hit reference `sseqid` (core equals that reference).
  2. core identical to ANY known sequence → that allele — independent of which allele won the
     BLAST best-hit and of `isin`. A core equal to a reference *is* that allele, never novel.
  3. `isin` only: core is a substring of a known allele → that allele (the read covers only
     part of the gene, with no internal variation — still known, not novel).
  4. otherwise → a hashed novel name (`unique_name`).

A core that strictly *contains* a known allele plus extra bases is deliberately left novel:
the extra bases are variation in the gene region (e.g. a junction insertion), not coverage.
`db_seqs` is the un-extended base reference as `(name, sequence)` pairs.
"""
function name_candidate(core::AbstractString, sseqid, core_aln_mismatch, db_seqs, isin::Bool)
    name_candidate(core, sseqid, core_aln_mismatch, seq_to_name_lookup(db_seqs), db_seqs, isin)
end

function name_candidate(core::AbstractString, sseqid, core_aln_mismatch,
                        exact_lookup::Dict{String, String},
                        db_seqs, isin::Bool)
    # Empty cores are trim failures, not known alleles. `occursin("", seq)` is true for every
    # reference, which previously named failed trims as the first DB allele.
    isempty(core) && return unique_name(sseqid, core)
    core_aln_mismatch == 0 && return String(sseqid)
    haskey(exact_lookup, core) && return exact_lookup[core]
    if isin
        for (name, seq) in db_seqs
            isempty(seq) && continue
            occursin(core, seq) && return name
        end
    end
    return unique_name(sseqid, core)
end

"""
    handle_blast(parsed_args, immunediscover_module, always_gz)

Handle the blast search pipeline. Extracted from real_main for modularity.
"""
function handle_blast(parsed_args, immunediscover_module, always_gz)
    using_cli = immunediscover_module.Cli
    if get(parsed_args["discover"]["blast"], "show-presets", false)
        using_cli.show_blast_presets()
        return
    end
    parsed_args = using_cli.apply_blast_presets!(parsed_args)
    gene = parsed_args["discover"]["blast"]["gene"]
    haskey(using_cli.BLAST_PRESETS, gene) && @info "Applied $gene gene preset (overrides logged above)"
    # Grouped, ordered parameter display (final values after presets/overrides).
    params_report(parsed_args["discover"]["blast"], using_cli.BLAST_PARAM_GROUPS;
                  title="discover blast — parameters")
    @info "Discovery with BLAST assignments"

    fasta_path = parsed_args["discover"]["blast"]["fasta"]
    output_path = parsed_args["discover"]["blast"]["output"]
    work_dir = resolve_work_dir(parsed_args["discover"]["blast"]["work-dir"])
    isdir(work_dir) || mkpath(work_dir)
    @info "Work directory (caches, intermediates): $work_dir"
    fasta_stem = file_stem(fasta_path)
    affixes_path = joinpath(work_dir, fasta_stem * ".affixes")
    DB = immunediscover_module.load_fasta(fasta_path)

    verbose = parsed_args["discover"]["blast"]["verbose"]
    overwrite = parsed_args["discover"]["blast"]["overwrite"]
    minquality = parsed_args["discover"]["blast"]["minquality"]
    forward_extension = parsed_args["discover"]["blast"]["forward"]
    reverse_extension = parsed_args["discover"]["blast"]["reverse"]
    keep_failed = parsed_args["discover"]["blast"]["keep-failed"]
    min_corecov = parsed_args["discover"]["blast"]["min-corecov"]
    isin = parsed_args["discover"]["blast"]["isin"]

    if (forward_extension < 7) && (forward_extension > 0)
        @warn "Forward extension $forward_extension is short and may lead to false positives"
    end
    if (reverse_extension < 7) && (reverse_extension > 0)
        @warn "Reverse extension $reverse_extension is short and may lead to false positives"
    end

    # Process pseudo genes if provided
    db_p = Vector{Tuple{String, String}}()
    pseudo = parsed_args["discover"]["blast"]["pseudo"]
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

    # Handle sequence extension
    if forward_extension == 0 && reverse_extension == 0
        ext_fasta_path = combined_fasta_path
        @info "No sequence extension requested, using original sequences"
        empty_affixes = [(name="", prefix="", suffix="")]
        CSV.write(affixes_path, DataFrame(empty_affixes), delim='\t')
    else
        ext_fasta_path = joinpath(work_dir, fasta_stem * "-combined-extended.fasta")
        if !isfile(ext_fasta_path) || overwrite
            demux = load_csv(parsed_args["discover"]["blast"]["input"])
            @info "Extending gene sequences by $forward_extension forward and $reverse_extension reverse nucleotides" reads=nrow(demux) references=length(db_p) julia_threads=nthreads()
            # Extend the COMBINED db (real + pseudo) so the BLAST database matches the
            # no-extension branch. Pseudo decoys absent in reads become unextended
            # singletons but remain in the DB; using DB (real only) here silently dropped
            # them, disabling -p in the default extension mode.
            extended = accumulate_affixes(db_p, demux,
                forward_extension=forward_extension,
                reverse_extension=reverse_extension)
            affixes = save_extended(extended, ext_fasta_path)
            CSV.write(affixes_path, DataFrame(affixes, [:name, :prefix, :suffix]), delim='\t')
            @info "Saved affixes in $affixes_path"
        else
            @info "Using existing extended sequences from $ext_fasta_path"
        end
    end

    # Run BLAST discovery
    blast_clusters = blast_discover(
        parsed_args["discover"]["blast"]["input"],
        ext_fasta_path;
        work_dir=work_dir,
        max_dist=parsed_args["discover"]["blast"]["max-blast-mismatch"],
        min_edge=parsed_args["discover"]["blast"]["edge"],
        min_scov=parsed_args["discover"]["blast"]["subjectcov"],
        args=parsed_args["discover"]["blast"]["args"],
        verbose=verbose,
        overwrite=overwrite,
        min_read_length=get(parsed_args["discover"]["blast"], "min-read-length", 0),
    )

    # Trim extensions if applicable
    if forward_extension != 0 || reverse_extension != 0
        if isfile(affixes_path)
            @info "Loading affixes from $affixes_path"
            affix_dict = Dict{String, Tuple{String, String}}()
            for row in eachrow(CSV.File(affixes_path, delim='\t') |> DataFrame)
                name = strip(ismissing(row.name) ? "" : String(row.name))
                prefix = ismissing(row.prefix) ? "" : String(row.prefix)
                suffix = ismissing(row.suffix) ? "" : String(row.suffix)
                isempty(name) && continue
                affix_dict[name] = (prefix, suffix)
            end

            stats = AlignmentStats()
            check_affix_quality_warning(forward_extension, minquality)
            check_affix_quality_warning(reverse_extension, minquality)

            # Normalize keys with strip() so FASTA headers with trailing space match BLAST sseqid (first token)
            db_keys = Set(strip(String(name)) for (name, seq) in DB)
            db_dict = Dict(strip(String(name)) => String(seq) for (name, seq) in DB)
            db_lengths = Dict(strip(String(name)) => length(seq) for (name, seq) in DB)

            # Parallel map with thread-local stats; merge into stats after (no concurrent mutation)
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

            # Core coverage filter (use base key for novel alleles with _S suffix)
            blast_clusters[:, :db_length] = [get(db_lengths, sseqid_to_db_key(strip(String(row.sseqid)), db_keys), 0) for row in eachrow(blast_clusters)]
            blast_clusters[:, :corecov] = map(eachrow(blast_clusters)) do row
                row.db_length > 0 ? length(row.aln_qseq) / row.db_length : 0.0
            end
            before_corecov = nrow(blast_clusters)
            section("BLAST discovery — trimming & core coverage")
            # corecov = trimmed core length (aln_qseq, the segment left after removing the
            # prefix/suffix affixes) ÷ length of the matched DB allele. Annotate, don't drop.
            mark_rejected!(blast_clusters, blast_clusters.corecov .< min_corecov,
                           "core coverage < $min_corecov (core÷DB len)", "core coverage")
            kept_corecov = count(isempty, blast_clusters.reject_reason)
            stage_report("core coverage ≥ $min_corecov (core÷DB len)", kept_corecov, before_corecov;
                         values=blast_clusters.corecov, histogram=true, hist_closed=:right)

            @info "Alignment stats: $(stats.total_attempts) attempts, $(stats.prefix_failures) prefix, $(stats.suffix_failures) suffix failures"
            if verbose
                for (gene_name, messages) in stats.failed_genes
                    for msg in messages
                        @info "  Failed [$gene_name]: $msg"
                    end
                end
            end
        end
    else
        blast_clusters[:, :aln_qseq] = blast_clusters[:, :qseq]
        blast_clusters[:, :core_aln_mismatch] = blast_clusters[:, :blast_mismatch]
    end

    add_blast_support_columns!(blast_clusters)

    # Candidate-quality metrics (composition + cross-donor recurrence + support) — computed
    # before the output filters so they can also serve as optional filter criteria.
    #  - gc_content / max_homopolymer: composition of the trimmed core,
    #  - n_donors: distinct donors sharing the exact core (cross-donor recurrence),
    #  - n_reads_total: reads backing the core across the run.
    blast_clusters[:, :gc_content] = gc_content.(blast_clusters.aln_qseq)
    blast_clusters[:, :max_homopolymer] = max_homopolymer.(blast_clusters.aln_qseq)
    transform!(groupby(blast_clusters, :aln_qseq), :case => (x -> length(unique(x))) => :n_donors)
    transform!(groupby(blast_clusters, :aln_qseq), :full_count => sum => :n_reads_total)
    # Peak per-donor allelic ratio for the core: the highest fraction of its gene's reads
    # it reaches in any single donor. A germline allele is a major allele in ≥1 carrier;
    # artifacts stay minor everywhere — the most discriminative recall-safe separator.
    transform!(groupby(blast_clusters, :aln_qseq), FULL_ALLELIC_RATIO => maximum => PEAK_ALLELIC_RATIO)

    blast_block = parsed_args["discover"]["blast"]
    criteria = build_blast_output_criteria(blast_block; keep_failed=keep_failed)

    # Apply each output filter in turn, annotating (not dropping) so the full table keeps
    # every candidate, and report how many each filter removes — separately.
    section("BLAST discovery — output filters")
    init_rejection_columns!(blast_clusters)
    for criterion in criteria
        before = count(isempty, blast_clusters.reject_reason)
        fail = Bool[!passes(row, criterion) for row in eachrow(blast_clusters)]
        mark_rejected!(blast_clusters, fail, criterion.label, "output filter")
        stage_report(criterion.label, count(isempty, blast_clusters.reject_reason), before)
    end
    # Name a candidate. "Novel" means genuine sequence variation in the (non-extended) gene
    # region relative to EVERY known allele — so an exact match to any known reference is
    # always a known call, never a hashed novel name. The match is judged on the trimmed
    # core (aln_qseq) against the un-extended base sequences (db_seqs).
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

    # Between-cluster separation (diagnostic): for each accepted core, the nearest
    # more-abundant core ("parent"). Small nn_dist + large parent_ratio = error satellite.
    @info "Computing neighbor statistics"
    add_neighbor_stats!(blast_clusters)

    # Two outputs: the filtered table (candidates that passed every stage) and a full table
    # holding every candidate plus reject_reason / reject_stage, for inspection and tuning.
    reason_cols = [:reject_reason, :reject_stage]
    output = always_gz(parsed_args["discover"]["blast"]["output"])
    full_arg = get(parsed_args["discover"]["blast"], "full-output", nothing)
    full_output = (full_arg === nothing || isempty(full_arg)) ?
        always_gz(replace(replace(output, r"\.gz$" => ""), r"\.tsv$" => "") * ".full.tsv") :
        always_gz(full_arg)

    section("BLAST discovery — summary")
    kept = accepted(blast_clusters)
    stage_report("accepted (passed all filters)", nrow(kept), nrow(blast_clusters))
    report_rejections(blast_clusters.reject_reason)
    report_recurrence(kept)
    # Distance alone is not suspicious — most genuine novel alleles are 1 bp from a known
    # parent. What flags a likely error is being a small fraction of a much more abundant
    # neighbour, i.e. a large parent_ratio. Report both, but weight by parent_ratio.
    n_close = count((kept.nn_dist .>= 0) .& (kept.nn_dist .<= 1))
    n_shadow = count(kept.likely_satellite)
    n_close > 0 && @info "$n_close accepted candidate(s) are within 1 bp of a more-abundant core; $n_shadow flagged likely_satellite (satellite_score ≥ 0.5). See nn_dist / parent_ratio / satellite_score / chimera_score — not a default filter."
    # Novel-vs-parent SNP heatmaps (non-satellites only): one panel per gene.
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
