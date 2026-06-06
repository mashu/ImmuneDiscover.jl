module Selftest
    # Evaluate how well `discover blast` recovers known-novel alleles. Inputs:
    #   - a discovery FULL table (every candidate + reject_reason / reject_stage),
    #   - a BASE reference FASTA (the reference used for discovery),
    #   - a TRUTH FASTA (known + novel).
    # Truth-novel = sequences in TRUTH not in BASE — the alleles the algorithm should discover.
    # We report recall (recovered / truth-novel), precision (true / accepted-novel), and for each
    # truth-novel allele whether it was recovered, rejected (and by which stage), or never seen.

    using DataFrames
    using CSV
    using Statistics: median
    using ..Data: load_fasta, barplot_if_available
    using ..Report: section, stage_report

    export handle_selftest, evaluate_recovery, classify_allele, is_novel, metric_separation,
           recall_safe_filters

    # Metric columns worth scanning for separation (intersected with what the table actually has).
    const DEFAULT_METRICS = ["max_full_ratio", "n_donors", "n_reads_total", "scov", "corecov",
                             "aln_mismatch", "mismatch", "gc_content", "max_homopolymer",
                             "nn_dist", "parent_ratio", "full_count"]

    "True if `seq` is absent from the base (reference) set."
    is_novel(seq::AbstractString, base::AbstractSet) = !(seq in base)

    "Match a truth allele against a candidate core: exact, or substring either way when `substring`."
    function seq_match(a::AbstractString, b::AbstractString; substring::Bool=true)
        a == b && return true
        substring && (occursin(a, b) || occursin(b, a)) && return true
        return false
    end

    """
        classify_allele(seq, accepted_set, accepted_list, rejected; substring) -> (status, stage)

    `status` ∈ ("recovered", "rejected", "missed"). Recovered if an accepted candidate core
    matches `seq`; else "rejected" (with the stage) if a rejected candidate matches it exactly;
    else "missed".
    """
    function classify_allele(seq::AbstractString, accepted_set::AbstractSet,
                             accepted_list::AbstractVector{<:AbstractString},
                             rejected::AbstractDict; substring::Bool=true)
        seq in accepted_set && return ("recovered", "")
        if substring
            for c in accepted_list
                seq_match(seq, c; substring=true) && return ("recovered", "")
            end
        end
        haskey(rejected, seq) && return ("rejected", rejected[seq])
        return ("missed", "")
    end

    """
        evaluate_recovery(discovery, base, truth; seq_col, substring)

    Returns `(summary, per_allele, false_positives)`. `truth` is a vector of (name, seq) pairs.
    """
    function evaluate_recovery(discovery::DataFrame, base::AbstractSet{String},
                               truth::AbstractVector{<:Tuple{<:AbstractString,<:AbstractString}};
                               seq_col::Symbol=:aln_qseq, substring::Bool=true)
        cols = names(discovery)
        String(seq_col) in cols ||
            error("Discovery table has no '$seq_col' column — pass --seq-col, or use the full table (<output>.full.tsv.gz).")
        # Empty fields round-trip through CSV as `missing`; coalesce so an accepted row's empty
        # reject_reason isn't read as "missing" (which would look rejected), and drop empty cores.
        getcol(name) = name in cols ? [ismissing(x) ? "" : String(x) for x in discovery[!, name]] :
                                      fill("", nrow(discovery))
        reason = getcol("reject_reason")
        stage  = getcol("reject_stage")
        seqs   = [ismissing(x) ? "" : String(x) for x in discovery[!, seq_col]]

        acc_mask = isempty.(reason)
        accepted_list = filter(!isempty, unique(seqs[acc_mask]))
        accepted_set = Set(accepted_list)
        rejected = Dict{String,String}()
        for i in eachindex(seqs)
            (acc_mask[i] || isempty(seqs[i])) && continue
            get!(rejected, seqs[i], stage[i])
        end

        truth_novel = [(String(n), String(s)) for (n, s) in truth if is_novel(s, base)]
        Row = NamedTuple{(:allele, :length, :status, :reject_stage), Tuple{String,Int,String,String}}
        rows = Row[]
        recovered = 0
        for (name, s) in truth_novel
            status, stg = classify_allele(s, accepted_set, accepted_list, rejected; substring=substring)
            status == "recovered" && (recovered += 1)
            push!(rows, (allele=name, length=length(s), status=status, reject_stage=stg))
        end
        per_allele = isempty(rows) ? DataFrame(allele=String[], length=Int[], status=String[], reject_stage=String[]) :
                                     DataFrame(rows)

        # Precision over accepted NOVEL candidate cores (cores not present in base).
        truth_set = Set(s for (_, s) in truth_novel)
        novel_cores = [c for c in accepted_list if is_novel(c, base)]
        tp = 0
        false_positives = String[]
        for c in novel_cores
            hit = c in truth_set || (substring && any(t -> seq_match(c, t; substring=true), truth_set))
            hit ? (tp += 1) : push!(false_positives, c)
        end

        n_truth = length(truth_novel)
        n_novel_acc = length(novel_cores)
        summary = (n_truth_novel = n_truth,
                   recovered = recovered,
                   recall = n_truth == 0 ? 0.0 : recovered / n_truth,
                   n_novel_accepted = n_novel_acc,
                   true_positive = tp,
                   false_positive = length(false_positives),
                   precision = n_novel_acc == 0 ? 0.0 : tp / n_novel_acc)
        return (summary = summary, per_allele = per_allele, false_positives = false_positives)
    end

    "Coerce a discovery column to Float64, mapping missing / unparseable entries to NaN."
    function numeric_col(df::DataFrame, name::AbstractString)
        raw = df[!, name]
        out = Vector{Float64}(undef, length(raw))
        @inbounds for (i, x) in enumerate(raw)
            out[i] = ismissing(x) ? NaN :
                     x isa Real   ? Float64(x) :
                     (v = tryparse(Float64, String(x)); v === nothing ? NaN : v)
        end
        return out
    end

    """
        best_threshold(tp, fp) -> (direction, threshold, youden, tp_kept, fp_kept)

    Over the observed cut points and both keep-directions, find the threshold that best
    separates true (`tp`) from false (`fp`) values by Youden's J (= TP-rate − FP-rate).
    `direction` is `:ge` (keep ≥ threshold) or `:le` (keep ≤ threshold).
    """
    function best_threshold(tp::AbstractVector{<:Real}, fp::AbstractVector{<:Real})
        cuts = sort!(unique(vcat(collect(tp), collect(fp))))
        ntp, nfp = length(tp), length(fp)
        best = (:ge, NaN, -Inf, 0, 0)
        for thr in cuts, dir in (:ge, :le)
            keep_tp = dir === :ge ? count(>=(thr), tp) : count(<=(thr), tp)
            keep_fp = dir === :ge ? count(>=(thr), fp) : count(<=(thr), fp)
            j = keep_tp / ntp - keep_fp / nfp
            j > best[3] && (best = (dir, Float64(thr), j, keep_tp, keep_fp))
        end
        return best
    end

    """
        metric_separation(discovery, base, truth; seq_col, substring, metrics) -> DataFrame

    Label every NOVEL candidate row (core absent from `base`) as true (`tp`, its core matches a
    truth-novel allele) or false (`fp`). For each numeric `metric` column present, report how
    well it separates true from false and the single threshold that does it best (Youden's J).
    Sorted by `youden` descending — the top row is the most discriminative metric to filter on.
    """
    function metric_separation(discovery::DataFrame, base::AbstractSet{String},
                               truth::AbstractVector{<:Tuple{<:AbstractString,<:AbstractString}};
                               seq_col::Symbol=:aln_qseq, substring::Bool=true,
                               metrics::AbstractVector{<:AbstractString}=DEFAULT_METRICS)
        cols = names(discovery)
        String(seq_col) in cols ||
            error("Discovery table has no '$seq_col' column — pass --seq-col, or use the full table (<output>.full.tsv.gz).")
        seqs = [ismissing(x) ? "" : String(x) for x in discovery[!, seq_col]]
        truth_novel = Set(String(s) for (_, s) in truth if is_novel(s, base))

        labels = fill(:none, length(seqs))
        for i in eachindex(seqs)
            c = seqs[i]
            (isempty(c) || !is_novel(c, base)) && continue
            hit = c in truth_novel || (substring && any(t -> seq_match(c, t; substring=true), truth_novel))
            labels[i] = hit ? :tp : :fp
        end
        tp_idx = findall(==(:tp), labels)
        fp_idx = findall(==(:fp), labels)

        Row = NamedTuple{(:metric, :n_tp, :n_fp, :tp_median, :fp_median, :direction,
                          :threshold, :youden, :tp_kept, :fp_removed),
                         Tuple{String,Int,Int,Float64,Float64,String,Float64,Float64,Int,Int}}
        rows = Row[]
        for m in metrics
            m in cols || continue
            vals = numeric_col(discovery, m)
            tp = filter(isfinite, vals[tp_idx])
            fp = filter(isfinite, vals[fp_idx])
            (isempty(tp) || isempty(fp)) && continue
            dir, thr, j, keep_tp, keep_fp = best_threshold(tp, fp)
            push!(rows, (metric=m, n_tp=length(tp), n_fp=length(fp),
                         tp_median=median(tp), fp_median=median(fp),
                         direction = dir === :ge ? "keep ≥" : "keep ≤",
                         threshold=thr, youden=j, tp_kept=keep_tp, fp_removed=length(fp) - keep_fp))
        end
        df = DataFrame(metric=String[], n_tp=Int[], n_fp=Int[], tp_median=Float64[],
                       fp_median=Float64[], direction=String[], threshold=Float64[],
                       youden=Float64[], tp_kept=Int[], fp_removed=Int[])
        isempty(rows) || (df = sort!(DataFrame(rows), :youden, rev=true))
        return df
    end

    """
        recall_safe_filters(discovery, base, truth; seq_col, substring, metrics) -> DataFrame

    The directly actionable companion to `metric_separation`. For each numeric metric, find the
    most aggressive single-sided threshold that can be applied to the ACCEPTED candidates while
    still recovering every truth-novel allele (recall stays 1.0), and report how many accepted
    false-positive novel cores that cut removes. "Apply this blast threshold to delete false
    positives at zero recall cost." Sorted by false positives removed, descending.
    """
    function recall_safe_filters(discovery::DataFrame, base::AbstractSet{String},
                                 truth::AbstractVector{<:Tuple{<:AbstractString,<:AbstractString}};
                                 seq_col::Symbol=:aln_qseq, substring::Bool=true,
                                 metrics::AbstractVector{<:AbstractString}=DEFAULT_METRICS)
        cols = names(discovery)
        String(seq_col) in cols ||
            error("Discovery table has no '$seq_col' column — pass --seq-col, or use the full table (<output>.full.tsv.gz).")
        getcol(name) = name in cols ? [ismissing(x) ? "" : String(x) for x in discovery[!, name]] :
                                      fill("", nrow(discovery))
        reason = getcol("reject_reason")
        seqs   = [ismissing(x) ? "" : String(x) for x in discovery[!, seq_col]]
        acc = isempty.(reason)

        truth_novel = [(String(n), String(s)) for (n, s) in truth if is_novel(s, base)]

        # Accepted rows grouped by novel core; label each core true/false vs truth-novel.
        core_rows = Dict{String,Vector{Int}}()
        for i in eachindex(seqs)
            (acc[i] && !isempty(seqs[i]) && is_novel(seqs[i], base)) || continue
            push!(get!(core_rows, seqs[i], Int[]), i)
        end
        isempty(core_rows) && return _empty_safe_df()
        cores = collect(keys(core_rows))
        is_tp = Dict(c => any(t -> seq_match(c, t[2]; substring=substring), truth_novel) for c in cores)
        # truth allele -> accepted cores that recover it (only recovered alleles constrain the cut)
        allele_cores = Dict{String,Vector{String}}()
        for (name, s) in truth_novel
            matched = [c for c in cores if seq_match(c, s; substring=substring)]
            isempty(matched) || (allele_cores[name] = matched)
        end
        fp_cores = [c for c in cores if !is_tp[c]]

        Row = NamedTuple{(:metric, :acc_fp, :direction, :threshold, :fp_removed),
                         Tuple{String,Int,String,Float64,Int}}
        rows = Row[]
        for m in metrics
            m in cols || continue
            vals = numeric_col(discovery, m)
            # per-core favourable values: max keeps a core under "keep ≥", min under "keep ≤".
            cmax = Dict{String,Float64}(); cmin = Dict{String,Float64}()
            for (c, idx) in core_rows
                fv = filter(isfinite, vals[idx])
                isempty(fv) && continue
                cmax[c] = maximum(fv); cmin[c] = minimum(fv)
            end
            (isempty(allele_cores) || isempty(fp_cores)) && continue
            # keep ≥ t: an allele survives if any of its cores has cmax ≥ t; safe t = min over
            # alleles of the allele's best cmax. keep ≤ t: symmetric with cmin.
            allele_ge = [maximum(cmax[c] for c in cs if haskey(cmax, c); init=-Inf) for cs in values(allele_cores)]
            allele_le = [minimum(cmin[c] for c in cs if haskey(cmin, c); init=Inf)  for cs in values(allele_cores)]
            safe_ge = minimum(allele_ge); safe_le = maximum(allele_le)
            fp_ge = count(c -> haskey(cmax, c) && cmax[c] < safe_ge, fp_cores)
            fp_le = count(c -> haskey(cmin, c) && cmin[c] > safe_le, fp_cores)
            if fp_ge >= fp_le
                push!(rows, (metric=m, acc_fp=length(fp_cores), direction="keep ≥", threshold=safe_ge, fp_removed=fp_ge))
            else
                push!(rows, (metric=m, acc_fp=length(fp_cores), direction="keep ≤", threshold=safe_le, fp_removed=fp_le))
            end
        end
        isempty(rows) && return _empty_safe_df()
        return sort!(DataFrame(rows), :fp_removed, rev=true)
    end

    _empty_safe_df() = DataFrame(metric=String[], acc_fp=Int[], direction=String[],
                                 threshold=Float64[], fp_removed=Int[])

    "Print the recall-safe filter table: cuts that drop false positives without losing alleles."
    function report_recall_safe(safe::DataFrame)
        section("Recall-safe filters — drop false positives without losing any truth-novel allele")
        if nrow(safe) == 0
            @info "No recall-safe cut found (need both recovered truth alleles and accepted false novel cores)."
            return
        end
        printstyled(rpad("metric", 16), rpad("safe rule", 18), "FP removed (recall stays 1.0)\n";
                    color=:cyan, bold=true)
        for r in eachrow(safe)
            rule = "$(r.direction) $(round(r.threshold; digits=4))"
            col = r.fp_removed > 0 ? :green : :light_black
            printstyled("  ", rpad(r.metric, 14); color=col, bold=true)
            println(rpad(rule, 18), "$(r.fp_removed)/$(r.acc_fp)")
        end
        println("  Apply the top rule as a blast threshold to cut false positives at zero recall cost.")
        return nothing
    end

    "Print the metric-separation table with a colored header and per-metric rows."
    function report_separation(sep::DataFrame)
        section("Metric separation — threshold that best splits true from false novel candidates")
        if nrow(sep) == 0
            @info "Not enough labelled candidates (need both true and false novel rows with metric columns)."
            return
        end
        printstyled(rpad("metric", 16), rpad("rule", 14), rpad("youden", 9),
                    rpad("true kept", 12), rpad("false dropped", 15), "median T | F\n";
                    color=:cyan, bold=true)
        for r in eachrow(sep)
            rule = "$(r.direction) $(round(r.threshold; digits=3))"
            col = r.youden >= 0.5 ? :green : r.youden >= 0.25 ? :yellow : :light_black
            printstyled("  ", rpad(r.metric, 14); color=col, bold=true)
            print(rpad(rule, 14),
                  rpad(round(r.youden; digits=3), 9),
                  rpad("$(r.tp_kept)/$(r.n_tp)", 12),
                  rpad("$(r.fp_removed)/$(r.n_fp)", 15),
                  "$(round(r.tp_median; digits=3)) | $(round(r.fp_median; digits=3))\n")
        end
        println("  Higher youden ⇒ better separation; apply the matching blast threshold to keep true and drop false.")
    end

    function handle_selftest(parsed_args)
        b = parsed_args["discover"]["selftest"]
        section("Self-test — recovery of novel alleles")
        discovery = CSV.File(b["discovery"], delim='\t') |> DataFrame
        base = Set(String(s) for (_, s) in load_fasta(b["base"], validate=false))
        truth = load_fasta(b["truth"], validate=false)
        seq_col = Symbol(get(b, "seq-col", "aln_qseq"))
        substring = !get(b, "no-substring", false)

        res = evaluate_recovery(discovery, base, truth; seq_col=seq_col, substring=substring)
        s = res.summary

        stage_report("recovered (recall)", s.recovered, s.n_truth_novel)
        printstyled("  recall    = "; color=:cyan, bold=true);  println(round(s.recall;    digits=3), "  ($(s.recovered)/$(s.n_truth_novel) truth-novel)")
        printstyled("  precision = "; color=:cyan, bold=true);  println(round(s.precision; digits=3), "  (TP=$(s.true_positive), FP=$(s.false_positive) of $(s.n_novel_accepted) accepted novel)")

        if nrow(res.per_allele) > 0
            tally = Dict{String,Int}()
            for st in res.per_allele.status
                tally[st] = get(tally, st, 0) + 1
            end
            ks = sort(collect(keys(tally)))
            println("  truth-novel alleles by outcome:")
            barplot_if_available(ks, [tally[k] for k in ks])
            # which stage rejected the ones we saw but dropped
            rej = filter(r -> r.status == "rejected", res.per_allele)
            if nrow(rej) > 0
                st = Dict{String,Int}()
                for r in eachrow(rej); st[r.reject_stage] = get(st, r.reject_stage, 0) + 1; end
                println("  rejected truth-novel by stage (which filter to relax):")
                rk = sort(collect(keys(st)))
                barplot_if_available(rk, [st[k] for k in rk])
            end
        end

        safe = recall_safe_filters(discovery, base, truth; seq_col=seq_col, substring=substring)
        report_recall_safe(safe)

        sep = metric_separation(discovery, base, truth; seq_col=seq_col, substring=substring)
        report_separation(sep)

        CSV.write(b["output"], res.per_allele, delim='\t')
        @info "Per-allele recovery table saved to $(b["output"]) ($(nrow(res.per_allele)) truth-novel alleles)"
        metrics_out = get(b, "metrics-output", nothing)
        if metrics_out !== nothing && !isempty(metrics_out)
            CSV.write(metrics_out, sep, delim='\t')
            @info "Metric-separation table saved to $metrics_out ($(nrow(sep)) metrics)"
        end
        return (recovery = res, separation = sep, recall_safe = safe)
    end
end
