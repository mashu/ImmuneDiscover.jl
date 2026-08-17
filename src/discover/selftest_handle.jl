function report_recall_safe(safe::DataFrame)
    section("Recall-safe filters — drop false positives without losing any truth-novel allele")
    if nrow(safe) == 0
        @info "No recall-safe cut found (need both recovered truth alleles and accepted false novel cores)."
        return
    end
    printstyled(rpad("metric", 16), rpad("blast flag", 28), "FP removed (recall stays 1.0)\n";
                color=:cyan, bold=true)
    for r in eachrow(safe)
        col = r.fp_removed > 0 ? :green : :light_black
        printstyled("  ", rpad(r.metric, 14); color=col, bold=true)
        println(rpad(r.cli_suggestion, 28), "$(r.fp_removed)/$(r.acc_fp)")
    end
    println("  Each metric tested in isolation on accepted novel cores (not the sequential filter cascade).")
    println("  Apply the top blast flag to cut false positives at zero recall cost.")
    return nothing
end

function report_marginal_shadow(shadow::Dict{String,Int})
    isempty(shadow) && return
    section("Marginal filter shadowing — other filters that would also fail (isolated, at run thresholds)")
    println("  discover blast applies output filters sequentially; reject_reason is the first failure only.")
    ks = sort(collect(keys(shadow)), by=k -> -shadow[k])
    barplot_if_available(ks, [shadow[k] for k in ks])
end

function report_separation(sep::DataFrame)
    section("Metric separation — threshold that best splits true from false novel candidates")
    if nrow(sep) == 0
        @info "Not enough labelled candidates (need both true and false novel rows with metric columns)."
        return
    end
    printstyled(rpad("metric", 16), rpad("blast flag", 22), rpad("youden", 9),
                rpad("true kept", 12), rpad("false dropped", 15), "median T | F\n";
                color=:cyan, bold=true)
    for r in eachrow(sep)
        col = r.youden >= 0.5 ? :green : r.youden >= 0.25 ? :yellow : :light_black
        printstyled("  ", rpad(r.metric, 14); color=col, bold=true)
        print(rpad(r.cli_suggestion, 22),
              rpad(round(r.youden; digits=3), 9),
              rpad("$(r.tp_kept)/$(r.n_tp)", 12),
              rpad("$(r.fp_removed)/$(r.n_fp)", 15),
              "$(round(r.tp_median; digits=3)) | $(round(r.fp_median; digits=3))\n")
    end
    println("  Each metric tested in isolation on all novel candidate rows (accepted + rejected).")
    println("  Higher youden ⇒ better separation; apply the blast flag on the top row to tune filters.")
end

function handle_selftest(parsed_args)
    b = parsed_args["discover"]["selftest"]
    section("Self-test — recovery of novel alleles")
    discovery = CSV.File(b["discovery"], delim='\t') |> DataFrame
    base = Set(String(s) for (_, s) in load_fasta(b["base"]))
    truth = load_fasta(b["truth"])
    seq_col = Symbol(get(b, "seq-col", "aln_qseq"))
    substring = !get(b, "no-substring", false)
    blast_block = selftest_blast_block(parsed_args)
    gene = get(b, "gene", "")
    n_metrics = length(blast_discoverable_metrics(discovery, blast_block))
    @info "Scanning $n_metrics blast metric column(s)$(isempty(gene) ? "" : " (gene preset $gene for marginal audit)")"

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
        rej = filter(r -> r.status == "rejected", res.per_allele)
        if nrow(rej) > 0
            st = Dict{String,Int}()
            for r in eachrow(rej); st[r.reject_stage] = get(st, r.reject_stage, 0) + 1; end
            println("  rejected truth-novel by stage (which filter to relax):")
            rk = sort(collect(keys(st)))
            barplot_if_available(rk, [st[k] for k in rk])
            nonempty = filter(r -> !isempty(r.reject_reason), rej)
            if nrow(nonempty) > 0
                println("  rejected truth-novel by filter (see per-allele reject_reason column):")
                rs = Dict{String,Int}()
                for r in eachrow(nonempty); rs[r.reject_reason] = get(rs, r.reject_reason, 0) + 1; end
                rk2 = sort(collect(keys(rs)))
                barplot_if_available(rk2, [rs[k] for k in rk2])
            end
        end
        missed = count(==( "missed"), res.per_allele.status)
        missed > 0 && println("  missed: $missed truth-novel allele(s) never reached the full table (e.g. failed --subjectcov during BLAST clustering).")
        if !isempty(gene)
            shadow = marginal_filter_shadowing(discovery, base, truth;
                                               seq_col=seq_col, substring=substring,
                                               blast_block=blast_block)
            report_marginal_shadow(shadow)
        elseif nrow(filter(r -> r.status == "rejected", res.per_allele)) > 0
            println("  Pass -g V|D|J to audit which other output filters would also fail (marginal / isolated).")
        end
    end

    safe = recall_safe_filters(discovery, base, truth; seq_col=seq_col, substring=substring,
                               blast_block=blast_block)
    report_recall_safe(safe)

    sep = metric_separation(discovery, base, truth; seq_col=seq_col, substring=substring,
                            blast_block=blast_block)
    report_separation(sep)

    CSV.write(b["output"], res.per_allele, delim='\t')
    @info "Per-allele recovery table saved to $(b["output"]) ($(nrow(res.per_allele)) truth-novel alleles)"
    write_metrics_table(optional(get(b, "metrics-output", nothing)), sep)
    return (recovery = res, separation = sep, recall_safe = safe)
end

write_metrics_table(::Absent, _) = nothing
function write_metrics_table(path::Present, sep)
    CSV.write(path.value, sep, delim='\t')
    @info "Metric-separation table saved to $(path.value) ($(nrow(sep)) metrics)"
    return nothing
end
