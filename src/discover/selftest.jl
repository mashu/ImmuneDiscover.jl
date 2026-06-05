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
    using ..Data: load_fasta, barplot_if_available
    using ..Report: section, stage_report

    export handle_selftest, evaluate_recovery, classify_allele, is_novel

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
        reason = "reject_reason" in cols ? string.(discovery.reject_reason) : fill("", nrow(discovery))
        stage  = "reject_stage"  in cols ? string.(discovery.reject_stage)  : fill("", nrow(discovery))
        seqs = String.(discovery[!, seq_col])

        acc_mask = isempty.(reason)
        accepted_list = unique(seqs[acc_mask])
        accepted_set = Set(accepted_list)
        rejected = Dict{String,String}()
        for i in eachindex(seqs)
            acc_mask[i] && continue
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

        CSV.write(b["output"], res.per_allele, delim='\t')
        @info "Per-allele recovery table saved to $(b["output"]) ($(nrow(res.per_allele)) truth-novel alleles)"
        return res
    end
end
