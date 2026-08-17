"""
    selftest_blast_block(parsed_args) -> Dict

Reconstruct the `discover blast` parameter block for marginal filter audit. Merges
`BLAST_DEFAULTS` with the optional `-g` preset (same keys as blast, no user overrides).
"""
function selftest_blast_block(parsed_args)
    b = parsed_args["discover"]["selftest"]
    block = Dict{String,Any}(k => v for (k, v) in BLAST_DEFAULTS)
    gene = get(b, "gene", "")
    if !isempty(gene) && haskey(BLAST_PRESETS, gene)
        for (k, v) in BLAST_PRESETS[gene]
            block[k] = v
        end
    end
    return block
end

"True if `seq` is absent from the base (reference) set."
is_novel(seq::AbstractString, base::AbstractSet) = !(seq in base)

"Match a truth allele against a candidate core: exact, or substring either way when `substring`."
function seq_match(a::AbstractString, b::AbstractString; substring::Bool=true)
    a == b && return true
    substring && (occursin(a, b) || occursin(b, a)) && return true
    return false
end

"""
    classify_allele(seq, accepted_set, accepted_list, rejected_stage, rejected_reason; substring)
        -> (status, reject_stage, reject_reason)

`status` ∈ ("recovered", "rejected", "missed"). Recovered if an accepted candidate core
matches `seq`; else "rejected" if a rejected candidate matches; else "missed".
Candidates that failed early cluster filters (e.g. `--subjectcov`) never appear and count as missed.
"""
function classify_allele(seq::AbstractString, accepted_set::AbstractSet,
                         accepted_list::AbstractVector{<:AbstractString},
                         rejected_stage::AbstractDict, rejected_reason::AbstractDict;
                         substring::Bool=true)
    seq in accepted_set && return ("recovered", "", "")
    if substring
        for c in accepted_list
            seq_match(seq, c; substring=true) && return ("recovered", "", "")
        end
    end
    haskey(rejected_stage, seq) &&
        return ("rejected", rejected_stage[seq], get(rejected_reason, seq, ""))
    return ("missed", "", "")
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
    getcol(name) = name in cols ? [ismissing(x) ? "" : String(x) for x in discovery[!, name]] :
                                  fill("", nrow(discovery))
    reason = getcol("reject_reason")
    stage  = getcol("reject_stage")
    seqs   = [ismissing(x) ? "" : String(x) for x in discovery[!, seq_col]]

    acc_mask = isempty.(reason)
    accepted_list = filter(!isempty, unique(seqs[acc_mask]))
    accepted_set = Set(accepted_list)
    rejected_stage = Dict{String,String}()
    rejected_reason = Dict{String,String}()
    for i in eachindex(seqs)
        (acc_mask[i] || isempty(seqs[i])) && continue
        get!(rejected_stage, seqs[i], stage[i])
        get!(rejected_reason, seqs[i], reason[i])
    end

    truth_novel = [(String(n), String(s)) for (n, s) in truth if is_novel(s, base)]
    Row = NamedTuple{(:allele, :length, :status, :reject_stage, :reject_reason),
                     Tuple{String,Int,String,String,String}}
    rows = Row[]
    recovered = 0
    for (name, s) in truth_novel
        status, stg, rsn = classify_allele(s, accepted_set, accepted_list,
                                           rejected_stage, rejected_reason; substring=substring)
        status == "recovered" && (recovered += 1)
        push!(rows, (allele=name, length=length(s), status=status, reject_stage=stg, reject_reason=rsn))
    end
    per_allele = isempty(rows) ? DataFrame(allele=String[], length=Int[], status=String[],
                                           reject_stage=String[], reject_reason=String[]) :
                                 DataFrame(rows)

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
