module HSMM

using Distributions
using CSV
using DataFrames
using Logging
using ProgressMeter
using Folds
using StringDistances

# Shared modules — avoids duplicate type definitions from repeated include()
using ..Data
using ..Exact
using ..Filters: FilterCriterion, MinThreshold, add_group_ratio!,
                 init_rejection_columns!, mark_rejected!, accepted, passes
using ..Report: section, stage_report, report_rejections
import ..Data: unique_name

@inline function dna_index(c::Char)::Int
    c == 'A' && return 1; c == 'C' && return 2; c == 'G' && return 3; c == 'T' && return 4; return 0
end
encode_dna(xs::AbstractString) = [dna_index(c) for c in xs]

abstract type AbstractEmissionModel end
struct IIDLogEmission <: AbstractEmissionModel; logp::NTuple{4,Float64}; end
@inline logemit(e::IIDLogEmission, base::Int)::Float64 = e.logp[base]

struct PWMLogEmission{L} <: AbstractEmissionModel; logp::Matrix{Float64}; end
@inline function logemit(e::PWMLogEmission{L}, pos::Int, base::Int)::Float64 where {L}
    @inbounds return e.logp[pos, base]
end

abstract type AbstractDurationModel end
struct DiscreteDuration <: AbstractDurationModel; min_len::Int; max_len::Int; logp::Vector{Float64}; end
@inline function logpmf(d::DiscreteDuration, len::Int)::Float64; @inbounds return d.logp[len - d.min_len + 1]; end

function discretize_truncated_normal(min_len::Int, max_len::Int, μ::Float64, σ::Float64)::DiscreteDuration
    @assert 0 < σ < Inf; @assert min_len <= max_len
    base = truncated(Normal(μ, σ), min_len, max_len)
    n = max_len - min_len + 1; pmf = Vector{Float64}(undef, n)
    for i in 1:n
        k = min_len + (i - 1); pmf[i] = max(cdf(base, k + 0.5) - cdf(base, k - 0.5), eps())
    end
    s = sum(pmf); @inbounds for i in eachindex(pmf); pmf[i] = log(pmf[i] / s); end
    return DiscreteDuration(min_len, max_len, pmf)
end

abstract type AbstractHSMMState end
struct FixedMotif{L,E<:AbstractEmissionModel} <: AbstractHSMMState; emission::E; end
Base.length(::FixedMotif{L}) where {L} = L
struct DGeneState{E<:AbstractEmissionModel,D<:AbstractDurationModel} <: AbstractHSMMState; emission::E; duration::D; end
@inline minduration(s::DGeneState) = s.duration.min_len
@inline maxduration(s::DGeneState) = s.duration.max_len

struct DGeneRSSHSMM{E9<:AbstractEmissionModel,E12<:AbstractEmissionModel,E7<:AbstractEmissionModel,EG<:AbstractEmissionModel,D<:AbstractDurationModel}
    pre_nonamer::FixedMotif{9,E9}; pre_spacer::FixedMotif{12,E12}; pre_heptamer::FixedMotif{7,E7}
    gene::DGeneState{EG,D}
    post_heptamer::FixedMotif{7,E7}; post_spacer::FixedMotif{12,E12}; post_nonamer::FixedMotif{9,E9}
end
@inline mingene(m::DGeneRSSHSMM) = minduration(m.gene)
@inline maxgene(m::DGeneRSSHSMM) = maxduration(m.gene)

function build_pwm_log(sequences::Vector{<:AbstractString}; pseudocount::Float64=0.01)
    @assert !isempty(sequences); L = length(sequences[1])
    @assert all(length(s) == L for s in sequences)
    counts = fill(pseudocount, L, 4)
    for seq in sequences; for (pos, ch) in enumerate(seq)
        idx = dna_index(ch); idx > 0 && (@inbounds counts[pos, idx] += 1.0)
    end; end
    logp = similar(counts)
    @inbounds for i in 1:L; s = sum(counts[i, :]); invs = 1.0 / s
        for b in 1:4; logp[i, b] = log(counts[i, b] * invs); end; end
    return logp
end

function estimate_iid_log_emission(sequences::Vector{<:AbstractString}; pseudocount::Float64=0.01)::IIDLogEmission
    counts = fill(Float64(pseudocount), 4)
    for s in sequences
        for ch in s
            idx = dna_index(ch)
            if 1 <= idx <= 4
                @inbounds counts[idx] += 1.0
            end
        end
    end
    total = sum(counts)
    invt = 1.0 / total
    return IIDLogEmission((log(counts[1] * invt), log(counts[2] * invt), log(counts[3] * invt), log(counts[4] * invt)))
end

function empirical_duration(seqs::Vector{<:AbstractString}, min_len::Int, max_len::Int; pseudocount::Float64=0.1)::DiscreteDuration
    n = max_len - min_len + 1; counts = fill(pseudocount, n)
    for s in seqs; L=length(s); min_len<=L<=max_len && (@inbounds counts[L-min_len+1]+=1.0); end
    total=sum(counts); invt=1.0/total; @inbounds for i in eachindex(counts); counts[i]=log(counts[i]*invt); end
    return DiscreteDuration(min_len, max_len, counts)
end

function fit_dgene_rss_hsmm(data_tuples::Vector{<:NTuple{7,AbstractString}}, min_gene::Int, max_gene::Int;
    pseudocount::Float64=0.01, duration::Symbol=:empirical,
    duration_mean::Float64=(min_gene+max_gene)/2, duration_std::Float64=max(1.0,(max_gene-min_gene)/6))
    pn=String[]; ps=String[]; ph=String[]; dg=String[]; oh=String[]; os=String[]; on=String[]
    for t in data_tuples; push!(pn,t[1]);push!(ps,t[2]);push!(ph,t[3]);push!(dg,t[4]);push!(oh,t[5]);push!(os,t[6]);push!(on,t[7]); end
    @assert all(length(s)==9 for s in pn); @assert all(length(s)==12 for s in ps); @assert all(length(s)==7 for s in ph)
    @assert all(length(s)==7 for s in oh); @assert all(length(s)==12 for s in os); @assert all(length(s)==9 for s in on)
    ge = estimate_iid_log_emission(dg; pseudocount)
    dur = duration===:empirical ? empirical_duration(dg,min_gene,max_gene;pseudocount=max(pseudocount,1e-6)) : discretize_truncated_normal(min_gene,max_gene,duration_mean,duration_std)
    return DGeneRSSHSMM(
        FixedMotif{9,PWMLogEmission{9}}(PWMLogEmission{9}(build_pwm_log(pn;pseudocount))),
        FixedMotif{12,PWMLogEmission{12}}(PWMLogEmission{12}(build_pwm_log(ps;pseudocount))),
        FixedMotif{7,PWMLogEmission{7}}(PWMLogEmission{7}(build_pwm_log(ph;pseudocount))),
        DGeneState{IIDLogEmission,DiscreteDuration}(ge, dur),
        FixedMotif{7,PWMLogEmission{7}}(PWMLogEmission{7}(build_pwm_log(oh;pseudocount))),
        FixedMotif{12,PWMLogEmission{12}}(PWMLogEmission{12}(build_pwm_log(os;pseudocount))),
        FixedMotif{9,PWMLogEmission{9}}(PWMLogEmission{9}(build_pwm_log(on;pseudocount))))
end

@inline function score_motif(obs::Vector{Int}, sp::Int, s::FixedMotif{L,PWMLogEmission{L}})::Float64 where {L}
    total=0.0; @inbounds for i in 1:L; b=obs[sp+i-1]; b==0 && return -Inf; total+=logemit(s.emission,i,b); end; total
end
@inline function score_gene_segment(obs::Vector{Int}, sp::Int, d::Int, s::DGeneState{IIDLogEmission})::Float64
    logp = s.emission.logp
    total = logpmf(s.duration, d)
    @inbounds for i in 0:(d-1)
        b = obs[sp + i]
        (b < 1 || b > 4) && return -Inf
        total += logp[b]
    end
    return total
end
@inline function logaddexp(a::Float64,b::Float64)::Float64; a==-Inf && return b; b==-Inf && return a; m=max(a,b); m+log1p(exp(min(a,b)-m)); end
@inline function pwm_logprob(seq::AbstractString, pwm::PWMLogEmission{L})::Float64 where {L}
    @assert length(seq)==L; total=0.0; @inbounds for (p,ch) in enumerate(seq); idx=dna_index(ch); idx==0 && return -Inf; total+=logemit(pwm,p,idx); end; total
end

function scan_best_and_total(sequence::AbstractString, model::DGeneRSSHSMM)
    obs=encode_dna(sequence); T=length(obs); k_pre=28; k_post=28
    best_log=-Inf; total_log=-Inf
    best=(prefix_start=0,prefix_end=0,gene_start=0,gene_end=0,suffix_start=0,suffix_end=0)
    min_g=mingene(model); max_g=maxgene(model)
    @inbounds for ps in 1:(T-(k_pre+min_g+k_post)+1)
        n9s=ps; s12s=n9s+9; h7s=s12s+12; pe=h7s+6
        lp=score_motif(obs,n9s,model.pre_nonamer); lp+=score_motif(obs,s12s,model.pre_spacer); lp+=score_motif(obs,h7s,model.pre_heptamer); lp==-Inf && continue
        gs=pe+1
        for d in min_g:max_g; ge=gs+d-1; ss=ge+1; se=ss+k_post-1; se>T && break
            lg=score_gene_segment(obs,gs,d,model.gene); lg==-Inf && continue
            h7p=ss; s12p=h7p+7; n9p=s12p+12
            ls=score_motif(obs,h7p,model.post_heptamer)+score_motif(obs,s12p,model.post_spacer)+score_motif(obs,n9p,model.post_nonamer); ls==-Inf && continue
            tot=lp+lg+ls; total_log=logaddexp(total_log,tot)
            if tot>best_log; best_log=tot; best=(prefix_start=n9s,prefix_end=pe,gene_start=gs,gene_end=ge,suffix_start=ss,suffix_end=se); end
        end
    end
    best_log==-Inf && return (gene_seq="",prefix_start=0,prefix_end=0,gene_start=0,gene_end=0,suffix_start=0,suffix_end=0,log_path_prob=-Inf,log_total_prob=-Inf,posterior_prob=0.0)
    gene_seq=sequence[best.gene_start:best.gene_end]; post=best_log-total_log; posterior=isfinite(post) ? exp(post) : 0.0
    return (gene_seq=gene_seq,prefix_start=best.prefix_start,prefix_end=best.prefix_end,gene_start=best.gene_start,gene_end=best.gene_end,suffix_start=best.suffix_start,suffix_end=best.suffix_end,log_path_prob=best_log,log_total_prob=total_log,posterior_prob=posterior)
end

function extract_dgene(sequence::AbstractString, model::DGeneRSSHSMM)
    r = scan_best_and_total(sequence, model)
    r.gene_seq == "" && return (gene_seq="", prefix_start=0, prefix_end=0, gene_start=0, gene_end=0, suffix_start=0, suffix_end=0, log_prob=-Inf)
    return (gene_seq=r.gene_seq, prefix_start=r.prefix_start, prefix_end=r.prefix_end, gene_start=r.gene_start, gene_end=r.gene_end, suffix_start=r.suffix_start, suffix_end=r.suffix_end, log_prob=r.log_path_prob)
end

    # Typed buffer row for HSMM scan results
    const HSMMScanRow = NamedTuple{
        (:well,:case,:sequence,:pre_nonamer,:pre_spacer,:pre_heptamer,
         :post_heptamer,:post_spacer,:post_nonamer,:heptamer_logp_pre,:heptamer_logp_post,
         :log_path_prob,:log_total_prob,:posterior_prob,:isin_db,:db_name,:nearest_db,:nearest_db_dist),
        Tuple{String,String,String,String,String,String,String,String,String,Float64,Float64,Float64,Float64,Float64,Bool,String,String,Int}}

    export AbstractEmissionModel, IIDLogEmission, PWMLogEmission, AbstractDurationModel, DiscreteDuration,
       AbstractHSMMState, FixedMotif, DGeneState, DGeneRSSHSMM,
       build_pwm_log, estimate_iid_log_emission, empirical_duration,
       discretize_truncated_normal, fit_dgene_rss_hsmm, extract_dgene,
       encode_dna, handle_hsmm

"""
    collapse_detections(res_df, min_posterior) -> DataFrame

Collapse per-detection HSMM rows to one row per (well, case, sequence). Each cluster is
represented by its BEST (max-posterior) detection, and `count` is the number of detections that
clear `min_posterior` — so the posterior threshold can be annotated downstream instead of
silently dropping detections, while `count` keeps meaning "confident detections".
"""
function collapse_detections(res_df::DataFrame, min_posterior::Real)
    combine(groupby(res_df, [:well, :case, :sequence])) do g
        bi = argmax(g.posterior_prob)
        (count = Base.count(>=(min_posterior), g.posterior_prob),
         pre_nonamer = g.pre_nonamer[bi], pre_spacer = g.pre_spacer[bi], pre_heptamer = g.pre_heptamer[bi],
         post_heptamer = g.post_heptamer[bi], post_spacer = g.post_spacer[bi], post_nonamer = g.post_nonamer[bi],
         heptamer_logp_pre = g.heptamer_logp_pre[bi], heptamer_logp_post = g.heptamer_logp_post[bi],
         log_path_prob = g.log_path_prob[bi], log_total_prob = g.log_total_prob[bi],
         posterior_prob = g.posterior_prob[bi],
         isin_db = any(g.isin_db), db_name = g.db_name[bi],
         nearest_db = g.nearest_db[bi], nearest_db_dist = minimum(g.nearest_db_dist))
    end
end

function run_hsmm(tsv::String, fasta_path::String, output::String;
    ratio::Float64=0.2, mincount::Int=5, min_gene_len::Int=0, max_gene_len::Int=0,
    limit::Int=0, min_posterior::Float64=0.7, out_mincount::Int=10, out_minratio::Float64=0.2,
    min_heptamer_prob_pre::Float64=0.0, min_heptamer_prob_post::Float64=0.0)
    @info "Loading demultiplex: $tsv"
    tbl = limit>0 ? Data.load_demultiplex(tsv,limit=limit) : Data.load_demultiplex(tsv)
    @info "Loaded $(nrow(tbl)) rows"
    @info "Loading D allele FASTA: $fasta_path"
    db = Data.load_fasta(fasta_path)
    db_seq_lookup = Dict{String,String}((seq=>name) for (name,seq) in db)
    db_names=first.(db); db_seqs=last.(db)
    @info "Searching known D alleles"
    # exact_search returns the unfiltered candidate table; keep D alleles with enough read
    # support for training (equivalent to the former mincount selection; no ratio filter here).
    known_df = Exact.exact_search(tbl, db, "D"; N=1000)
    nrow(known_df)==0 && (@warn "No exact D matches"; return DataFrame())
    filter!(r -> r.full_count >= mincount, known_df)
    nrow(known_df)==0 && (@warn "No D alleles passed mincount=$mincount"; return DataFrame())
    agg = combine(groupby(known_df, [:case,:gene,:sequence,:db_name]), :count=>sum=>:case_count)
    transform!(groupby(agg, [:case,:gene]), :case_count=>(x->x./maximum(x))=>:case_ratio)
    filter!(x->x.case_ratio>=ratio, agg)
    train_df = combine(groupby(agg, [:case,:gene])) do g
        sorted = sort(g, :case_count, rev=true); us = unique(sorted.sequence); kept = us[1:min(2,length(us))]
        filter(x->x.sequence in kept, sorted)
    end
    flank_cols = [:pre_nonamer,:pre_spacer,:pre_heptamer,:post_heptamer,:post_spacer,:post_nonamer]
    avail = intersect(flank_cols, Symbol.(names(known_df)))
    length(avail)<6 && (@warn "Not all RSS flank columns found"; return DataFrame())
    flank_df = unique(known_df[:, vcat(:db_name,:sequence,avail...)])
    train_df = leftjoin(train_df, flank_df, on=[:db_name,:sequence], makeunique=true)
    dropmissing!(train_df, avail)
    tuples = [(String(r.pre_nonamer),String(r.pre_spacer),String(r.pre_heptamer),String(r.sequence),String(r.post_heptamer),String(r.post_spacer),String(r.post_nonamer)) for r in eachrow(train_df)]
    lengths = length.(train_df.sequence)
    minL = min_gene_len>0 ? min_gene_len : minimum(lengths); maxL = max_gene_len>0 ? max_gene_len : maximum(lengths)
    @info "Training HSMM [$minL, $maxL] from $(length(tuples)) examples"
    model = fit_dgene_rss_hsmm(tuples, minL, maxL)
    pb=(n9=9,s12=12,h7=7); qb=(h7=7,s12=12,n9=9)
    @info "Scanning reads with HSMM"
    N_rows = nrow(tbl); chunk = max(10_000, min(100_000, Int(cld(N_rows,50))))
    starts = collect(1:chunk:N_rows); prog = Progress(length(starts), desc="Scanning chunks")
    buffer = HSMMScanRow[]
    for s in starts
        e = min(s+chunk-1, N_rows); sub = view(tbl, s:e, :)
        chunk_results = Folds.map(eachrow(sub)) do row
            seq = row.genomic_sequence; det = scan_best_and_total(seq, model); det.gene_seq=="" && return nothing
            ps_=det.prefix_start;pe_=det.prefix_end;ss_=det.suffix_start
            pn_=String(seq[ps_:ps_+pb.n9-1]); sp_=String(seq[ps_+pb.n9:ps_+pb.n9+pb.s12-1]); hp_=String(seq[ps_+pb.n9+pb.s12:pe_])
            oh_=String(seq[ss_:ss_+qb.h7-1]); os_=String(seq[ss_+qb.h7:ss_+qb.h7+qb.s12-1]); on_=String(seq[ss_+qb.h7+qb.s12:ss_+qb.h7+qb.s12+qb.n9-1])
            dseq=String(det.gene_seq); isin=haskey(db_seq_lookup,dseq); dn=isin ? db_seq_lookup[dseq] : "Novel"
            bd=typemax(Int); bc=""
            for (i,rs) in enumerate(db_seqs); d=evaluate(Levenshtein(),dseq,rs); d<bd && (bd=d;bc=db_names[i]; bd==0 && break); end
            (well=String(row.well),case=String(row.case),sequence=dseq,pre_nonamer=pn_,pre_spacer=sp_,pre_heptamer=hp_,
             post_heptamer=oh_,post_spacer=os_,post_nonamer=on_,heptamer_logp_pre=pwm_logprob(hp_,model.pre_heptamer.emission),
             heptamer_logp_post=pwm_logprob(oh_,model.post_heptamer.emission),log_path_prob=det.log_path_prob,
             log_total_prob=det.log_total_prob,posterior_prob=det.posterior_prob,isin_db=isin,db_name=dn,nearest_db=bc,nearest_db_dist=bd)
        end
        for r in chunk_results; r!==nothing && push!(buffer, r); end; next!(prog)
    end; finish!(prog)
    isempty(buffer) && (@warn "No D segments detected"; return DataFrame())
    res_df = DataFrame(buffer)
    @info "HSMM detections: $(nrow(res_df)) before collapse"
    # Collapse detections per (well, case, sequence). The posterior threshold is NOT applied here
    # as a hard drop; instead each cluster is represented by its BEST (max-posterior) detection and
    # carries count = number of detections clearing --min-posterior. The posterior is then annotated
    # downstream as a detection-stage filter, so every cluster — including low-posterior ones —
    # reaches the full table with a reason. count keeps its old meaning (confident detections).
    collapsed = collapse_detections(res_df, min_posterior)
    collapsed[!,:allele_name] = map(eachrow(collapsed)) do row
        is_known = (row.nearest_db_dist==0)||row.isin_db; bn = row.nearest_db!="" ? String(row.nearest_db) : String(row.db_name)
        is_known ? bn : unique_name(bn, String(row.sequence))
    end
    collapsed[!,:heptamer_prob_pre] = map(x->isfinite(x) ? exp(x) : 0.0, collapsed.heptamer_logp_pre)
    collapsed[!,:heptamer_prob_post] = map(x->isfinite(x) ? exp(x) : 0.0, collapsed.heptamer_logp_post)
    collapsed[!,:gene] = map(r->(n=String(r.nearest_db); n=="" ? "" : first(split(n,'*'))), eachrow(collapsed))
    init_rejection_columns!(collapsed)

    # Annotate (not drop) so the full table records why each detection was rejected. First reason
    # wins, so the HSMM-specific posterior filter is applied before the generic output filters.
    function annotate!(crits, stage)
        for criterion in crits
            before = count(isempty, collapsed.reject_reason)
            fail = Bool[!passes(row, criterion) for row in eachrow(collapsed)]
            mark_rejected!(collapsed, fail, criterion.label, stage)
            stage_report(criterion.label, count(isempty, collapsed.reject_reason), before)
        end
    end

    section("HSMM D detection — filters")
    # Detection-quality filter, specialized to the HSMM: the posterior probability of each
    # sequence's best detection must clear --min-posterior.
    annotate!(FilterCriterion[MinThreshold(:posterior_prob, min_posterior,
                              "min posterior (--min-posterior $min_posterior)")], "detection filter")

    # Always apply the count/ratio output filters. Grouping by gene handles an all-novel batch
    # (gene == "" collapses to one bucket) instead of silently skipping the filters.
    add_group_ratio!(collapsed, :count, [:well,:case,:gene], :ratio)
    out_criteria = FilterCriterion[
        MinThreshold(:count, Float64(out_mincount), "min output count (--out-mincount $out_mincount)"),
        MinThreshold(:ratio, out_minratio, "min output ratio (--out-minratio $out_minratio)"),
    ]
    min_heptamer_prob_pre > 0 && push!(out_criteria, MinThreshold(:heptamer_prob_pre, min_heptamer_prob_pre, "min pre-heptamer prob (--min-heptamer-prob-pre $min_heptamer_prob_pre)"))
    min_heptamer_prob_post > 0 && push!(out_criteria, MinThreshold(:heptamer_prob_post, min_heptamer_prob_post, "min post-heptamer prob (--min-heptamer-prob-post $min_heptamer_prob_post)"))
    annotate!(out_criteria, "output filter")

    reason_cols = [:reject_reason, :reject_stage]
    kept = accepted(collapsed)
    outpath = endswith(output, ".gz") ? output : output * ".gz"
    full_output = replace(replace(outpath, r"\.gz$" => ""), r"\.tsv$" => "") * ".full.tsv.gz"
    section("HSMM D detection — summary")
    stage_report("accepted (passed all filters)", nrow(kept), nrow(collapsed))
    report_rejections(collapsed.reject_reason)
    CSV.write(outpath, select(kept, Not(reason_cols)), compress=true, delim='\t')
    @info "Filtered D detections ($(nrow(kept)) rows) saved to $outpath"
    CSV.write(full_output, collapsed, compress=true, delim='\t')
    @info "Full annotated table ($(nrow(collapsed)) candidates + reject reason) saved to $full_output"
    return kept
end

function handle_hsmm(parsed_args)
    @info "HSMM D detection"
    b = parsed_args["discover"]["hsmm"]
    run_hsmm(b["tsv"], b["fasta"], b["output"]; ratio=b["ratio"], mincount=b["mincount"],
        min_gene_len=b["min-gene-len"], max_gene_len=b["max-gene-len"], limit=b["limit"],
        min_posterior=b["min-posterior"], out_mincount=b["out-mincount"], out_minratio=b["out-minratio"],
        min_heptamer_prob_pre=b["min-heptamer-prob-pre"], min_heptamer_prob_post=b["min-heptamer-prob-post"])
end

end # module HSMM
