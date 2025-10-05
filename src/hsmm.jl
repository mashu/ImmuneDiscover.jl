module HSMM

using Distributions
using CSV
using DataFrames
using Logging
using ProgressMeter
using Folds
using StringDistances

# Bring in local utilities as submodules to avoid cross-module coupling
include("data.jl")
include("exact.jl")
using .data
using .exact

# -----------------------------------------------------------------------------
# DNA encoding utilities (kept local for module independence)
# -----------------------------------------------------------------------------

@inline function dna_index(c::Char)::Int
    c == 'A' && return 1
    c == 'C' && return 2
    c == 'G' && return 3
    c == 'T' && return 4
    return 0
end

encode_dna(xs::AbstractString) = [dna_index(c) for c in xs]

# -----------------------------------------------------------------------------
# Emission models
# -----------------------------------------------------------------------------

abstract type AbstractEmissionModel end

"""
    IIDLogEmission

Independent and identically distributed emission model with precomputed log-probabilities
for the DNA alphabet order (A, C, G, T).
"""
struct IIDLogEmission <: AbstractEmissionModel
    logp::NTuple{4,Float64}
end

@inline logemit(e::IIDLogEmission, base::Int)::Float64 = e.logp[base]

"""
    PWMLogEmission{L}

Position-specific emission model (PWM) with precomputed log-probabilities.
Rows correspond to positions (1..L), columns to DNA alphabet (A,C,G,T).
"""
struct PWMLogEmission{L} <: AbstractEmissionModel
    logp::Matrix{Float64}  # size (L × 4)
end

@inline function logemit(e::PWMLogEmission{L}, pos::Int, base::Int)::Float64 where {L}
    @inbounds return e.logp[pos, base]
end

# -----------------------------------------------------------------------------
# Duration models
# -----------------------------------------------------------------------------

abstract type AbstractDurationModel end

"""
    DiscreteDuration

Discrete duration distribution over integers in [min_len, max_len] with precomputed logpmf.
"""
struct DiscreteDuration <: AbstractDurationModel
    min_len::Int
    max_len::Int
    logp::Vector{Float64}  # length max_len - min_len + 1
end

@inline function logpmf(d::DiscreteDuration, len::Int)::Float64
    @inbounds return d.logp[len - d.min_len + 1]
end

"""
    discretize_truncated_normal(min_len, max_len, μ, σ)

Build a `DiscreteDuration` by discretizing a truncated normal distribution using midpoint bins.
"""
function discretize_truncated_normal(min_len::Int, max_len::Int, μ::Float64, σ::Float64)::DiscreteDuration
    @assert 0 < σ < Inf "σ must be positive"
    @assert min_len <= max_len "min_len must be ≤ max_len"
    base = truncated(Normal(μ, σ), min_len, max_len)
    n = max_len - min_len + 1
    pmf = Vector{Float64}(undef, n)
    # Use mid-point discretization
    for i in 1:n
        k = min_len + (i - 1)
        left = cdf(base, k - 0.5)
        right = cdf(base, k + 0.5)
        p = max(right - left, eps())
        pmf[i] = p
    end
    s = sum(pmf)
    @assert s > 0.0 "Invalid discretization: zero mass"
    @inbounds for i in eachindex(pmf)
        pmf[i] = log(pmf[i] / s)
    end
    return DiscreteDuration(min_len, max_len, pmf)
end

# -----------------------------------------------------------------------------
# HSMM states
# -----------------------------------------------------------------------------

abstract type AbstractHSMMState end

"""
    FixedMotif{L,E}

Fixed-length emitting state (e.g., nonamer, spacer, heptamer) with position-specific emissions.
"""
struct FixedMotif{L,E<:AbstractEmissionModel} <: AbstractHSMMState
    emission::E
end

Base.length(::FixedMotif{L}) where {L} = L

"""
    DGeneState{E,D}

Variable-length gene state with IID emissions and a duration distribution.
"""
struct DGeneState{E<:AbstractEmissionModel,D<:AbstractDurationModel} <: AbstractHSMMState
    emission::E
    duration::D
end

@inline minduration(s::DGeneState) = s.duration.min_len
@inline maxduration(s::DGeneState) = s.duration.max_len

# -----------------------------------------------------------------------------
# Complete model: RSS (5') → D gene → RSS (3')
# Pattern: nonamer-spacer-heptamer-dgene-heptamer-spacer-nonamer
# -----------------------------------------------------------------------------

struct DGeneRSSHSMM{E9<:AbstractEmissionModel,E12<:AbstractEmissionModel,E7<:AbstractEmissionModel,EG<:AbstractEmissionModel,D<:AbstractDurationModel}
    pre_nonamer::FixedMotif{9,E9}
    pre_spacer::FixedMotif{12,E12}
    pre_heptamer::FixedMotif{7,E7}
    gene::DGeneState{EG,D}
    post_heptamer::FixedMotif{7,E7}
    post_spacer::FixedMotif{12,E12}
    post_nonamer::FixedMotif{9,E9}
end

@inline mingene(m::DGeneRSSHSMM) = minduration(m.gene)
@inline maxgene(m::DGeneRSSHSMM) = maxduration(m.gene)

# -----------------------------------------------------------------------------
# Estimation utilities from labeled tuples
# -----------------------------------------------------------------------------

"""
    build_pwm_log(sequences; pseudocount=0.01)

Build log-PWM (rows × 4) from aligned sequences of the same length using a Dirichlet
pseudocount. Returns a matrix of log-probabilities.
"""
function build_pwm_log(sequences::Vector{<:AbstractString}; pseudocount::Float64=0.01)
    @assert !isempty(sequences) "Need at least one sequence"
    L = length(sequences[1])
    @assert all(length(s) == L for s in sequences) "All sequences must be the same length"
    counts = fill(pseudocount, L, 4)
    for seq in sequences
        @assert length(seq) == L
        for (pos, ch) in enumerate(seq)
            idx = dna_index(ch)
            if idx > 0
                @inbounds counts[pos, idx] += 1.0
            end
        end
    end
    # Normalize rows to probabilities and take logs
    logp = similar(counts)
    @inbounds for i in 1:L
        s = sum(counts[i, :])
        invs = 1.0 / s
        for b in 1:4
            logp[i, b] = log(counts[i, b] * invs)
        end
    end
    return logp
end

"""
    estimate_iid_log_emission(sequences; pseudocount=0.01)

Estimate IID emission log-probabilities from sequences by counting base frequencies
with Dirichlet pseudocounts.
"""
function estimate_iid_log_emission(sequences::Vector{<:AbstractString}; pseudocount::Float64=0.01)::IIDLogEmission
    counts = (pseudocount, pseudocount, pseudocount, pseudocount)
    A = Float64(counts[1]); C = Float64(counts[2]); G = Float64(counts[3]); T = Float64(counts[4])
    for s in sequences
        for ch in s
            idx = dna_index(ch)
            if idx == 1
                A += 1.0
            elseif idx == 2
                C += 1.0
            elseif idx == 3
                G += 1.0
            elseif idx == 4
                T += 1.0
            end
        end
    end
    total = A + C + G + T
    invt = 1.0 / total
    return IIDLogEmission((log(A * invt), log(C * invt), log(G * invt), log(T * invt)))
end

"""
    empirical_duration(sequences, min_len, max_len; pseudocount=0.1)

Construct a `DiscreteDuration` model from observed lengths, clamped to [min_len,max_len].
"""
function empirical_duration(dgene_sequences::Vector{<:AbstractString}, min_len::Int, max_len::Int; pseudocount::Float64=0.1)::DiscreteDuration
    @assert min_len <= max_len
    n = max_len - min_len + 1
    counts = fill(pseudocount, n)
    for s in dgene_sequences
        L = length(s)
        if min_len <= L <= max_len
            @inbounds counts[L - min_len + 1] += 1.0
        end
    end
    total = sum(counts)
    invt = 1.0 / total
    @inbounds for i in eachindex(counts)
        counts[i] = log(counts[i] * invt)
    end
    return DiscreteDuration(min_len, max_len, counts)
end

"""
    fit_dgene_rss_hsmm(data, min_gene, max_gene; kwargs...)

Fit model parameters from labeled 7-tuples: (pre_nonamer, pre_spacer, pre_heptamer,
dgene, post_heptamer, post_spacer, post_nonamer).

Options:
- `pseudocount::Float64=0.01`: Dirichlet pseudocount for emissions
- `duration=:empirical`: build discrete duration from observed lengths
- `duration_mean::Float64`, `duration_std::Float64`: if provided and `duration=:normal`,
   discretize a truncated normal in [min_gene, max_gene]
"""
function fit_dgene_rss_hsmm(
    data::Vector{<:NTuple{7,AbstractString}},
    min_gene::Int,
    max_gene::Int;
    pseudocount::Float64=0.01,
    duration::Symbol=:empirical,
    duration_mean::Float64=(min_gene + max_gene) / 2,
    duration_std::Float64= max(1.0, (max_gene - min_gene) / 6),
)
    # Split tuples into components
    pre_nonamers = String[]
    pre_spacers = String[]
    pre_heptamers = String[]
    dgenes = String[]
    post_heptamers = String[]
    post_spacers = String[]
    post_nonamers = String[]

    sizehint!(pre_nonamers, length(data))
    sizehint!(pre_spacers, length(data))
    sizehint!(pre_heptamers, length(data))
    sizehint!(dgenes, length(data))
    sizehint!(post_heptamers, length(data))
    sizehint!(post_spacers, length(data))
    sizehint!(post_nonamers, length(data))

    for tup in data
        push!(pre_nonamers, tup[1])
        push!(pre_spacers, tup[2])
        push!(pre_heptamers, tup[3])
        push!(dgenes, tup[4])
        push!(post_heptamers, tup[5])
        push!(post_spacers, tup[6])
        push!(post_nonamers, tup[7])
    end

    @assert all(length(s) == 9 for s in pre_nonamers) "pre_nonamer must be length 9"
    @assert all(length(s) == 12 for s in pre_spacers) "pre_spacer must be length 12"
    @assert all(length(s) == 7 for s in pre_heptamers) "pre_heptamer must be length 7"
    @assert all(length(s) == 7 for s in post_heptamers) "post_heptamer must be length 7"
    @assert all(length(s) == 12 for s in post_spacers) "post_spacer must be length 12"
    @assert all(length(s) == 9 for s in post_nonamers) "post_nonamer must be length 9"

    # Emissions (log-space)
    pre_nonamer_log = build_pwm_log(pre_nonamers; pseudocount)
    pre_spacer_log = build_pwm_log(pre_spacers; pseudocount)
    pre_heptamer_log = build_pwm_log(pre_heptamers; pseudocount)
    post_heptamer_log = build_pwm_log(post_heptamers; pseudocount)
    post_spacer_log = build_pwm_log(post_spacers; pseudocount)
    post_nonamer_log = build_pwm_log(post_nonamers; pseudocount)

    gene_emission = estimate_iid_log_emission(dgenes; pseudocount)

    dur_model = duration === :empirical ?
        empirical_duration(dgenes, min_gene, max_gene; pseudocount=max(pseudocount, 1e-6)) :
        discretize_truncated_normal(min_gene, max_gene, duration_mean, duration_std)

    return DGeneRSSHSMM(
        FixedMotif{9,PWMLogEmission{9}}(PWMLogEmission{9}(pre_nonamer_log)),
        FixedMotif{12,PWMLogEmission{12}}(PWMLogEmission{12}(pre_spacer_log)),
        FixedMotif{7,PWMLogEmission{7}}(PWMLogEmission{7}(pre_heptamer_log)),
        DGeneState{IIDLogEmission,DiscreteDuration}(gene_emission, dur_model),
        FixedMotif{7,PWMLogEmission{7}}(PWMLogEmission{7}(post_heptamer_log)),
        FixedMotif{12,PWMLogEmission{12}}(PWMLogEmission{12}(post_spacer_log)),
        FixedMotif{9,PWMLogEmission{9}}(PWMLogEmission{9}(post_nonamer_log)),
    )
end

# -----------------------------------------------------------------------------
# Scoring and Viterbi-like local detection (no insertions/deletions within motifs)
# -----------------------------------------------------------------------------

@inline function score_motif(obs::Vector{Int}, start_pos::Int, s::FixedMotif{L,PWMLogEmission{L}})::Float64 where {L}
    total = 0.0
    @inbounds for i in 1:L
        base = obs[start_pos + i - 1]
        if base == 0
            return -Inf
        end
        total += logemit(s.emission, i, base)
    end
    return total
end

@inline function score_gene_segment(obs::Vector{Int}, start_pos::Int, d::Int, s::DGeneState{IIDLogEmission})::Float64
    @inbounds begin
        logpA, logpC, logpG, logpT = s.emission.logp
        total = logpmf(s.duration, d)
        for i in 0:(d-1)
            base = obs[start_pos + i]
            if base == 1
                total += logpA
            elseif base == 2
                total += logpC
            elseif base == 3
                total += logpG
            elseif base == 4
                total += logpT
            else
                return -Inf
            end
        end
        return total
    end
end

"""
    extract_dgene(sequence, model)

Find the best-scoring occurrence of the RSS-D-RSS pattern inside `sequence`.
Returns a named tuple with positions (1-based, inclusive) and the extracted D gene.
"""
function extract_dgene(sequence::AbstractString, model::DGeneRSSHSMM)
    obs = encode_dna(sequence)
    T = length(obs)

    k_pre = 9 + 12 + 7
    k_post = 7 + 12 + 9

    best_log = -Inf
    best = (prefix_start=0, prefix_end=0, gene_start=0, gene_end=0, suffix_start=0, suffix_end=0)

    min_g = mingene(model)
    max_g = maxgene(model)

    @inbounds for prefix_start in 1:(T - (k_pre + min_g + k_post) + 1)
        # Pre RSS blocks
        n9_start = prefix_start
        s12_start = n9_start + 9
        h7_start = s12_start + 12
        prefix_end = h7_start + 7 - 1

        lp = score_motif(obs, n9_start, model.pre_nonamer)
        lp += score_motif(obs, s12_start, model.pre_spacer)
        lp += score_motif(obs, h7_start, model.pre_heptamer)
        if lp == -Inf
            continue
        end

        gene_start = prefix_end + 1

        # Try all durations
        for d in min_g:max_g
            gene_end = gene_start + d - 1
            suffix_start = gene_end + 1
            suffix_end = suffix_start + k_post - 1
            if suffix_end > T
                break
            end

            lg = score_gene_segment(obs, gene_start, d, model.gene)
            if lg == -Inf
                continue
            end

            # Post RSS blocks
            h7p_start = suffix_start
            s12p_start = h7p_start + 7
            n9p_start = s12p_start + 12

            ls = score_motif(obs, h7p_start, model.post_heptamer)
            ls += score_motif(obs, s12p_start, model.post_spacer)
            ls += score_motif(obs, n9p_start, model.post_nonamer)
            if ls == -Inf
                continue
            end

            total = lp + lg + ls
            if total > best_log
                best_log = total
                best = (
                    prefix_start=n9_start,
                    prefix_end=prefix_end,
                    gene_start=gene_start,
                    gene_end=gene_end,
                    suffix_start=suffix_start,
                    suffix_end=suffix_end,
                )
            end
        end
    end

    if best_log == -Inf
        return (gene_seq = "", prefix_start = 0, prefix_end = 0, gene_start = 0, gene_end = 0,
                suffix_start = 0, suffix_end = 0, log_prob = -Inf)
    end

    gene_seq = sequence[best.gene_start:best.gene_end]
    return (gene_seq = gene_seq,
            prefix_start = best.prefix_start,
            prefix_end = best.prefix_end,
            gene_start = best.gene_start,
            gene_end = best.gene_end,
            suffix_start = best.suffix_start,
            suffix_end = best.suffix_end,
            log_prob = best_log)
end

# -----------------------------------------------------------------------------
# Convenience constructors
# -----------------------------------------------------------------------------

"""
    from_pwm_and_duration(pre_nonamer_pwm, pre_spacer_pwm, pre_heptamer_pwm,
                          post_heptamer_pwm, post_spacer_pwm, post_nonamer_pwm,
                          gene_emission_logp, duration)

Build a `DGeneRSSHSMM` from precomputed log-PWMs, gene IID log-emission and duration model.
"""
function from_pwm_and_duration(
    pre_nonamer_log::Matrix{Float64},
    pre_spacer_log::Matrix{Float64},
    pre_heptamer_log::Matrix{Float64},
    post_heptamer_log::Matrix{Float64},
    post_spacer_log::Matrix{Float64},
    post_nonamer_log::Matrix{Float64},
    gene_emission_logp::NTuple{4,Float64},
    duration::DiscreteDuration,
)
    @assert size(pre_nonamer_log) == (9, 4)
    @assert size(pre_spacer_log) == (12, 4)
    @assert size(pre_heptamer_log) == (7, 4)
    @assert size(post_heptamer_log) == (7, 4)
    @assert size(post_spacer_log) == (12, 4)
    @assert size(post_nonamer_log) == (9, 4)
    return DGeneRSSHSMM(
        FixedMotif{9,PWMLogEmission{9}}(PWMLogEmission{9}(pre_nonamer_log)),
        FixedMotif{12,PWMLogEmission{12}}(PWMLogEmission{12}(pre_spacer_log)),
        FixedMotif{7,PWMLogEmission{7}}(PWMLogEmission{7}(pre_heptamer_log)),
        DGeneState{IIDLogEmission,DiscreteDuration}(IIDLogEmission(gene_emission_logp), duration),
        FixedMotif{7,PWMLogEmission{7}}(PWMLogEmission{7}(post_heptamer_log)),
        FixedMotif{12,PWMLogEmission{12}}(PWMLogEmission{12}(post_spacer_log)),
        FixedMotif{9,PWMLogEmission{9}}(PWMLogEmission{9}(post_nonamer_log)),
    )
end

export AbstractEmissionModel, IIDLogEmission, PWMLogEmission, AbstractDurationModel, DiscreteDuration,
       AbstractHSMMState, FixedMotif, DGeneState, DGeneRSSHSMM,
       build_pwm_log, estimate_iid_log_emission, empirical_duration,
       discretize_truncated_normal, fit_dgene_rss_hsmm, extract_dgene,
       from_pwm_and_duration, encode_dna

@inline function logaddexp(a::Float64, b::Float64)::Float64
    if a == -Inf
        return b
    elseif b == -Inf
        return a
    end
    m = max(a, b)
    return m + log1p(exp(min(a, b) - m))
end

@inline function pwm_logprob(seq::AbstractString, pwm::PWMLogEmission{L})::Float64 where {L}
    @assert length(seq) == L
    total = 0.0
    @inbounds for (pos, ch) in enumerate(seq)
        idx = dna_index(ch)
        if idx == 0
            return -Inf
        end
        total += logemit(pwm, pos, idx)
    end
    return total
end

"""
    scan_best_and_total(sequence::AbstractString, model::DGeneRSSHSMM)

Compute the best-scoring placement of the RSS–D–RSS pattern and the total
log-probability mass (log-sum-exp) across all valid placements.
"""
function scan_best_and_total(sequence::AbstractString, model::DGeneRSSHSMM)
    obs = encode_dna(sequence)
    T = length(obs)

    k_pre = 9 + 12 + 7
    k_post = 7 + 12 + 9

    best_log = -Inf
    best = (prefix_start=0, prefix_end=0, gene_start=0, gene_end=0, suffix_start=0, suffix_end=0)
    total_log = -Inf

    min_g = mingene(model)
    max_g = maxgene(model)

    @inbounds for prefix_start in 1:(T - (k_pre + min_g + k_post) + 1)
        # Pre RSS blocks
        n9_start = prefix_start
        s12_start = n9_start + 9
        h7_start = s12_start + 12
        prefix_end = h7_start + 7 - 1

        lp = score_motif(obs, n9_start, model.pre_nonamer)
        lp += score_motif(obs, s12_start, model.pre_spacer)
        lp += score_motif(obs, h7_start, model.pre_heptamer)
        if lp == -Inf
            continue
        end

        gene_start = prefix_end + 1

        # Try all durations
        for d in min_g:max_g
            gene_end = gene_start + d - 1
            suffix_start = gene_end + 1
            suffix_end = suffix_start + k_post - 1
            if suffix_end > T
                break
            end

            lg = score_gene_segment(obs, gene_start, d, model.gene)
            if lg == -Inf
                continue
            end

            # Post RSS blocks
            h7p_start = suffix_start
            s12p_start = h7p_start + 7
            n9p_start = s12p_start + 12

            ls = score_motif(obs, h7p_start, model.post_heptamer)
            ls += score_motif(obs, s12p_start, model.post_spacer)
            ls += score_motif(obs, n9p_start, model.post_nonamer)
            if ls == -Inf
                continue
            end

            total = lp + lg + ls
            total_log = logaddexp(total_log, total)
            if total > best_log
                best_log = total
                best = (
                    prefix_start=n9_start,
                    prefix_end=prefix_end,
                    gene_start=gene_start,
                    gene_end=gene_end,
                    suffix_start=suffix_start,
                    suffix_end=suffix_end,
                )
            end
        end
    end

    if best_log == -Inf
        return (gene_seq = "", prefix_start = 0, prefix_end = 0, gene_start = 0, gene_end = 0,
                suffix_start = 0, suffix_end = 0, log_path_prob = -Inf, log_total_prob = -Inf,
                posterior_prob = 0.0)
    end

    gene_seq = sequence[best.gene_start:best.gene_end]
    post = best_log - total_log
    posterior = isfinite(post) ? exp(post) : 0.0
    return (gene_seq = gene_seq,
            prefix_start = best.prefix_start,
            prefix_end = best.prefix_end,
            gene_start = best.gene_start,
            gene_end = best.gene_end,
            suffix_start = best.suffix_start,
            suffix_end = best.suffix_end,
            log_path_prob = best_log,
            log_total_prob = total_log,
            posterior_prob = posterior)
end

# -----------------------------------------------------------------------------
# HSMM pipeline: known D selection → flank training → detection
# -----------------------------------------------------------------------------

"""
    run_hsmm(tsv::String, fasta::String, output::String;
             ratio::Float64=0.2, mincount::Int=5,
             min_gene_len::Int=0, max_gene_len::Int=0,
             limit::Int=0)

Detect D genes using an HSMM trained on recombination signal sequences (RSS) flanks.

Steps:
- Search for known D alleles per donor (case) and extract RSS flanks on both sides.
- Apply per-donor, per-gene allelic ratio filter ≥ `ratio` and retain at most 2 alleles
  for the same gene per donor, based only on allele sequence (ignoring flanks).
- Train an HSMM on all flanks from the retained known D alleles.
- Apply the model to all reads to extract candidate D segments with flanks; mark if present in DB.

Writes a TSV.GZ with detected candidates and returns the resulting DataFrame.
"""
function run_hsmm(tsv::String, fasta::String, output::String;
                  ratio::Float64=0.2, mincount::Int=5,
                  min_gene_len::Int=0, max_gene_len::Int=0,
                  limit::Int=0, min_posterior::Float64=0.7,
                  out_mincount::Int=10, out_minratio::Float64=0.2,
                  min_heptamer_prob_pre::Float64=0.0,
                  min_heptamer_prob_post::Float64=0.0)
    # Load inputs
    @info "Loading demultiplex: $tsv"
    table = limit > 0 ? data.load_demultiplex(tsv, limit=limit) : data.load_demultiplex(tsv)
    @info "Loaded $(nrow(table)) rows"

    @info "Loading D allele FASTA: $fasta"
    db = data.load_fasta(fasta, validate=false)
    db_seq_lookup = Dict{String,String}((seq => name) for (name, seq) in db)
    db_names = first.(db)
    db_seqs = last.(db)

    # Known D exact search to collect flanks; do not filter by ratio here
    @info "Searching known D alleles and extracting RSS flanks"
    known_df = exact.exact_search(table, db, "D"; mincount=mincount, minratio=0.0, N=1000)
    if nrow(known_df) == 0
        @warn "No exact D matches found; HSMM training cannot proceed"
        return DataFrame()
    end

    # Aggregate by (case, gene, sequence) across wells and compute per-donor (case) ratios
    agg = combine(groupby(known_df, [:case, :gene, :sequence, :db_name]),
                  :count => sum => :case_count)
    transform!(groupby(agg, [:case, :gene]), :case_count => (x -> x ./ maximum(x)) => :case_ratio)

    # Apply allelic ratio filter (per donor/gene) and mincount on aggregated counts
    filter!(r -> (r.case_count >= mincount) && (r.case_ratio >= ratio), agg)
    if nrow(agg) == 0
        @warn "No alleles passed ratio >= $(ratio) and mincount >= $(mincount); HSMM training cannot proceed"
        return DataFrame()
    end

    kept = Vector{Tuple{String,String,String}}()  # (case, gene, sequence)
    for grp in groupby(agg, [:case, :gene])
        # Sort by aggregated count then aggregated ratio, descending
        ord = sort(grp, [:case_count, :case_ratio], rev=[true, true])
        n = min(2, nrow(ord))
        for i in 1:n
            push!(kept, (ord[i, :case], ord[i, :gene], ord[i, :sequence]))
        end
    end
    kept_set = Set(kept)

    # Filter known_df rows to training set sequences
    train_df = filter(r -> (r.case, r.gene, r.sequence) in kept_set, known_df)

    # Build labeled tuples for PWM/HSMM training; keep only complete-length flanks
    function complete_flanks(row)
        try
            return length(row.pre_nonamer) == 9 &&
                   length(row.pre_spacer) == 12 &&
                   length(row.pre_heptamer) == 7 &&
                   length(row.post_heptamer) == 7 &&
                   length(row.post_spacer) == 12 &&
                   length(row.post_nonamer) == 9
        catch
            return false
        end
    end

    train_df = filter(complete_flanks, train_df)
    if nrow(train_df) == 0
        @warn "No complete-length RSS flanks available for training"
        return DataFrame()
    end

    tuples = [(row.pre_nonamer, row.pre_spacer, row.pre_heptamer,
               row.sequence,
               row.post_heptamer, row.post_spacer, row.post_nonamer) for row in eachrow(train_df)]

    # Determine gene length bounds if not provided
    lengths = length.(train_df.sequence)
    minL = min_gene_len > 0 ? min_gene_len : minimum(lengths)
    maxL = max_gene_len > 0 ? max_gene_len : maximum(lengths)
    @info "Training HSMM with D length range [$minL, $maxL] from $(length(tuples)) examples"

    model = fit_dgene_rss_hsmm(tuples, minL, maxL)

    # Apply model to all reads to detect candidate Ds (parallel)
    pre_blocks = (n9=9, s12=12, h7=7)
    post_blocks = (h7=7, s12=12, n9=9)

    @info "Scanning reads with HSMM (parallel with chunked progress)"
    N = nrow(table)
    chunk = max(10_000, min(100_000, Int(cld(N, 50))))
    starts = collect(1:chunk:N)
    prog = Progress(length(starts), desc="Scanning chunks")
    buffer = NamedTuple[]
    for s in starts
        e = min(s + chunk - 1, N)
        sub = view(table, s:e, :)
        chunk_results = Folds.map(eachrow(sub)) do row
            seq = row.genomic_sequence
            det = scan_best_and_total(seq, model)
            if det.gene_seq == ""
                return nothing
            end
            # Reconstruct flanks from positions
            ps = det.prefix_start
            pe = det.prefix_end
            gs = det.gene_start
            ge = det.gene_end
            ss = det.suffix_start
            # Pre RSS
            pre_nonamer = seq[ps:ps + pre_blocks.n9 - 1]
            pre_spacer = seq[ps + pre_blocks.n9 : ps + pre_blocks.n9 + pre_blocks.s12 - 1]
            pre_heptamer = seq[ps + pre_blocks.n9 + pre_blocks.s12 : pe]
            # Post RSS
            post_heptamer = seq[ss : ss + post_blocks.h7 - 1]
            post_spacer = seq[ss + post_blocks.h7 : ss + post_blocks.h7 + post_blocks.s12 - 1]
            post_nonamer = seq[ss + post_blocks.h7 + post_blocks.s12 : ss + post_blocks.h7 + post_blocks.s12 + post_blocks.n9 - 1]

            dseq = det.gene_seq
            isin = haskey(db_seq_lookup, dseq)
            db_name = isin ? db_seq_lookup[dseq] : "Novel"

            # Nearest DB allele by Levenshtein distance
            best_dist = typemax(Int)
            best_call = ""
            if !isempty(db_seqs)
                for (i, refseq) in enumerate(db_seqs)
                    dist = evaluate(Levenshtein(), dseq, refseq)
                    if dist < best_dist
                        best_dist = dist
                        best_call = db_names[i]
                        if best_dist == 0
                            break
                        end
                    end
                end
            end

        return (
            well = string(row.well),
            case = string(row.case),
            sequence = dseq,
            pre_nonamer = pre_nonamer,
            pre_spacer = pre_spacer,
            pre_heptamer = pre_heptamer,
            post_heptamer = post_heptamer,
            post_spacer = post_spacer,
            post_nonamer = post_nonamer,
            heptamer_logp_pre = pwm_logprob(pre_heptamer, model.pre_heptamer.emission),
            heptamer_logp_post = pwm_logprob(post_heptamer, model.post_heptamer.emission),
            log_path_prob = det.log_path_prob,
            log_total_prob = det.log_total_prob,
            posterior_prob = det.posterior_prob,
            isin_db = isin,
            db_name = db_name,
            nearest_db = best_call,
            nearest_db_dist = best_dist,
        )
        end
        # Append valid results from this chunk
        for r in chunk_results
            if r !== nothing
                push!(buffer, r)
            end
        end
        next!(prog)
    end
    finish!(prog)

    valid_results = buffer
    if isempty(valid_results)
        @warn "No D segments detected by HSMM"
        return DataFrame()
    end

    res_df = DataFrame(valid_results)
    # Filter detections by minimum posterior probability
    before = nrow(res_df)
    filter!(x -> x.posterior_prob >= min_posterior, res_df)
    @info "Posterior filter >= $(min_posterior): kept $(nrow(res_df)) / $(before) detections"

    # Collapse to counts per (well, case, sequence) and keep one set of flanks
    grouped = groupby(res_df, [:well, :case, :sequence])
    collapsed = combine(grouped, nrow => :count,
        :pre_nonamer => first => :pre_nonamer,
        :pre_spacer => first => :pre_spacer,
        :pre_heptamer => first => :pre_heptamer,
        :post_heptamer => first => :post_heptamer,
        :post_spacer => first => :post_spacer,
        :post_nonamer => first => :post_nonamer,
        :heptamer_logp_pre => first => :heptamer_logp_pre,
        :heptamer_logp_post => first => :heptamer_logp_post,
        :log_path_prob => first => :log_path_prob,
        :log_total_prob => first => :log_total_prob,
        :posterior_prob => first => :posterior_prob,
        :isin_db => any => :isin_db,
        :db_name => first => :db_name,
        :nearest_db => first => :nearest_db,
        :nearest_db_dist => minimum => :nearest_db_dist)

    # Allele naming based on best match: keep original name if known; otherwise add _S hash from sequence
    collapsed[!, :allele_name] = map(eachrow(collapsed)) do row
        is_known = (row.nearest_db_dist == 0) || row.isin_db
        base_name = (row.nearest_db != "") ? String(row.nearest_db) : String(row.db_name)
        is_known ? base_name : unique_name(base_name, String(row.sequence))
    end

    # Convert heptamer log-probabilities to probabilities for user-friendly thresholds
    collapsed[!, :heptamer_prob_pre] = map(x -> isfinite(x) ? exp(x) : 0.0, collapsed.heptamer_logp_pre)
    collapsed[!, :heptamer_prob_post] = map(x -> isfinite(x) ? exp(x) : 0.0, collapsed.heptamer_logp_post)

    # Output-level mincount and allelic ratio filters per donor and gene
    # First, derive gene from nearest_db if present else empty; fallback: keep empty gene label
    collapsed[!, :gene] = map(row -> begin
        name = String(row.nearest_db)
        if name == ""
            ""
        else
            first(split(name, '*'))
        end
    end, eachrow(collapsed))

    # Compute ratio within (well, case, gene)
    if any(x->x != "", collapsed.gene)
        transform!(groupby(collapsed, [:well, :case, :gene]), :count => (x->x./maximum(x)) => :ratio)
        before_out = nrow(collapsed)
        filter!(x -> x.count >= out_mincount, collapsed)
        filter!(x -> x.ratio >= out_minratio, collapsed)
        @info "Output filters: kept $(nrow(collapsed)) / $(before_out) with count >= $(out_mincount) and ratio >= $(out_minratio)"
    else
        @info "Skipping allelic ratio filter: gene names unavailable (nearest_db empty)"
    end

    # Heptamer probability-based filtering if requested
    before_hep = nrow(collapsed)
    filter!(x -> (x.heptamer_prob_pre >= min_heptamer_prob_pre) & (x.heptamer_prob_post >= min_heptamer_prob_post), collapsed)
    if (min_heptamer_prob_pre > 0.0) || (min_heptamer_prob_post > 0.0)
        @info "Heptamer PWM prob filter: kept $(nrow(collapsed)) / $(before_hep)"
    end

    # Save
    outpath = endswith(output, ".gz") ? output : output * ".gz"
    CSV.write(outpath, collapsed, compress=true, delim='\t')
    @info "HSMM detection results saved to $outpath ($(nrow(collapsed)) rows)"

    return collapsed
end

end # module HSMM



