# Trim consensus affixes off a BLAST HSP, leaving the gene core. Failure is `absent`;
# success is `Present` of the ungapped core. Prefix vs Suffix extraction is dispatch.

side_label(::Prefix) = "prefix"
side_label(::Suffix) = "suffix"

increment_failures!(stats, ::Prefix) = (stats.prefix_failures += 1)
increment_failures!(stats, ::Suffix) = (stats.suffix_failures += 1)

function record_trim_failure!(stats, side::AffixSide, sseqid, msg)
    push!(stats.failed_genes[sseqid], msg)
    increment_failures!(stats, side)
    return absent
end

function check_affix_quality_warning(affix_length::Int, quality_threshold::Float64)
    if affix_length > 0 && affix_length <= 20 && quality_threshold > 0.5
        @warn "Quality threshold $(round(quality_threshold * 100, digits=1))% might be too strict for short affixes ($affix_length nt). Consider lowering --minquality."
    end
end

function safe_pairalign(seq1::LongDNA{4}, seq2::LongDNA{4}, scoremodel::AffineGapScoreModel)
    (length(seq1) == 0 || length(seq2) == 0) && return absent
    return Present(pairalign(SemiGlobalAlignment(), seq1, seq2, scoremodel))
end

function remove_gaps(query::LongDNA{4})::LongDNA{4}
    return LongDNA{4}(filter(nt -> nt != DNA_Gap, query))
end

"""For prefix: core is everything after the last aligned affix position."""
function extract_core_range(aligned_affix, aligned_query, ::Prefix)
    boundary = findlast(x -> x != DNA_Gap, aligned_affix)
    boundary === nothing && return absent, "No match found"
    start = boundary + 1
    start > length(aligned_query) && return absent, "Query too short after trimming"
    return Present(aligned_query[start:end]), ""
end

"""For suffix: core is everything before the first aligned affix position."""
function extract_core_range(aligned_affix, aligned_query, ::Suffix)
    boundary = findfirst(x -> x != DNA_Gap, aligned_affix)
    boundary === nothing && return absent, "No match found"
    boundary <= 1 && return absent, "Starts too early"
    return Present(aligned_query[1:boundary-1]), ""
end

"""
    AffixTrim(side, scoremodel, min_quality)

Callable trimmer: align `affix` to `query` and keep the opposite side as core.
"""
struct AffixTrim{S<:AffixSide,M}
    side::S
    scoremodel::M
    min_quality::Float64
end

function (t::AffixTrim)(affix::LongDNA{4}, query::LongDNA{4}, stats; sseqid="")
    length(affix) == 0 && return Present(query)
    return finish_affix_align(t, stats, sseqid, safe_pairalign(affix, query, t.scoremodel))
end

finish_affix_align(t::AffixTrim, stats, sseqid, ::Absent) =
    record_trim_failure!(stats, t.side, sseqid, "$(titlecase(side_label(t.side))) alignment failed")

function finish_affix_align(t::AffixTrim, stats, sseqid, aln::Present)
    label = side_label(t.side)
    pairs = collect(alignment(aln.value))
    positions = findall(p -> first(p) != DNA_Gap, pairs)
    isempty(positions) && return record_trim_failure!(stats, t.side, sseqid, "No $label content in alignment")
    matches = sum(first(pairs[i]) == last(pairs[i]) for i in positions)
    quality = matches / length(positions)
    quality < t.min_quality && return record_trim_failure!(stats, t.side, sseqid,
        "Poor $label alignment quality ($(round(quality * 100, digits=1))% match)")
    aligned_affix = first.(pairs)
    aligned_query = last.(pairs)
    core, msg = extract_core_range(aligned_affix, aligned_query, t.side)
    return finish_core_range(t, stats, sseqid, label, core, msg)
end

finish_core_range(t::AffixTrim, stats, sseqid, label, ::Absent, msg) =
    record_trim_failure!(stats, t.side, sseqid, "$(titlecase(label)): $msg")

function finish_core_range(t::AffixTrim, stats, sseqid, label, core::Present, _)
    result = remove_gaps(LongDNA{4}(join(core.value)))
    length(result) == 0 && return record_trim_failure!(stats, t.side, sseqid, "Empty core after $label trimming")
    return Present(result)
end

function trim_sequence(query::LongDNA{4}, prefix::LongDNA{4}, suffix::LongDNA{4}, stats,
    scoremodel::AffineGapScoreModel=AffineGapScoreModel(EDNAFULL, gap_open=-10, gap_extend=-1);
    min_quality=0.75, sseqid="")

    if length(query) == 0
        push!(stats.failed_genes[sseqid], "Empty query sequence")
        return absent
    end
    stats.total_attempts += 1
    pre = AffixTrim(Prefix(), scoremodel, Float64(min_quality))
    suf = AffixTrim(Suffix(), scoremodel, Float64(min_quality))
    return finish_both_trims(pre(prefix, query, stats; sseqid=sseqid), suf, suffix, stats, sseqid)
end

finish_both_trims(::Absent, _, _, _, _) = absent
function finish_both_trims(partial::Present, suf::AffixTrim, suffix, stats, sseqid)
    return finish_core_length(suf(suffix, partial.value, stats; sseqid=sseqid), stats, sseqid)
end

finish_core_length(::Absent, _, _) = absent
function finish_core_length(q::Present, stats, sseqid)
    if length(q.value) < 10
        push!(stats.failed_genes[sseqid], "Core sequence too short ($(length(q.value)) < 10 nt)")
        stats.prefix_failures += 1
        return absent
    end
    push!(stats.successful_trims[sseqid], "Trimmed successfully")
    return q
end

function trim_and_align_sequence(query::String, prefix::String, suffix::String, reference::String, stats;
                                 min_quality=0.75, sseqid="")
    query_dna = BioSequences.LongDNA{4}(query)
    prefix_dna = BioSequences.LongDNA{4}(prefix)
    suffix_dna = BioSequences.LongDNA{4}(suffix)
    trimmed = trim_sequence(query_dna, prefix_dna, suffix_dna, stats, min_quality=min_quality, sseqid=sseqid)
    return core_distance(trimmed, reference)
end

core_distance(::Absent, _) = ("", -1)
function core_distance(trimmed::Present, reference)
    t = replace(String(trimmed.value), '-' => "")
    return t, Align.core_edit_distance(t, reference)
end
