module Align
    # Core-vs-germline diff using the same Levenshtein alignment as BLAST trim (`core_aln_mismatch`).
    using BioAlignments
    using BioSequences

    export core_edit_distance, core_mismatch_row

    "Edit distance between trimmed core and full DB reference (same as `core_aln_mismatch`)."
    function core_edit_distance(core::AbstractString, reference::AbstractString)
        aln = pairalign(LevenshteinDistance(), LongDNA{4}(reference), LongDNA{4}(core))
        return score(aln)
    end

    """
        core_mismatch_row(core, reference; weight=1.0) -> Vector{Float64}

    Per-position SNP weights along `core` after global alignment to `reference` (the matched
    germline allele). Index `j` follows the trimmed core (`aln_qseq`), not BLAST coordinates.
    """
    function core_mismatch_row(core::AbstractString, reference::AbstractString, weight::Real=1.0)
        Lc = length(core)
        Lc == 0 && return Float64[]
        w = Float64(weight)
        row = zeros(Float64, Lc)
        core == reference && return row
        pa = alignment(pairalign(LevenshteinDistance(), LongDNA{4}(reference), LongDNA{4}(core)))
        gap = DNA_Gap
        ci = 0
        for (a, b) in pa
            b == gap && continue
            ci += 1
            ci > Lc && break
            (a == gap || a != b) && (row[ci] = w)
        end
        return row
    end
end
