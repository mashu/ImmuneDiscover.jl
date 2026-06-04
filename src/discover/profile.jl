module Profile
    # Branch-based nucleotide → index (1..4, 0 for gaps/unknown). Allocation-free and
    # type-stable, matching HSMM.dna_index; avoids a Dict lookup in the counting loop.
    @inline function dna2ind(c::Char)::Int
        c == 'A' && return 1
        c == 'C' && return 2
        c == 'G' && return 3
        c == 'T' && return 4
        return 0
    end

    """
        counts(motifs)

    Compute counts given motifs of equal length
    """
    function counts(motifs)
        cols = length(motifs[1])
        rows = length(motifs)
        c = zeros(Int64, 4, cols)
        for col in 1:cols
            for row in 1:rows
                ind = dna2ind(motifs[row][col])
                if ind > 0  # Gaps have zero count and probability
                    c[ind,col] += 1
                end
            end
        end
        return c
    end

    """
        motif_profile(motifs)

    Compute profile from motifs
    """
    function motif_profile(motifs)
        c = counts(motifs) .+ 0.0001
        return c ./ sum(c, dims=1)
    end
end
