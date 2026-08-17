module Profile
    using ..DNA: dna_index

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
                ind = dna_index(motifs[row][col])
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
