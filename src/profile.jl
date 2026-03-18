module profile
    const dna2ind = Dict{Char,Int}('A'=>1,'C'=>2,'G'=>3,'T'=>4)
    const ind2dna = Dict{Int,Char}(1=>'A',2=>'C',3=>'G',4=>'T')

    """
        counts(motifs)

    Compute counts given motifs of equal length
    """
    function counts(motifs)
        cols = length(motifs[1])
        rows = length(motifs)
        c = zeros(Int64, length(dna2ind), cols)
        for col in 1:cols
            for row in 1:rows
                ind = get(dna2ind,motifs[row][col],0)
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

    """
        motif_prob(motif, prof)

    Calculate the product of all probabilties from nucleotides in a motif given the profile.
    Returns probability of a motif.
    """
    function motif_prob(motif, prof)
        accu = 0.0
        for i in eachindex(motif)
            if motif[i] == '-'
                return 0.0  # Gaps have zero probability (Float64 for type stability)
            end
            ind = get(dna2ind, motif[i], 0)
            if ind > 0
                accu += log2(prof[ind, i])
            end
        end
        return 2^(accu)  # Summing of logs more numerically stable
    end
end
