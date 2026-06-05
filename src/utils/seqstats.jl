module SeqStats
    # Pure sequence-quality / variability statistics used to separate genuine novel alleles
    # from artifacts (composition, cross-read agreement, diversity). No external deps.

    export gc_content, max_homopolymer, n_content, shannon_entropy, consensus_fraction

    "Fraction of G/C bases (0 for empty)."
    function gc_content(seq::AbstractString)
        isempty(seq) && return 0.0
        gc = count(c -> c == 'G' || c == 'C' || c == 'g' || c == 'c', seq)
        return gc / length(seq)
    end

    "Length of the longest run of a single repeated base (0 for empty)."
    function max_homopolymer(seq::AbstractString)
        isempty(seq) && return 0
        best = 1
        run = 1
        prev = first(seq)
        for c in Iterators.drop(seq, 1)
            if c == prev
                run += 1
                best = max(best, run)
            else
                run = 1
                prev = c
            end
        end
        return best
    end

    "Fraction of N bases (0 for empty)."
    function n_content(seq::AbstractString)
        isempty(seq) && return 0.0
        return count(c -> c == 'N' || c == 'n', seq) / length(seq)
    end

    _require_equal_length(seqs) =
        all(s -> length(s) == length(first(seqs)), seqs) ||
            throw(ArgumentError("sequences must be equal length"))

    """
        shannon_entropy(seqs) -> Float64

    Mean per-position Shannon entropy (bits) over a set of equal-length sequences. 0 when all
    sequences are identical (a clean consensus); higher when positions disagree (noisy cluster).
    """
    function shannon_entropy(seqs)
        isempty(seqs) && return 0.0
        _require_equal_length(seqs)
        L = length(first(seqs))
        L == 0 && return 0.0
        n = length(seqs)
        total = 0.0
        counts = Dict{Char,Int}()
        for j in 1:L
            empty!(counts)
            for s in seqs
                c = s[j]
                counts[c] = get(counts, c, 0) + 1
            end
            h = 0.0
            for cnt in values(counts)
                p = cnt / n
                h -= p * log2(p)
            end
            total += h
        end
        return total / L
    end

    """
        consensus_fraction(seqs) -> Float64

    Mean fraction of sequences matching the per-position consensus base (1.0 when identical).
    """
    function consensus_fraction(seqs)
        isempty(seqs) && return 1.0
        _require_equal_length(seqs)
        L = length(first(seqs))
        L == 0 && return 1.0
        n = length(seqs)
        agree = 0
        counts = Dict{Char,Int}()
        for j in 1:L
            empty!(counts)
            for s in seqs
                c = s[j]
                counts[c] = get(counts, c, 0) + 1
            end
            agree += maximum(values(counts))
        end
        return agree / (L * n)
    end
end
