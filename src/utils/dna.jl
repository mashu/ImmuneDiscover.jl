module DNA
    # Nucleotide → 1..4 index used by motif profiles and the HSMM. One definition so Profile
    # and HSMM cannot drift.

    export dna_index, encode_dna

    @inline function dna_index(c::Char)::Int
        c == 'A' && return 1
        c == 'C' && return 2
        c == 'G' && return 3
        c == 'T' && return 4
        return 0
    end

    encode_dna(xs::AbstractString) = [dna_index(c) for c in xs]
end
