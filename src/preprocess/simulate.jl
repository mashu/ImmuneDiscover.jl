module Simulate
    using FASTX
    using CSV
    using DataFrames
    using Random
    using MD5

    using ..Data: unique_name, sequence_hash
    using ..Exact: GeneType, VGene, DGene, JGene

    const NUCLEOTIDES = ['A', 'C', 'G', 'T']
    const RSS_SIGNAL = "CACAGTG"
    const LEADER_SIGNAL = "GTTTTTGT"
    const BARCODE_LENGTH = 10

    # Canonical recombination signal sequence (RSS) building blocks for synthetic
    # V/D/J reads. Real RSS = heptamer(7) – spacer(12 or 23) – nonamer(9). These let
    # the simulator emit reads with the flank architecture that the exact/heptamer
    # search and the D-gene HSMM expect, so detection can be tested end-to-end.
    const HEPTAMER = RSS_SIGNAL          # 7 nt conserved heptamer (CACAGTG)
    const NONAMER  = "ACAAAAACC"         # 9 nt conserved nonamer
    const SPACER12 = 12                  # 12 nt spacer (D genes; one V/J side)
    const SPACER23 = 23                  # 23 nt spacer (V/J genes)

    export generate_fasta_with_mutations, unique_name, sequence_hash
    export rss_prefix, rss_suffix, assemble_read, simulate_d_read,
           v_end_variant, short_d_variant, decoy_read, invalid_d_read

    function hamming_distance(s1::String, s2::String)
        @assert length(s1) == length(s2)
        return sum(c1 != c2 for (c1, c2) in zip(s1, s2))
    end

    function generate_distant_barcode(existing_barcodes::Vector{String}, min_distance::Int=4)
        max_attempts = 1000
        for _ in 1:max_attempts
            candidate = join(rand(NUCLEOTIDES, BARCODE_LENGTH))
            if isempty(existing_barcodes) ||
               all(hamming_distance(candidate, existing) >= min_distance for existing in existing_barcodes)
                return candidate
            end
        end
        error("Could not generate barcode with required minimum distance after $max_attempts attempts")
    end

    function generate_barcode_indices(output_file::String="indices.tsv")
        forward_barcodes = String[]
        for _ in 1:2
            push!(forward_barcodes, generate_distant_barcode(forward_barcodes))
        end
        reverse_barcodes = String[]
        all_existing = vcat(forward_barcodes, reverse_barcodes)
        for _ in 1:2
            push!(reverse_barcodes, generate_distant_barcode(all_existing))
            push!(all_existing, reverse_barcodes[end])
        end
        cases = ["Donor1", "Donor2"]
        indices = DataFrame(forward_index=forward_barcodes, reverse_index=reverse_barcodes, case=cases)
        CSV.write(output_file, indices, delim='\t')
        return indices
    end

    function random_sequence(min_length::Int, max_length::Int)
        len = rand(min_length:max_length)
        return join(rand(NUCLEOTIDES, len))
    end

    function insert_mutation(seq::String, pos::Int, mutation_length::Int)
        insertion = random_sequence(mutation_length, mutation_length)
        return string(seq[1:pos], insertion, seq[pos+1:end])
    end

    function delete_mutation(seq::String, pos::Int, mutation_length::Int)
        return string(seq[1:pos-1], seq[pos+mutation_length:end])
    end

    function substitute_mutation(seq::String, pos::Int, mutation_length::Int)
        substitution = Char[]
        for i in 1:mutation_length
            original = seq[pos+i-1]
            new_base = rand(setdiff(NUCLEOTIDES, [original]))
            push!(substitution, new_base)
        end
        return string(seq[1:pos-1], join(substitution), seq[pos+mutation_length:end])
    end

    const MUTATION_APPLY = Dict{String, Function}(
        "insertion" => insert_mutation,
        "deletion" => delete_mutation,
        "substitution" => substitute_mutation,
    )

    function apply_random_mutation(seq::String, mutation_type::String, mutation_length::Int)
        safe_start = 50  # legacy generator keeps mutations away from both ends
        safe_end = min(length(seq) - 50, length(seq) - mutation_length)
        if safe_start >= safe_end
            safe_start = max(1, div(length(seq), 4))
            safe_end = min(length(seq) - mutation_length, 3 * div(length(seq), 4))
        end
        pos = rand(safe_start:safe_end)
        f = get(MUTATION_APPLY, mutation_type, nothing)
        f !== nothing && return f(seq, pos, mutation_length)
        error("Unknown mutation type: $mutation_type")
    end

    function append_unique_novel!(
        novel_records::Vector{FASTARecord},
        used_names::Set{String},
        used_seqs::Set{String},
        name::String,
        seq::String,
    )
        name in used_names && throw(ArgumentError("duplicate novel allele name: $name"))
        seq in used_seqs && throw(ArgumentError("duplicate novel sequence under name $name"))
        push!(novel_records, FASTARecord(name, seq))
        push!(used_names, name)
        push!(used_seqs, seq)
        return nothing
    end

    function unique_mutation(
        reference_seq::String,
        used_seqs::Set{String},
        mutation_types,
        mutation_lengths;
        max_tries::Int=200,
    )
        for _ in 1:max_tries
            variant = apply_random_mutation(reference_seq, rand(mutation_types), rand(mutation_lengths))
            variant in used_seqs || return variant
        end
        throw(ArgumentError("Could not generate a unique mutation sequence after $max_tries attempts"))
    end

    function unique_mutation(
        reference_seq::String,
        used_seqs::Set{String},
        mutation_type::String,
        mutation_length::Int;
        max_tries::Int=200,
    )
        for _ in 1:max_tries
            variant = apply_random_mutation(reference_seq, mutation_type, mutation_length)
            variant in used_seqs || return variant
        end
        throw(ArgumentError(
            "Could not generate a unique $mutation_type mutation of length $mutation_length " *
            "after $max_tries attempts"))
    end

    function generate_fasta_with_mutations(fasta_output::String, indices_output::String,
                                           reference_output::String, novel_output::String;
                                           n_reads::Int=100, base_length::Int=400)
        indices = generate_barcode_indices(indices_output)
        reference_seq = random_sequence(base_length, base_length)
        open(FASTA.Writer, reference_output) do writer
            write(writer, FASTARecord("REF*01", reference_seq))
        end

        records = FASTARecord[]
        novel_records = FASTARecord[]
        used_novel_names = Set{String}()
        used_novel_seqs = Set{String}()
        mutation_types = ["insertion", "deletion", "substitution"]
        mutation_lengths = [1, 3, 5]

        for (donor_idx, donor_row) in enumerate(eachrow(indices))
            forward = donor_row.forward_index
            reverse_bc = donor_row.reverse_index
            case_name = donor_row.case

            for mut_type in mutation_types
                for mut_len in mutation_lengths
                    mutated = unique_mutation(reference_seq, used_novel_seqs, mut_type, mut_len)
                    allele_name = "REF*$(donor_idx)_$(mut_type)_$(mut_len)"
                    append_unique_novel!(novel_records, used_novel_names, used_novel_seqs, allele_name, mutated)

                    for read_idx in 1:n_reads
                        prefix = random_sequence(20, 30)
                        suffix = random_sequence(20, 30)
                        full_read = forward * prefix * LEADER_SIGNAL * mutated * RSS_SIGNAL * suffix * reverse_bc
                        read_name = "$(case_name)_$(mut_type)_$(mut_len)_read$(read_idx)"
                        push!(records, FASTARecord(read_name, full_read))
                    end
                end
            end

            for read_idx in 1:n_reads
                prefix = random_sequence(20, 30)
                suffix = random_sequence(20, 30)
                full_read = forward * prefix * LEADER_SIGNAL * reference_seq * RSS_SIGNAL * suffix * reverse_bc
                read_name = "$(case_name)_reference_read$(read_idx)"
                push!(records, FASTARecord(read_name, full_read))
            end
        end

        # Add reference-derived novel records
        for (donor_idx, _) in enumerate(eachrow(indices))
            for i in 1:17
                variant = unique_mutation(reference_seq, used_novel_seqs, mutation_types, mutation_lengths)
                append_unique_novel!(
                    novel_records, used_novel_names, used_novel_seqs,
                    "REF*$(donor_idx)_extra_$(i)", variant)
            end
        end

        open(FASTA.Writer, fasta_output) do writer
            for record in records
                write(writer, record)
            end
        end
        open(FASTA.Writer, novel_output) do writer
            for record in novel_records
                write(writer, record)
            end
        end

        return records, indices
    end

    # ===================== Gene-type-aware synthetic reads (dispatch) =====================
    #
    # The legacy `generate_fasta_with_mutations` above only mutates the middle of a single
    # random reference. The functions below model real V/D/J architecture so the simulator
    # can produce the two cases the detection pipeline must handle: variants at the 3' end of
    # V genes, and very short D genes. Construction mirrors the layout the consumers expect:
    #   V read:  [prefix]                      gene  HEPTAMER spacer23 NONAMER  [suffix]
    #   D read:  [prefix] NONAMER spacer12 HEPTAMER gene  HEPTAMER spacer12 NONAMER [suffix]
    #   J read:  [prefix] NONAMER spacer23 HEPTAMER gene                          [suffix]

    """5' RSS flank emitted before the gene segment for a given gene type."""
    rss_prefix(::VGene) = ""
    rss_prefix(::DGene) = string(NONAMER, random_sequence(SPACER12, SPACER12), HEPTAMER)
    rss_prefix(::JGene) = string(NONAMER, random_sequence(SPACER23, SPACER23), HEPTAMER)

    """3' RSS flank emitted after the gene segment for a given gene type."""
    rss_suffix(::VGene) = string(HEPTAMER, random_sequence(SPACER23, SPACER23), NONAMER)
    rss_suffix(::DGene) = string(HEPTAMER, random_sequence(SPACER12, SPACER12), NONAMER)
    rss_suffix(::JGene) = ""

    """
        assemble_read(gene_type, gene_seq; barcode_f, barcode_r, flank)

    Build a synthetic read embedding `gene_seq` with the RSS architecture for `gene_type`,
    flanked by random padding (and optional barcodes). The gene appears as an exact substring
    so exact/heptamer search can recover it.
    """
    function assemble_read(gt::GeneType, gene_seq::AbstractString;
                           barcode_f::AbstractString="", barcode_r::AbstractString="", flank::Int=20)
        prefix = random_sequence(flank, flank)
        suffix = random_sequence(flank, flank)
        return string(barcode_f, prefix, rss_prefix(gt), gene_seq, rss_suffix(gt), suffix, barcode_r)
    end

    """
        simulate_d_read(gene_seq; flank) -> (read, flanks)

    Build a D-gene read with full pre/post RSS (nonamer–spacer12–heptamer on each side) and
    return both the read and the exact flank components, so the HSMM can be trained on the same
    motifs it is asked to detect.
    """
    function simulate_d_read(gene_seq::AbstractString; flank::Int=10)
        pre_spacer  = random_sequence(SPACER12, SPACER12)
        post_spacer = random_sequence(SPACER12, SPACER12)
        prefix = random_sequence(flank, flank)
        suffix = random_sequence(flank, flank)
        read = string(prefix, NONAMER, pre_spacer, HEPTAMER, gene_seq,
                      HEPTAMER, post_spacer, NONAMER, suffix)
        flanks = (pre_nonamer=NONAMER, pre_spacer=pre_spacer, pre_heptamer=HEPTAMER,
                  gene=String(gene_seq),
                  post_heptamer=HEPTAMER, post_spacer=post_spacer, post_nonamer=NONAMER)
        return read, flanks
    end

    """
        v_end_variant(seq; n_end)

    Return a V allele variant that differs from `seq` only within its last `n_end` nucleotides
    (each substituted to a different base). Mimics alleles separated only at the 3' V border,
    the hardest case for border/heptamer logic.
    """
    function v_end_variant(seq::AbstractString; n_end::Int=10)
        L = length(seq)
        n = min(n_end, L)
        head = seq[1:L-n]
        tail = String([rand(setdiff(NUCLEOTIDES, [seq[i]])) for i in (L-n+1):L])
        return string(head, tail)
    end

    """
        short_d_variant(; len)
        short_d_variant(germline; len)

    Generate a very short D segment (default 10 nt). With a `germline` argument, truncate it to
    `len` to mimic a shorter allelic variant; otherwise emit a fresh random short D.
    """
    short_d_variant(; len::Int=10) = random_sequence(len, len)
    function short_d_variant(germline::AbstractString; len::Int=10)
        L = length(germline)
        len >= L && return String(germline)
        return String(germline[1:len])
    end

    # ===================== Negative controls (false-positive checks) =====================

    """
        decoy_read(; len)

    A purely random read with no RSS architecture at all. Represents background sequence
    that is not a D gene; a correct detector must not call a confident D here.
    """
    decoy_read(; len::Int=100) = random_sequence(len, len)

    """
        invalid_d_read(gene_seq; flank)

    A read shaped like a D context (gene flanked by nonamer–spacer on each side) but with the
    conserved heptamers replaced by random 7-mers, i.e. a broken/invalid RSS. The recombination
    signal a real D depends on is absent, so it should score far below a genuine D read.
    """
    function invalid_d_read(gene_seq::AbstractString; flank::Int=10)
        scramble() = random_sequence(length(HEPTAMER), length(HEPTAMER))
        pre_spacer  = random_sequence(SPACER12, SPACER12)
        post_spacer = random_sequence(SPACER12, SPACER12)
        prefix = random_sequence(flank, flank)
        suffix = random_sequence(flank, flank)
        return string(prefix, NONAMER, pre_spacer, scramble(), gene_seq,
                      scramble(), post_spacer, NONAMER, suffix)
    end
end
