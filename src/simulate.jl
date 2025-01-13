module simulate
    using FASTX
    using CSV
    using DataFrames
    using Random
    using MD5

    const NUCLEOTIDES = ['A', 'C', 'G', 'T']
    const RSS_SIGNAL = "CACAGTG"  # Example RSS signal
    const LEADER_SIGNAL = "GTTTTTGT"  # Example leader signal
    const BARCODE_LENGTH = 10  # Length of barcodes

    """
        hamming_distance(s1::String, s2::String)

    Calculate the Hamming distance between two strings
    """
    function hamming_distance(s1::String, s2::String)
        @assert length(s1) == length(s2)
        return sum(c1 != c2 for (c1, c2) in zip(s1, s2))
    end

    """
        generate_distant_barcode(existing_barcodes::Vector{String}, min_distance::Int=4)

    Generate a new barcode with minimum Hamming distance from existing ones
    """
    function generate_distant_barcode(existing_barcodes::Vector{String}, min_distance::Int=4)
        max_attempts = 1000
        for _ in 1:max_attempts
            candidate = join(rand(NUCLEOTIDES, BARCODE_LENGTH))
            # Check if candidate maintains minimum distance from all existing barcodes
            if isempty(existing_barcodes) ||
               all(hamming_distance(candidate, existing) >= min_distance
                   for existing in existing_barcodes)
                return candidate
            end
        end
        error("Could not generate barcode with required minimum distance after $max_attempts attempts")
    end

    """
        generate_barcode_indices(output_file::String="indices.tsv")

    Generate indices.tsv file with predefined barcodes for two donors
    """
    function generate_barcode_indices(output_file::String="indices.tsv")
        # Generate forward barcodes
        forward_barcodes = String[]
        for _ in 1:2
            push!(forward_barcodes, generate_distant_barcode(forward_barcodes))
        end

        # Generate reverse barcodes (different from forward ones)
        reverse_barcodes = String[]
        all_existing = vcat(forward_barcodes, reverse_barcodes)
        for _ in 1:2
            push!(reverse_barcodes, generate_distant_barcode(all_existing))
            push!(all_existing, reverse_barcodes[end])
        end

        cases = ["Donor1", "Donor2"]

        indices = DataFrame(
            forward_index = forward_barcodes,
            reverse_index = reverse_barcodes,
            case = cases
        )

        CSV.write(output_file, indices, delim='\t')
        return indices
    end

    """
        sequence_hash(seq; digits=4)

    Helper function for hashing the sequence
    Note this is not guaranteed to be unique but is simply visual indicator and this encoding is inherited from IgDiscover
    """
    function sequence_hash(seq; digits=4)
        "S"*lpad(string(parse(Int,bytes2hex(MD5.md5(seq))[(end-(digits-1)):end],base=16) % 10^digits), digits, '0')
    end

    """
        unique_name(name, sequence; digits=4)

    Wrapper function to update allele names
    """
    function unique_name(name, sequence; digits=4)
        "$(first(rsplit(name,"_")))_$(sequence_hash(sequence,digits=digits))"
    end

    """
        random_sequence(min_length::Int, max_length::Int)

    Generate a random DNA sequence of length between min_length and max_length
    """
    function random_sequence(min_length::Int, max_length::Int)
        length = rand(min_length:max_length)
        return join(rand(NUCLEOTIDES, length))
    end

    """
        insert_mutation(seq::String, pos::Int, mutation_length::Int)

    Insert random nucleotides at specified position
    """
    function insert_mutation(seq::String, pos::Int, mutation_length::Int)
        insertion = random_sequence(mutation_length, mutation_length)
        return string(seq[1:pos], insertion, seq[pos+1:end])
    end

    """
        delete_mutation(seq::String, pos::Int, mutation_length::Int)

    Delete nucleotides starting at specified position (inclusive)
    """
    function delete_mutation(seq::String, pos::Int, mutation_length::Int)
        return string(seq[1:pos-1], seq[pos+mutation_length:end])
    end

    """
        substitute_mutation(seq::String, pos::Int, mutation_length::Int)

    Replace nucleotides at specified position with random nucleotides
    """
    function substitute_mutation(seq::String, pos::Int, mutation_length::Int)
        # Generate random substitution different from original
        substitution = Char[]  # Using Char array for nucleotides
        for i in 1:mutation_length
            original = seq[pos+i-1]
            # Choose from nucleotides excluding the original
            new_base = rand(setdiff(NUCLEOTIDES, [original]))
            push!(substitution, new_base)
        end
        return string(seq[1:pos-1], join(substitution), seq[pos+mutation_length:end])
    end

    """
        apply_mutation(seq::String, mutation_type::Symbol, mutation_length::Int)

    Apply specified mutation type with specific length to sequence
    """
    function apply_mutation(seq::String, mutation_type::Symbol, mutation_length::Int)
        mid_point = length(seq) ÷ 2

        mutations = Dict(
            :mid_insert => () -> insert_mutation(seq, mid_point, mutation_length),
            :mid_delete => () -> delete_mutation(seq, mid_point, mutation_length),
            :mid_substitute => () -> substitute_mutation(seq, mid_point, mutation_length),
            :five_insert => () -> insert_mutation(seq, 1, mutation_length),
            # For 5' deletion, we want to delete the first n nucleotides
            :five_delete => () -> seq[mutation_length+1:end],
            :five_substitute => () -> substitute_mutation(seq, 1, mutation_length),
            :three_insert => () -> insert_mutation(seq, length(seq) - mutation_length + 1, mutation_length),
            :three_delete => () -> seq[1:end-mutation_length],
            :three_substitute => () -> substitute_mutation(seq, length(seq) - mutation_length + 1, mutation_length)
        )

        return mutations[mutation_type](), mutation_type, mutation_length
    end

    """
        generate_reference_sequence(base_seq::String)

    Generate reference sequence (original unmutated sequence without any affixes)
    """
    function generate_reference_sequence(base_seq::String, output_file::String="reference.fasta")
        record = FASTARecord("reference_unmutated", base_seq)
        writer = FASTA.Writer(open(output_file, "w"))
        write(writer, record)
        close(writer)
        return base_seq
    end

    """
        generate_fasta_with_mutations(input_fasta::String, indices_tsv::String, reference_fasta::String, novel_fasta::String)

    Generate input FASTA file with mutated sequences, reference FASTA, and novel FASTA.
    Also generates indices.tsv with barcodes for two donors.
    """
    function generate_fasta_with_mutations(input_fasta::String, indices_tsv::String, reference_fasta::String="reference.fasta", novel_fasta::String="novel.fasta")
        # First generate the barcode indices file
        indices = generate_barcode_indices(indices_tsv)

        # Generate original base sequence (290-310 nt)
        base_seq = random_sequence(290, 310)

        # Save original reference sequence before any modifications
        ref_seq = generate_reference_sequence(base_seq, reference_fasta)

        sequences = String[]
        labels = String[]

        # Store novel sequences (single copies before signals and replication)
        novel_records = Vector{FASTARecord}()

        mutation_types = [:mid_insert, :mid_delete, :mid_substitute,
                         :five_insert, :five_delete, :five_substitute,
                         :three_insert, :three_delete, :three_substitute]
        mutation_lengths = [1, 3, 6]  # Specific mutation lengths

        # Add signals to reference sequence to calculate padding
        base_with_signals = string(LEADER_SIGNAL, ref_seq, RSS_SIGNAL)
        target_length = 600

        # For each donor, calculate their specific padding based on their barcodes
        for (row_idx, row) in enumerate(eachrow(indices))
            prefix_barcode = row.forward_index
            suffix_barcode = row.reverse_index
            case_id = row.case

            # Calculate remaining padding needed for this donor's barcodes
            total_padding = target_length - (length(base_with_signals) +
                                          length(prefix_barcode) +
                                          length(suffix_barcode))
            prefix_padding = total_padding ÷ 2
            suffix_padding = total_padding - prefix_padding

            # Generate random padding specific to this donor
            donor_prefix_random = join(rand(NUCLEOTIDES, prefix_padding))
            donor_suffix_random = join(rand(NUCLEOTIDES, suffix_padding))

            # Generate mutations for this donor
            for mutation_type in mutation_types
                for mutation_length in mutation_lengths
                    # First generate single mutated sequence and save to novel.fasta
                    mutated_seq, mut_type, mut_len = apply_mutation(ref_seq, mutation_type, mutation_length)
                    novel_name = "novel_$(mut_type)_$(mut_len)nt"
                    novel_name_hash = unique_name(novel_name, mutated_seq)
                    push!(novel_records, FASTARecord(novel_name_hash, mutated_seq))

                    # Now replicate with donor's barcodes and padding
                    for _ in 1:34  # ~100 sequences per mutation type (34 * 3 ≈ 100)
                        final_seq = string(prefix_barcode, donor_prefix_random,
                                        LEADER_SIGNAL, mutated_seq, RSS_SIGNAL,
                                        donor_suffix_random, suffix_barcode)

                        push!(sequences, final_seq)
                        push!(labels, "$(case_id)_$(mut_type)_$(mut_len)nt")
                    end
                end
            end
        end

        # Save novel sequences
        writer = FASTA.Writer(open(novel_fasta, "w"))
        for record in novel_records
            write(writer, record)
        end
        close(writer)

        # Create and save FASTA records
        records = Vector{FASTARecord}()
        for (i, (seq, label)) in enumerate(zip(sequences, labels))
            record = FASTARecord("seq$(i)_$(label)", seq)
            push!(records, record)
        end

        writer = FASTA.Writer(open(input_fasta, "w"))
        for record in records
            write(writer, record)
        end
        close(writer)

        return records, indices
    end
end
