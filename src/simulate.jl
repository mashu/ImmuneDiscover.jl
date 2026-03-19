module Simulate
    using FASTX
    using CSV
    using DataFrames
    using Random
    using MD5

    using ..Data: unique_name, sequence_hash

    const NUCLEOTIDES = ['A', 'C', 'G', 'T']
    const RSS_SIGNAL = "CACAGTG"
    const LEADER_SIGNAL = "GTTTTTGT"
    const BARCODE_LENGTH = 10

    export generate_fasta_with_mutations, unique_name, sequence_hash

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

    function apply_random_mutation(seq::String, mutation_type::String, mutation_length::Int)
        safe_start = max(1, 50)
        safe_end = min(length(seq) - 50, length(seq) - mutation_length)
        if safe_start >= safe_end
            safe_start = max(1, div(length(seq), 4))
            safe_end = min(length(seq) - mutation_length, 3 * div(length(seq), 4))
        end
        pos = rand(safe_start:safe_end)
        if mutation_type == "insertion"
            return insert_mutation(seq, pos, mutation_length)
        elseif mutation_type == "deletion"
            return delete_mutation(seq, pos, mutation_length)
        elseif mutation_type == "substitution"
            return substitute_mutation(seq, pos, mutation_length)
        else
            error("Unknown mutation type: $mutation_type")
        end
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
        mutation_types = ["insertion", "deletion", "substitution"]
        mutation_lengths = [1, 3, 5]

        for (donor_idx, donor_row) in enumerate(eachrow(indices))
            forward = donor_row.forward_index
            reverse_bc = donor_row.reverse_index
            case_name = donor_row.case

            for mut_type in mutation_types
                for mut_len in mutation_lengths
                    mutated = apply_random_mutation(reference_seq, mut_type, mut_len)
                    allele_name = "REF*$(donor_idx)_$(mut_type)_$(mut_len)"
                    push!(novel_records, FASTARecord(allele_name, mutated))

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
                variant = apply_random_mutation(reference_seq, rand(mutation_types), rand(mutation_lengths))
                push!(novel_records, FASTARecord("REF*$(donor_idx)_extra_$(i)", variant))
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
end
