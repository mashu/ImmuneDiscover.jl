module bwa
    using BurrowsWheelerAligner
    using BurrowsWheelerAligner.LibBWA: mem_aln_t
    using BurrowsWheelerAligner: Aligner
    using FASTX
    import FASTX.description
    using ProgressMeter
    using DataStructures
    using DataFrames
    using CSV
    using BioSequences

    function description(aln::BurrowsWheelerAligner.LibBWA.mem_aln_t, aligner::BurrowsWheelerAligner.Aligner)
        anns = BurrowsWheelerAligner.LibBWA.unsafe_load(aligner.index.bns).anns
        anno = BurrowsWheelerAligner.LibBWA.unsafe_load(anns, aln.rid+1).anno
        BurrowsWheelerAligner.LibBWA.unsafe_string(anno)
    end

    """
        is_reverse_strand(aln::mem_aln_t) -> Bool

    Check if the alignment is on the reverse strand.
    Bit 0 of is_rev_is_alt_mapq_NM indicates reverse strand.
    """
    function is_reverse_strand(aln::BurrowsWheelerAligner.LibBWA.mem_aln_t)
        return (aln.is_rev_is_alt_mapq_NM >> 0 & 0x01) == 0x00000001
    end

    """
        reverse_complement_seq(seq::String) -> String

    Return the reverse complement of a DNA sequence.
    """
    function reverse_complement_seq(seq::String)
        if isempty(seq)
            return seq
        end
        try
            dna_seq = LongDNA{4}(seq)
            return String(reverse_complement(dna_seq))
        catch
            # If sequence contains invalid characters, return as-is
            return seq
        end
    end

    function create_aligner(genome_path)::Vector{Tuple{String, BurrowsWheelerAligner.Aligner}}
        return [(path, BurrowsWheelerAligner.Aligner(path)) for path in genome_path]
    end

    """
        load_reference_sequences(fasta_path::String) -> Dict{String, String}

    Load reference sequences from a FASTA file into a dictionary.
    Keys are sequence identifiers, values are the actual sequences.
    """
    function load_reference_sequences(fasta_path::String)
        ref_seqs = Dict{String, String}()
        @info "Loading reference sequences from $fasta_path"
        open(FASTA.Reader, fasta_path) do reader
            for record in reader
                seq_id = FASTA.identifier(record)
                seq = FASTA.sequence(String, record)
                ref_seqs[seq_id] = seq
            end
        end
        @info "Loaded $(length(ref_seqs)) reference sequences"
        return ref_seqs
    end

    """
        extract_ref_subsequence(ref_seqs::Dict, ref_name::String, pos::Int, query_length::Int) -> String

    Extract a subsequence from the reference at the given position.
    pos is 0-based (BWA convention), query_length is the query sequence length.
    """
    function extract_ref_subsequence(ref_seqs::Dict{String, String}, ref_name::String, pos::Int, query_length::Int)
        if !haskey(ref_seqs, ref_name)
            return ""
        end
        ref_seq = ref_seqs[ref_name]
        # pos is 0-based in BWA, convert to 1-based for Julia
        start_pos = pos + 1
        end_pos = min(start_pos + query_length - 1, length(ref_seq))
        
        if start_pos < 1 || start_pos > length(ref_seq)
            return ""
        end
        
        return ref_seq[start_pos:end_pos]
    end

    function bwa_sequences(genome_path, sequences, chromosome_name; tag="Primary Assembly")
        result = zeros(Bool, length(sequences))
        aligners = create_aligner(genome_path)
        discard = Accumulator{Tuple{String,String,String}, Int}()
        position = fill("", length(sequences))
        edit_distance = fill(-1, length(sequences))  # -1 indicates no alignment
        ref_sequence = fill("", length(sequences))
        orientation = fill("", length(sequences))
        
        # Load reference sequences from FASTA files
        ref_seqs_dict = Dict{String, Dict{String, String}}()
        for genome_file in genome_path
            ref_seqs_dict[genome_file] = load_reference_sequences(genome_file)
        end
        
        # Compile tag regex once before the loop
        tag_regex = Regex(tag)
        
        @showprogress for (nallele, (name, sequence)) in enumerate(sequences)
            record = FASTA.Record(name, sequence)
            for (genome_file, aligner) in aligners
                alns = BurrowsWheelerAligner.align(aligner, record)
                # Pick alignments with the highest score
                if length(alns) == 0
                    @info "$name does **NOT** align to the $genome_file uniquely or at all (skipping)"
                    result[nallele] = false
                    continue
                end
                max_score = maximum(map(x->(x.score), alns))
                alns = filter(x->x.score == max_score, alns)
                descriptions = map(x->description(x, aligner), alns)

                index = map(x->occursin(chromosome_name, x) && occursin(tag_regex, x), descriptions)
                match = descriptions[index]
                nomatch = descriptions[.!index]

                if length(match) > 0
                    result[nallele] = true
                    best_aln = first(alns[index])
                    ref_name = BurrowsWheelerAligner.refname(best_aln, aligner)
                    is_reverse = is_reverse_strand(best_aln)
                    
                    position[nallele] = "$(ref_name):$(best_aln.pos)"
                    # Extract NM (edit distance) from packed field: bits 10-31 contain NM
                    # Format: is_rev:1, is_alt:1, mapq:8, NM:22
                    edit_distance[nallele] = best_aln.is_rev_is_alt_mapq_NM >> 10 & 0x003fffff
                    
                    # Set orientation: + for forward, - for reverse
                    orientation[nallele] = is_reverse ? "-" : "+"
                    
                    # Extract reference sequence at alignment position
                    extracted_seq = extract_ref_subsequence(
                        ref_seqs_dict[genome_file], 
                        ref_name, 
                        best_aln.pos, 
                        length(sequence)
                    )
                    
                    # If alignment is on reverse strand, reverse complement the reference
                    # to match the query orientation
                    if is_reverse
                        ref_sequence[nallele] = reverse_complement_seq(extracted_seq)
                    else
                        ref_sequence[nallele] = extracted_seq
                    end
                else
                    bad_chromosome = first(nomatch)
                    push!(discard, (genome_file, name, bad_chromosome))
                    result[nallele] = false 
                end
            end
        end
        discarded = Vector{String}()
        for ((genome_file, name, chr),n) in discard
            @info "Discarded $name matching $chr in $genome_file (total $n)"
            push!(discarded, name)
        end
        CSV.write("/tmp/discarded.tsv", DataFrame(name=discarded), delim='\t')
        @info "Discarded sequence names written to /tmp/discarded.tsv"
        return result, position, edit_distance, ref_sequence, orientation
    end

    """
        handle_bwa(parsed_args, immunediscover_module, always_gz)

    Handle BWA command from CLI arguments
    """
    function handle_bwa(parsed_args, immunediscover_module, always_gz)
        @info "BWA search to filter candidates if they match correct chromosome"
        df = CSV.File(parsed_args["analyze"]["bwa"]["tsv"], delim='\t') |> DataFrame
        chromosome_name = parsed_args["analyze"]["bwa"]["chromosome"]
        outtsv = parsed_args["analyze"]["bwa"]["output"]
        genome = parsed_args["analyze"]["bwa"]["genome"]
        colname = parsed_args["analyze"]["bwa"]["colname"]
        colseq = parsed_args["analyze"]["bwa"]["colseq"]
        tag = parsed_args["analyze"]["bwa"]["tag"]
        @info "Using columns: $colname, $colseq"
        @info "Using genome: $genome"
        @info "Filter by chromosome: $chromosome_name and $tag tag"
        
        sequences = map(eachrow(df)) do row
            concatenated_sequence = immunediscover_module.concatenate_columns(row, colseq)
            (row[colname], concatenated_sequence)
        end
        
        indices, position, edit_distance, ref_sequence, orientation = bwa_sequences(genome, sequences, chromosome_name, tag=tag)
        df[:, :position] = position
        df[:, :edit_distance] = edit_distance
        df[:, :ref_sequence] = ref_sequence
        df[:, :orientation] = orientation
        @info "$(length(sequences))"
        @info "$(nrow(df)) sequences matched chromosome $chromosome_name"
        @info "$(length(indices)) sequences matched chromosome $chromosome_name"
        @info "$(indices[1:10])"
        @info "$(df[1:10,:])"
        
        CSV.write(outtsv, df[indices,:], compress=true, delim='\t')
        @info "Filtered result saved in $outtsv"
    end

    export bwa_sequences, handle_bwa
end