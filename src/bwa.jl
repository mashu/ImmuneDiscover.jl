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
    using CodecZlib

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

    """
        get_cigar_string(aln::mem_aln_t) -> String

    Extract CIGAR string from BWA alignment structure.
    BWA stores CIGAR as array of uint32 where lower 4 bits = operation, upper 28 bits = length.
    """
    function get_cigar_string(aln::BurrowsWheelerAligner.LibBWA.mem_aln_t)
        if aln.n_cigar == 0
            return ""
        end
        
        # CIGAR operations: 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7=X, 8==
        cigar_ops = ['M', 'I', 'D', 'N', 'S', 'H', 'P', 'X', '=']
        
        cigar_str = ""
        for i in 1:aln.n_cigar
            cigar_val = unsafe_load(aln.cigar, i)
            # Lower 4 bits = operation, upper 28 bits = length
            op_code = cigar_val & 0x0f
            op_len = cigar_val >> 4
            
            if op_code + 1 <= length(cigar_ops)
                cigar_str *= string(op_len) * cigar_ops[op_code + 1]
            end
        end
        
        return cigar_str
    end

    """
        calculate_leading_n_from_cigar(cigar::String) -> Int

    Calculate the number of leading 'N' operations in the CIGAR string.
    These represent skipped regions at the start that should be skipped in extraction.
    """
    function calculate_leading_n_from_cigar(cigar::String)
        if isempty(cigar)
            return 0
        end
        
        num_str = ""
        for ch in cigar
            if '0' <= ch <= '9'
                num_str *= ch
            else
                if !isempty(num_str)
                    len = parse(Int, num_str)
                    # If first operation is N, return its length
                    if ch == 'N'
                        return len
                    else
                        # First operation is not N, so no leading N
                        return 0
                    end
                end
            end
        end
        return 0
    end

    """
        calculate_ref_span_from_cigar(cigar::String) -> Int

    Calculate the total span on the reference sequence from a CIGAR string.
    Operations that consume reference: M, D, N, X, =
    Operations that don't consume reference: I, S, H, P
    Note: N (skipped region) IS INCLUDED to get the full genomic span
    """
    function calculate_ref_span_from_cigar(cigar::String)
        if isempty(cigar)
            return 0
        end
        
        ref_span = 0
        num_str = ""
        
        for ch in cigar
            if '0' <= ch <= '9'
                num_str *= ch
            else
                if !isempty(num_str)
                    len = parse(Int, num_str)
                    # Operations that consume reference bases (including N for skipped regions)
                    if ch in ('M', 'D', 'N', 'X', '=')
                        ref_span += len
                    end
                    num_str = ""
                end
            end
        end
        
        return ref_span
    end

    function create_aligner(genome_path)::Vector{Tuple{String, BurrowsWheelerAligner.Aligner}}
        return [(path, BurrowsWheelerAligner.Aligner(path)) for path in genome_path]
    end

    """
        load_reference_sequences(fasta_path::String) -> Dict{String, String}

    Load reference sequences from a FASTA file into a dictionary.
    Keys are sequence identifiers, values are the actual sequences.
    Supports both compressed (.gz) and uncompressed FASTA files.
    """
    function load_reference_sequences(fasta_path::String)
        ref_seqs = Dict{String, String}()
        @info "Loading reference sequences from $fasta_path"
        
        # Check if file is gzip-compressed by extension
        if endswith(fasta_path, ".gz")
            @info "Detected gzip-compressed file, decompressing on-the-fly"
            open(fasta_path) do file
                stream = GzipDecompressorStream(file)
                reader = FASTA.Reader(stream)
                try
                    for record in reader
                        seq_id = FASTA.identifier(record)
                        seq = FASTA.sequence(String, record)
                        ref_seqs[seq_id] = seq
                    end
                finally
                    close(reader)
                    close(stream)
                end
            end
        else
            # Uncompressed file
            open(FASTA.Reader, fasta_path) do reader
                for record in reader
                    seq_id = FASTA.identifier(record)
                    seq = FASTA.sequence(String, record)
                    ref_seqs[seq_id] = seq
                end
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
        cigar = fill("", length(sequences))
        
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
                    
                    # Get CIGAR string and calculate actual reference span
                    cigar_str = get_cigar_string(best_aln)
                    cigar[nallele] = cigar_str
                    ref_span = calculate_ref_span_from_cigar(cigar_str)
                    leading_n = calculate_leading_n_from_cigar(cigar_str)
                    
                    # Use ref_span if available, otherwise fall back to query length
                    extraction_length = ref_span > 0 ? ref_span : length(sequence)
                    
                    # Adjust start position for leading N operations
                    # For forward strand: leading N means skip those reference bases at the start
                    # For reverse strand: leading N in CIGAR (which is in query orientation after RC)
                    # represents query bases at the END of the original query (before RC)
                    # These map to reference bases BEFORE pos, so we need to extract from earlier
                    if is_reverse
                        # For reverse strand, go back by leading_n to include those bases
                        adjusted_pos = best_aln.pos - leading_n
                    else
                        # For forward strand, skip forward by leading_n
                        adjusted_pos = best_aln.pos + leading_n
                    end
                    
                    # Extract reference sequence at alignment position
                    extracted_seq = extract_ref_subsequence(
                        ref_seqs_dict[genome_file], 
                        ref_name, 
                        adjusted_pos, 
                        extraction_length
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
        return result, position, edit_distance, ref_sequence, orientation, cigar
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
        
        indices, position, edit_distance, ref_sequence, orientation, cigar_str = bwa_sequences(genome, sequences, chromosome_name, tag=tag)
        df[:, :position] = position
        df[:, :edit_distance] = edit_distance
        df[:, :ref_sequence] = ref_sequence
        df[:, :orientation] = orientation
        df[:, :cigar] = cigar_str
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