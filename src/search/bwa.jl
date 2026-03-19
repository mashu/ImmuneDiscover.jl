module Bwa
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

    function is_reverse_strand(aln::BurrowsWheelerAligner.LibBWA.mem_aln_t)
        return (aln.is_rev_is_alt_mapq_NM >> 0 & 0x01) == 0x00000001
    end

    # Valid DNA characters for BioSequences LongDNA{4}
    const VALID_DNA_CHARS = Set(['A', 'C', 'G', 'T', 'N', 'a', 'c', 'g', 't', 'n',
                                  'M', 'R', 'W', 'S', 'Y', 'K', 'V', 'H', 'D', 'B',
                                  'm', 'r', 'w', 's', 'y', 'k', 'v', 'h', 'd', 'b'])

    """
        reverse_complement_seq(seq::String) -> String

    Return the reverse complement of a DNA sequence.
    Validates input characters before conversion.
    """
    function reverse_complement_seq(seq::String)
        isempty(seq) && return seq
        # Validate all characters are valid DNA before conversion
        if all(c -> c in VALID_DNA_CHARS, seq)
            dna_seq = LongDNA{4}(seq)
            return String(reverse_complement(dna_seq))
        else
            # If sequence contains invalid characters, return as-is
            return seq
        end
    end

    function get_cigar_string(aln::BurrowsWheelerAligner.LibBWA.mem_aln_t)
        aln.n_cigar == 0 && return ""
        cigar_ops = ['M', 'I', 'D', 'N', 'S', 'H', 'P', 'X', '=']
        buf = IOBuffer()
        for i in 1:aln.n_cigar
            cigar_val = unsafe_load(aln.cigar, i)
            op_code = cigar_val & 0x0f
            op_len = cigar_val >> 4
            if op_code + 1 <= length(cigar_ops)
                print(buf, op_len, cigar_ops[op_code + 1])
            end
        end
        return String(take!(buf))
    end

    function calculate_leading_n_from_cigar(cigar::String)
        isempty(cigar) && return 0
        num_str = ""
        for ch in cigar
            if '0' <= ch <= '9'
                num_str *= ch
            else
                if !isempty(num_str)
                    return ch == 'N' ? parse(Int, num_str) : 0
                end
            end
        end
        return 0
    end

    function calculate_ref_span_from_cigar(cigar::String)
        isempty(cigar) && return 0
        ref_span = 0
        num_str = ""
        for ch in cigar
            if '0' <= ch <= '9'
                num_str *= ch
            else
                if !isempty(num_str)
                    len = parse(Int, num_str)
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

    function load_reference_sequences(fasta_path::String)
        ref_seqs = Dict{String, String}()
        @info "Loading reference sequences from $fasta_path"
        if endswith(fasta_path, ".gz")
            open(fasta_path) do file
                stream = GzipDecompressorStream(file)
                reader = FASTA.Reader(stream)
                for record in reader
                    ref_seqs[FASTA.identifier(record)] = FASTA.sequence(String, record)
                end
                close(reader)
                close(stream)
            end
        else
            open(FASTA.Reader, fasta_path) do reader
                for record in reader
                    ref_seqs[FASTA.identifier(record)] = FASTA.sequence(String, record)
                end
            end
        end
        @info "Loaded $(length(ref_seqs)) reference sequences"
        return ref_seqs
    end

    function extract_ref_subsequence(ref_seqs::Dict{String, String}, ref_name::String, pos::Int, query_length::Int)
        haskey(ref_seqs, ref_name) || return ""
        ref_seq = ref_seqs[ref_name]
        start_pos = pos + 1
        end_pos = min(start_pos + query_length - 1, length(ref_seq))
        (start_pos < 1 || start_pos > length(ref_seq)) && return ""
        return ref_seq[start_pos:end_pos]
    end

    function hamming_distance(seq1::String, seq2::String)
        length(seq1) != length(seq2) && return -1
        return sum(c1 != c2 for (c1, c2) in zip(seq1, seq2))
    end

    function bwa_sequences(genome_path, sequences, chromosome_name; tag="Primary Assembly")
        result = zeros(Bool, length(sequences))
        aligners = create_aligner(genome_path)
        discard = Accumulator{Tuple{String,String,String}, Int}()
        position = fill("", length(sequences))
        edit_distance = fill(-1, length(sequences))
        ref_sequence = fill("", length(sequences))
        orientation = fill("", length(sequences))
        cigar = fill("", length(sequences))

        ref_seqs_dict = Dict{String, Dict{String, String}}()
        for genome_file in genome_path
            ref_seqs_dict[genome_file] = load_reference_sequences(genome_file)
        end

        tag_regex = Regex(tag)

        @showprogress for (nallele, (name, sequence)) in enumerate(sequences)
            record = FASTA.Record(name, sequence)
            for (genome_file, aligner) in aligners
                alns = BurrowsWheelerAligner.align(aligner, record)
                if length(alns) == 0
                    @info "$name does **NOT** align to the $genome_file uniquely or at all (skipping)"
                    result[nallele] = false
                    continue
                end
                max_score = maximum(map(x -> x.score, alns))
                alns = filter(x -> x.score == max_score, alns)
                descriptions = map(x -> description(x, aligner), alns)
                index = map(x -> occursin(chromosome_name, x) && occursin(tag_regex, x), descriptions)
                match_descs = descriptions[index]
                nomatch = descriptions[.!index]

                if length(match_descs) > 0
                    result[nallele] = true
                    best_aln = first(alns[index])
                    ref_name = BurrowsWheelerAligner.refname(best_aln, aligner)
                    is_reverse = is_reverse_strand(best_aln)
                    position[nallele] = "$(ref_name):$(best_aln.pos)"
                    orientation[nallele] = is_reverse ? "-" : "+"
                    cigar_str = get_cigar_string(best_aln)
                    cigar[nallele] = cigar_str
                    ref_span = calculate_ref_span_from_cigar(cigar_str)
                    leading_n = calculate_leading_n_from_cigar(cigar_str)
                    extraction_length = ref_span > 0 ? ref_span : length(sequence)
                    adjusted_pos = is_reverse ? best_aln.pos - leading_n : best_aln.pos + leading_n
                    extracted_seq = extract_ref_subsequence(ref_seqs_dict[genome_file], ref_name, adjusted_pos, extraction_length)
                    ref_sequence[nallele] = is_reverse ? reverse_complement_seq(extracted_seq) : extracted_seq
                    if length(sequence) == length(ref_sequence[nallele])
                        edit_distance[nallele] = hamming_distance(sequence, ref_sequence[nallele])
                    else
                        edit_distance[nallele] = best_aln.is_rev_is_alt_mapq_NM >> 10 & 0x003fffff
                    end
                else
                    bad_chromosome = first(nomatch)
                    push!(discard, (genome_file, name, bad_chromosome))
                    result[nallele] = false
                end
            end
        end
        discarded = Vector{String}()
        for ((genome_file, name, chr), n) in discard
            @info "Discarded $name matching $chr in $genome_file (total $n)"
            push!(discarded, name)
        end
        CSV.write("/tmp/discarded.tsv", DataFrame(name=discarded), delim='\t')
        @info "Discarded sequence names written to /tmp/discarded.tsv"
        return result, position, edit_distance, ref_sequence, orientation, cigar
    end

    function handle_bwa(parsed_args, immunediscover_module, always_gz)
        @info "BWA search to filter candidates if they match correct chromosome"
        df = CSV.File(parsed_args["search"]["bwa"]["tsv"], delim='\t') |> DataFrame
        chromosome_name = parsed_args["search"]["bwa"]["chromosome"]
        outtsv = parsed_args["search"]["bwa"]["output"]
        genome = parsed_args["search"]["bwa"]["genome"]
        colname = parsed_args["search"]["bwa"]["colname"]
        colseq = parsed_args["search"]["bwa"]["colseq"]
        tag = parsed_args["search"]["bwa"]["tag"]

        sequences = map(eachrow(df)) do row
            concatenated_sequence = immunediscover_module.concatenate_columns(row, colseq)
            (row[colname], concatenated_sequence)
        end

        indices, position, edit_dist, ref_seq, orient, cigar_str = bwa_sequences(genome, sequences, chromosome_name, tag=tag)
        df[:, :position] = position
        df[:, :edit_distance] = edit_dist
        df[:, :ref_sequence] = ref_seq
        df[:, :orientation] = orient
        df[:, :cigar] = cigar_str
        @info "$(nrow(df)) sequences processed, $(sum(indices)) matched chromosome $chromosome_name"

        CSV.write(outtsv, df[indices, :], compress=true, delim='\t')
        @info "Filtered result saved in $outtsv"
    end

    export bwa_sequences, handle_bwa
end
