module blast
    using CSV
    using DataFrames
    using FASTX
    using DataStructures
    using BioAlignments
    using BioSequences
    using Folds
    include("data.jl")
    export blast_discover, save_to_fasta, accumulate_affixes, save_extended
    const columns = ["qseqid", "sseqid", "pident", "nident", "length", "mismatch", "gapopen", "qcovs", "qcovhsp", "qstart", "qend", "sstart", "send", "qlen", "slen", "evalue", "bitscore", "sstrand", "qseq"]

    """
        load_fasta(filepath::String) -> Vector{Tuple{String, String}}

    Load a FASTA file and return a vector of tuples, each containing the sequence name and sequence.

    # Arguments
    - `filepath::String`: Path to the FASTA file.

    # Returns
    - `Vector{Tuple{String, String}}`: A vector of tuples, where each tuple contains the sequence name and sequence.
    """
    function read_fasta(filepath::String)
        fasta_records = Vector{Tuple{String, String}}()
        open(filepath, "r") do stream
            for record in FASTA.Reader(stream)
                name = string(FASTA.identifier(record))
                sequence = string(FASTA.sequence(record))
                push!(fasta_records, (name, sequence))
            end
        end
        return fasta_records
    end

    load_csv(path::String; delim::Char='\t') = CSV.File(path, delim=delim) |> DataFrame
    select_columns(df::DataFrame, cols::Vector{Symbol}) = unique(select(df, cols))


    """
    Saves a collection of sequences to a FASTA file.

    # Arguments
    - `records::Vector{Tuple{String, String}}`: A vector of tuples, where each tuple contains a sequence name and the sequence itself.
    - `output_file::String`: The path to the output FASTA file.
    """
    function save_to_fasta(records, output_file::String)
        open(output_file, "w") do io
            for (well, case, name, sequence) in records
                record = ">$name $well $case\n$sequence\n"
                write(io, record)
            end
        end
    end

    function save_to_fasta(records::Vector{Tuple{String, String}}, output_file::String)
        open(output_file, "w") do io
            for (name, sequence) in records
                record = ">$name \n$sequence\n"
                write(io, record)
            end
        end
    end

    function time_command(cmd::Cmd)
        start_time = time()
        run(cmd)
        elapsed = time() - start_time
        return elapsed
    end

    """
    Performs a BLASTn search using the specified query and database, then saves the results in the specified output format.

    # Arguments
    - `query_file::String`: Path to the query file.
    - `database::String`: Path to the BLAST database.
    - `output_file::String`: Path where the output will be saved.
    - `output_format::String`: Output format (e.g., "6" for tabular).
    """
    function blastn(query_file::String, database::String, output_file::String;
        output_format::String = "6 $(join(columns, ' '))", args="")
        nthreads = Threads.nthreads()
        makeblastdb_cmd = `makeblastdb -in $database -dbtype nucl -parse_seqids`
        blastn_cmd = `blastn $(split(args)) -num_threads $nthreads -query $query_file -db $database -out $output_file -outfmt $output_format`
        run(makeblastdb_cmd)
        elapsed = time_command(blastn_cmd)
        @info "BLASTn search completed in $(elapsed) seconds."
    end

    """
        replace_extension(filename, new_extension)

    Replace the extension of a file with a new one
    """
    function replace_extension(filename, new_extension; tag="")
        splitext = split(filename, '.')
        return first(splitext) * tag * '.' * new_extension
    end

    function accumulate_affixes(alleles, dmux; forward_extension=20, reverse_extension=20)
        # Group alleles by gene
        gene_groups = Dict{String, Vector{Tuple{String,String}}}()
        for (name, sequence) in alleles
            gene = split(name, '*')[1]
            if !haskey(gene_groups, gene)
                gene_groups[gene] = Vector{Tuple{String,String}}()
            end
            push!(gene_groups[gene], (name, sequence))
        end

        # Process genes in parallel using Folds.map
        results = Folds.map(collect(gene_groups)) do (gene, gene_alleles)
            prefix = Accumulator{String,Int}()
            suffix = Accumulator{String,Int}()

            # Accumulate affixes for all alleles of this gene
            for (name, sequence) in gene_alleles
                for i in eachrow(dmux)
                    pos = findfirst(sequence, i.genomic_sequence)
                    if pos === nothing
                        continue
                    end

                    start = minimum(pos)
                    stop = maximum(pos)

                    # Handle forward extension
                    if forward_extension > 0
                        if start-forward_extension < 1
                            continue
                        end
                        pre = @view(i.genomic_sequence[start-forward_extension:start-1])
                        push!(prefix, String(pre))
                    else
                        push!(prefix, "")
                    end

                    # Handle reverse extension
                    if reverse_extension > 0
                        if stop+reverse_extension > length(i.genomic_sequence)
                            continue
                        end
                        suf = @view(i.genomic_sequence[stop+1:stop+reverse_extension])
                        push!(suffix, String(suf))
                    else
                        push!(suffix, "")
                    end
                end
            end

            # Find most common prefix and suffix for this gene
            if isempty(prefix) || isempty(suffix)
                @warn "Gene $gene has no flanking sequences found in any of its $(length(gene_alleles)) alleles!"
                return (gene, "", "", 0, 0)
            end

            common_prefix, common_prefix_count = last(sort(collect(prefix), by=x->x[2]))
            common_suffix, common_suffix_count = last(sort(collect(suffix), by=x->x[2]))

            if forward_extension > 0 || reverse_extension > 0
                extension_info = []
                if forward_extension > 0
                    push!(extension_info, "prefix $common_prefix ($common_prefix_count)")
                end
                if reverse_extension > 0
                    push!(extension_info, "suffix $common_suffix ($common_suffix_count)")
                end
                @info "Extending gene $gene with " * join(extension_info, " and ")
            end

            return (gene, common_prefix, common_suffix, common_prefix_count, common_suffix_count)
        end

        # Process results and create final output
        singleton = Vector{Tuple{String, String, String, String}}()

        # Apply extensions to ALL alleles using the gene's affixes
        for (gene, common_prefix, common_suffix, prefix_count, suffix_count) in results
            # Apply extensions to all alleles of this gene
            for (name, sequence) in gene_groups[gene]
                extended_sequence = string(common_prefix, sequence, common_suffix)
                push!(singleton, (name, extended_sequence, common_prefix, common_suffix))
            end
        end

        return singleton
    end

    function save_extended(extended_Ds, fasta_path)
        base_affixes = Vector{Tuple{String, String, String}}()
        open(fasta_path, "w") do io
            for (name, sequence, prefix, suffix) in extended_Ds
                record = ">$name\n$sequence\n"
                write(io, record)
                push!(base_affixes, (name, prefix, suffix))
            end
        end
        return base_affixes
    end

    # Helper funciton to remove BLAST gaps before re-aligning
    nogaps(s) = replace(s,'-'=>"")

    function check_affix_quality_warning(affix_length::Int, quality_threshold::Float64)
        """
        Check if quality threshold might be too strict for short affixes and issue warning.
        """
        if affix_length > 0 && affix_length <= 20 && quality_threshold > 0.5
            @warn "Quality threshold $(round(quality_threshold * 100, digits=1))% might be too strict for short affixes ($affix_length nt). Consider lowering --minquality for better sensitivity."
        end
    end

    function calculate_alignment_quality(pairs, affix_is_first::Bool=true)
        # Find the first and last non-gap positions in the affix
        if affix_is_first
            # Affix is in the first position of each pair
            affix_start = findfirst(p -> first(p) != DNA_Gap, pairs)
            affix_end = findlast(p -> first(p) != DNA_Gap, pairs)
        else
            # Affix is in the second position of each pair (reversed alignment)
            affix_start = findfirst(p -> last(p) != DNA_Gap, pairs)
            affix_end = findlast(p -> last(p) != DNA_Gap, pairs)
        end

        if affix_start === nothing || affix_end === nothing
            return 0.0
        end

        # Only consider the region where the affix aligns (no gaps in affix)
        aligned_region = pairs[affix_start:affix_end]

        # Count matches in the aligned region
        matches = sum(first.(aligned_region) .== last.(aligned_region))
        total = length(aligned_region)

        return matches / total
    end

    function safe_pairalign(seq1::LongDNA{4}, seq2::LongDNA{4}, scoremodel::AffineGapScoreModel)
        try
            return pairalign(SemiGlobalAlignment(), seq1, seq2, scoremodel)
        catch e
            if isa(e, BoundsError)
                @warn "Alignment failed due to bounds error with sequences: \nseq1: $seq1\nseq2: $seq2"
                return nothing
            end
            rethrow(e)
        end
    end

    function remove_gaps(query::LongDNA{4})::LongDNA{4}
        return LongDNA{4}(filter(nt -> nt != DNA_Gap, query))
    end

    function trim_sequence(query::LongDNA{4}, prefix::LongDNA{4}, suffix::LongDNA{4}, stats,
        scoremodel::AffineGapScoreModel=AffineGapScoreModel(EDNAFULL, gap_open=-10, gap_extend=-1);
        min_quality=0.75, sseqid="")

        # Input validation
        if length(query) == 0
            push!(stats.failed_genes[sseqid], "Empty query sequence")
            return nothing
        end

        partial_query = query
        stats.total_attempts += 1


        # The key insight: Extensions help find matches, but we trim back to the original core
        # Use semi-global alignment but with better quality calculation
        
        # Only trim prefix if it has length
        if length(prefix) > 0
            # Perform semi-global alignment: prefix against query
            prefix_aln_result = safe_pairalign(prefix, query, scoremodel)

            if prefix_aln_result === nothing
                push!(stats.failed_genes[sseqid], "Prefix alignment failed")
                stats.prefix_failures += 1
                return nothing
            end

            prefix_aln = alignment(prefix_aln_result)
            prefix_pairs = collect(prefix_aln)


            # Calculate quality: only count positions where prefix is not a gap
            prefix_positions = findall(p -> first(p) != DNA_Gap, prefix_pairs)
            if isempty(prefix_positions)
                push!(stats.failed_genes[sseqid], "No prefix content in alignment")
                stats.prefix_failures += 1
                return nothing
            end
            
            # Count matches only in prefix positions (no gaps in prefix)
            matches = sum(first(prefix_pairs[i]) == last(prefix_pairs[i]) for i in prefix_positions)
            prefix_quality = matches / length(prefix_positions)
            
            
            if prefix_quality < min_quality
                push!(stats.failed_genes[sseqid], "Poor prefix alignment quality ($(round(prefix_quality * 100, digits=1))% match)")
                stats.prefix_failures += 1
                return nothing
            end

            # Find where the prefix ends (last non-gap in prefix)
            aligned_prefix = first.(prefix_pairs)
            aligned_query = last.(prefix_pairs)
            
            prefix_end = findlast(x -> x != DNA_Gap, aligned_prefix)
            if prefix_end === nothing
                push!(stats.failed_genes[sseqid], "No prefix match found")
                stats.prefix_failures += 1
                return nothing
            end
            
            # Extract core sequence after prefix
            query_core_start = prefix_end + 1
            if query_core_start > length(aligned_query)
                push!(stats.failed_genes[sseqid], "Query too short after prefix")
                stats.prefix_failures += 1
                return nothing
            end
            
            core_part = aligned_query[query_core_start:end]
            partial_query = remove_gaps(LongDNA{4}(join(core_part)))
            if length(partial_query) == 0
                push!(stats.failed_genes[sseqid], "Empty core after prefix trimming")
                stats.prefix_failures += 1
                return nothing
            end
        end

        # Only trim suffix if it has length
        if length(suffix) > 0
            # Perform semi-global alignment: suffix against partial_query
            suffix_aln_result = safe_pairalign(suffix, partial_query, scoremodel)

            if suffix_aln_result === nothing
                push!(stats.failed_genes[sseqid], "Suffix alignment failed")
                stats.suffix_failures += 1
                return nothing
            end

            suffix_aln = alignment(suffix_aln_result)
            suffix_pairs = collect(suffix_aln)


            # Calculate quality: only count positions where suffix is not a gap
            suffix_positions = findall(p -> first(p) != DNA_Gap, suffix_pairs)
            if isempty(suffix_positions)
                push!(stats.failed_genes[sseqid], "No suffix content in alignment")
                stats.suffix_failures += 1
                return nothing
            end
            
            # Count matches only in suffix positions (no gaps in suffix)
            matches = sum(first(suffix_pairs[i]) == last(suffix_pairs[i]) for i in suffix_positions)
            suffix_quality = matches / length(suffix_positions)
            
            
            if suffix_quality < min_quality
                push!(stats.failed_genes[sseqid], "Poor suffix alignment quality ($(round(suffix_quality * 100, digits=1))% match)")
                stats.suffix_failures += 1
                return nothing
            end

            # Find where the suffix starts (first non-gap in suffix)
            aligned_suffix = first.(suffix_pairs)
            aligned_query = last.(suffix_pairs)
            
            suffix_start = findfirst(x -> x != DNA_Gap, aligned_suffix)
            if suffix_start === nothing
                push!(stats.failed_genes[sseqid], "No suffix match found")
                stats.suffix_failures += 1
                return nothing
            end
            
            if suffix_start <= 1
                push!(stats.failed_genes[sseqid], "Suffix starts too early")
                stats.suffix_failures += 1
                return nothing
            end
            
            # Extract core sequence before suffix
            core_part = aligned_query[1:suffix_start-1]
            partial_query = remove_gaps(LongDNA{4}(join(core_part)))
            if length(partial_query) == 0
                push!(stats.failed_genes[sseqid], "Empty core after suffix trimming")
                stats.suffix_failures += 1
                return nothing
            end
        end

        # Check if remaining core sequence is long enough for specificity
        core_length = length(partial_query)
        if core_length < 10
            push!(stats.failed_genes[sseqid], "Core sequence too short ($core_length < 10 nt)")
            stats.prefix_failures += 1
            return nothing
        end


        push!(stats.successful_trims[sseqid], "Trimmed successfully")
        return partial_query
    end

    function compute_edit_distance(query::String, reference::String)
       aln = pairalign(LevenshteinDistance(), reference, query)
       return score(aln)
    end

    function trim_and_align_sequence(query::String, prefix::String, suffix::String, reference::String, stats; min_quality=0.75, sseqid="")
        nogaps(s) = replace(s, '-' => "")

        # Convert strings to DNA sequences
        query_dna = BioSequences.LongDNA{4}(query)
        prefix_dna = BioSequences.LongDNA{4}(prefix)
        suffix_dna = BioSequences.LongDNA{4}(suffix)

        trimmed_dna = trim_sequence(query_dna, prefix_dna, suffix_dna, stats, min_quality=min_quality, sseqid=sseqid)
        # If trimming failed, return empty sequence and negative distance
        if trimmed_dna === nothing
            return "", -1
        end
        trimmed = nogaps(String(trimmed_dna))
        distance = compute_edit_distance(trimmed, reference)

        return trimmed, distance
    end

    function verify_blastn_version(min_version::VersionNumber, max_version::VersionNumber)
        try
            # Run blastn to get version info
            cmd_output = read(`blastn -version`, String)

            # Extract version number from output
            # Typical output looks like: "blastn: 2.13.0+"
            version_match = match(r"blastn:\s+(\d+\.\d+\.\d+)", cmd_output)

            if isnothing(version_match)
                return (false, "Could not parse blastn version")
            end

            current_version = VersionNumber(version_match[1])

            if min_version ≤ current_version ≤ max_version
                return (true, "blastn version $current_version is within acceptable range")
            else
                return (false, "blastn version $current_version is outside acceptable range ($min_version - $max_version)")
            end

        catch e
            if isa(e, SystemError)
                return (false, "blastn command not found in PATH")
            else
                return (false, "Error checking blastn version: $(string(e))")
            end
        end
    end

    function edge(qseq, read)
        m = findfirst(nogaps(qseq), read)
        if isnothing(m)
            return missing, missing
        end
        five_prime, three_prime = extrema(m)
        return five_prime, length(read)-three_prime
    end

    """
    Peform assignments and discovery of alleles based on BLAST results
    """
    function blast_discover(tsv_path, combined_db_fasta; max_dist=10, min_edge=10, min_scov=0.1, args="", verbose=false, overwrite=false)
        min_ver = v"2.15.0"
        max_ver = v"2.17.0"
        is_valid, message = verify_blastn_version(min_ver, max_ver)
        @info message
        if !is_valid
            @warn "Please install BLAST version $min_ver - $max_ver"
            exit(1)
        end

        # Save query sequences to a FASTA file
        df = load_csv(tsv_path)
        @info "Read $(nrow(df)) rows from $tsv_path file before BLAST assignment"
        query_fasta = replace_extension(tsv_path, "fasta")
        blast_tsv = replace_extension(tsv_path, "blast")
        unique_df = unique(df, :name)
        @info "Unique reads: $(nrow(unique_df))"
        if nrow(df) != nrow(unique_df)
            @error "Duplicated names found in the input demultiplex file. Using unique ones but this is likely user error!"
        end
        # If there are duplicated names, use the unique ones
        df = unique_df

        # Save the genomic sequences to a FASTA file
        query_sequences = collect.(eachrow(select_columns(df, [:well, :case, :name, :genomic_sequence])))
        save_to_fasta(query_sequences, query_fasta)

        # If file does not exist run blast
        blast_file = blast_tsv * ".gz"
        file_exists = isfile(blast_file)
        @info "BLAST cache check: file='$blast_file', exists=$file_exists, overwrite=$overwrite"
        if file_exists && !overwrite
            @info "BLASTn results already exist $blast_file. Skipping BLASTn."
        else
            @info "Running BLASTn. This may take a while."
            # Run BLAST
            blastn(query_fasta, combined_db_fasta, blast_tsv, args=args)
            # Compress the BLAST results
            run(`gzip -f $blast_tsv`)
        end

        # Remove fasta reads file
        rm(query_fasta)

        # Read the BLAST results
        blast = CSV.File(blast_tsv*".gz", delim='\t', header=columns) |> DataFrame
        @info "BLASTn results read from $(blast_tsv).gz: $(nrow(blast)) rows"
        if nrow(blast) == 0
            @error "No BLASTn results found (wrong BLAST parameters?). Exiting."
            exit(1)
        end

        # Take best hits from BLAST based on coverage, identity, and bitscore
        blast = combine(groupby(blast, :qseqid), x -> first(sort(x, [:pident, :qcovhsp, :qcovs, :bitscore], rev=true)))

        # Add edge information as well as original genomic sequence
        leftjoin!(blast, df, on=:qseqid => :name)
        transform!(blast, [:qseq, :genomic_sequence] => ByRow(edge) => [:five_prime_edge, :three_prime_edge])

        @info "BLASTn results after taking best hits: $(nrow(blast)) rows"
        # Filter out targets too close to edge
        filter!(x->x.five_prime_edge > min_edge && x.three_prime_edge > min_edge, blast)
        @info "BLASTn results after filtering edge > $min_edge: $(nrow(blast)) rows"

        # Add subject coverage to detect truncations
        transform!(blast, [:length, :slen] => ByRow((len, slen) -> len/slen) => :scov)
        filter!(x->x.scov > min_scov, blast)
        @info "BLASTn results after filtering subject (database) coverage > $min_scov: $(nrow(blast)) rows"

        # Remove gaps from the query sequence
        transform!(blast, :qseq => ByRow(x -> replace(x, "-" => "")) => :qseq)

        # Discard pseudo genes
        filter!(x->!startswith(x.sseqid, "P"), blast)
        @info "BLASTn results after filtering pseudo genes: $(nrow(blast)) rows"
        if verbose
            @info "Saving BLAST results after filtering pseudo genes to $(blast_tsv)-pseudo.tsv"
            CSV.write(blast_tsv*"-pseudo.tsv", blast)
        end

        # Add the well and case columns
        read_name = Dict([(r.name,(r.well, r.case)) for r in eachrow(df)])
        blast[:,:well] = [read_name[x.qseqid][1] for x in eachrow(blast)]
        blast[:,:case] = [read_name[x.qseqid][2] for x in eachrow(blast)]

        # Count the number of matches
        clusters = combine(groupby(blast, [:well, :case, :sseqid, :qseq, :mismatch]), :qseqid => length => :full_count)

        @info "Clusters after grouping by well, case, sseqid, qseq, mismatch: $(nrow(clusters)) rows"
        if verbose
            @info "Saving clusters to $(blast_tsv)-clusters.tsv"
            CSV.write(blast_tsv*"-clusters.tsv", clusters)
        end

        # Mandatory filter to reduce the number of clusters by maximum distance
        filter!(x->x.mismatch <= max_dist, clusters)
        @info "BLASTn results after filtering mismatches <= $max_dist: $(nrow(clusters)) rows"
        if verbose
            @info "Saving BLAST results after filtering mismatches to $(blast_tsv)-clusters-mismatch.tsv"
            CSV.write(blast_tsv*"-clusters-mismatch.tsv", clusters)
        end

        # Add the gene column
        transform!(clusters, :sseqid => ByRow(x -> split(x, "*")[1]) => :gene)

        # Compute the frequency of each gene in each well and case
        transform!(groupby(clusters, [:well, :case, :gene]), :full_count => (x -> x ./ maximum(x)) => :full_frequency)

        # Sort the results
        clusters_sorted = sort(clusters, [:well,:case, :sseqid], rev=false)

        # Save the results
        return clusters_sorted
    end
end
