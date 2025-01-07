module blast
    using CSV
    using DataFrames
    using FASTX
    using DataStructures
    using BioAlignments
    include("data.jl")
    export blast_discover, save_to_fasta, accumulate_Ds, save_Ds

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
    """
    Performs a BLASTn search using the specified query and database, then saves the results in the specified output format.

    # Arguments
    - `query_file::String`: Path to the query file.
    - `database::String`: Path to the BLAST database.
    - `output_file::String`: Path where the output will be saved.
    - `output_format::String`: Output format (e.g., "6" for tabular).
    """
    function blastn(query_file::String, database::String, output_file::String;
        output_format::String = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand qseq", args="")
        nthreads = Threads.nthreads()
        makeblastdb_cmd = `makeblastdb -in $database -dbtype nucl -parse_seqids`
        blastn_cmd = `blastn $(split(args)) -num_threads $nthreads -query $query_file -db $database -out $output_file -outfmt $output_format`
        run(makeblastdb_cmd)
        run(blastn_cmd)
    end

    """
        replace_extension(filename, new_extension)

    Replace the extension of a file with a new one
    """
    function replace_extension(filename, new_extension; tag="")
        splitext = split(filename, '.')
        return first(splitext) * tag * '.' * new_extension
    end

    function accumulate_Ds(nameDs, dmux; ext_size=20)
        singleton = Vector{Tuple{String, String}}()
        for (name, sequence) in nameDs
            prefix = Accumulator{String,Int}()
            suffix = Accumulator{String,Int}()
            offset = ext_size
            for i in eachrow(dmux)
                pos = findfirst(sequence, i.genomic_sequence)
                genomic_sequence = i.genomic_sequence
                if pos === nothing
                    continue
                end
                start = minimum(pos)
                stop = maximum(pos)
                if start-offset < 1 || stop+offset > length(i.genomic_sequence)
                    continue
                end
                pre = genomic_sequence[start-offset:start-1]
                suf = genomic_sequence[stop+1:stop+offset]
                push!(prefix, pre)
                push!(suffix, suf)
            end
            if length(prefix) == 0 || length(suffix) == 0
                @warn "No common prefix or suffix for $name"
                continue
            end
            common_prefix, common_prefix_count = last(sort(collect(prefix),by=x->x[2]))
            common_suffix, common_suffix_count = last(sort(collect(suffix),by=x->x[2]))
            println("Extending $name with top $common_prefix ($common_prefix_count) prefixes and $common_suffix ($common_suffix_count) suffixes")
            push!(singleton, (name, "$(common_prefix)$sequence$(common_suffix)"))
        end
        return singleton
    end

    function save_Ds(extended_Ds, fasta_path)
        open(fasta_path, "w") do io
            for (name, sequence) in extended_Ds
                record = ">$name\n$sequence\n"
                write(io, record)
            end
        end
    end


    # Helper funciton to remove BLAST gaps before re-aligning
    nogaps(s) = replace(s,'-'=>"")

    function semilocal_alignment(query::String, reference::String)
        scoremodel = AffineGapScoreModel(EDNAFULL, gap_open=-20, gap_extend=-1)
        aln = alignment(pairalign(SemiGlobalAlignment(), nogaps(reference), nogaps(query), scoremodel))
        aligned_reference = last(split(string(aln.a)))
        aligned_query = last(split(string(aln.b)))
        m = match(r"[^-]+.*[^-]+", aligned_reference)
        pos = findfirst(m.match, aligned_reference)
        if maximum(pos) > length(aligned_query)
            @warn "Skipping outside of read result\n\t$aligned_query\n\t$aligned_reference"
            return "", -1
        end
        distance = sum([x != y for (x,y) in zip(aligned_query[pos], aligned_reference[pos])])
        return aligned_query[pos], distance
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

    """
        blast_discover(tsv_path, db_fasta; max_dist=10, min_count=5, min_frequency=0.1)

    Peform assignments and discovery of alleles based on BLAST results
    """
    function blast_discover(tsv_path, db_fasta; max_dist=10, min_count=5, min_frequency=0.1, min_length=290, pseudo="", args="", verbose=false, overwrite=false)
        min_ver = v"2.15.0"
        max_ver = v"2.16.0"
        is_valid, message = verify_blastn_version(min_ver, max_ver)
        @info message
        if !is_valid
            @warn "Please install BLAST version $min_ver - $max_ver"
            exit(1)
        end
        # Alter db_fasta
        db_p = Vector{Tuple{String, String}}()
        if pseudo != ""
            for (name, seq) in read_fasta(pseudo)
                push!(db_p, ("P"*name, seq))
            end
        end
        for (name, seq) in read_fasta(db_fasta)
            push!(db_p, (name, seq))
        end

        # Write the combined database to a FASTA file
        combined_db_fasta = replace_extension(db_fasta, "fasta", tag="-combined")
        save_to_fasta(db_p, combined_db_fasta)

        # Save query sequences to a FASTA file
        df = load_csv(tsv_path)
        @info "Read $(nrow(df)) rows from $tsv_path file before BLAST assignment"
        query_fasta = replace_extension(tsv_path, "fasta")
        blast_tsv = replace_extension(tsv_path, "blast")

        # Save the genomic sequences to a FASTA file
        query_sequences = collect.(eachrow(select_columns(df, [:well, :case, :name, :genomic_sequence])))
        save_to_fasta(query_sequences, query_fasta)

        # If file does not exist run blast
        if isfile(blast_tsv*".gz") && !overwrite
            @info "BLASTn results already exist $(blast_tsv).gz. Skipping BLASTn."
        else
            @info "Running BLASTn. This may take a while."
            # Run BLAST
            blastn(query_fasta, combined_db_fasta, blast_tsv, args=args)
            # Compress the BLAST results
            run(`gzip -f $blast_tsv`)
        end

        # Remove fasta reads file
        rm(query_fasta)
        rm(combined_db_fasta)

        # Read the BLAST results
        blast = CSV.File(blast_tsv*".gz", delim='\t') |> DataFrame
        @info "BLASTn results read from $(blast_tsv).gz: $(nrow(blast)) rows"
        # Rename the columns
        rename!(blast, [:qseqid, :sseqid, :pident, :length, :mismatch, :gapopen, :qstart, :qend, :sstart, :send, :evalue, :bitscore, :sstrand, :qseq])
        # Remove gaps from the query sequence
        blast.qseq = replace.(blast.qseq, "-" => "")

        # Discard pseudo genes
        filter!(x->!startswith(x.sseqid, "P"), blast)

        @info "BLASTn results after filtering pseudo genes: $(nrow(blast)) rows"
        if verbose
            @info "Saving BLAST results after filtering pseudo genes to $(blast_tsv)-pseudo.tsv"
            CSV.write(blast_tsv*"-pseudo.tsv", blast)
        end
        # Filter the results
        filter!(x->x.mismatch <= max_dist, blast)
        @info "BLASTn results after filtering mismatches: $(nrow(blast)) rows"
        if verbose
            @info "Saving BLAST results after filtering mismatches to $(blast_tsv)-mismatch.tsv"
            CSV.write(blast_tsv*"-mismatch.tsv", blast)
        end
        # Add the well and case columns
        read_name = Dict([(r.name,(r.well, r.case)) for r in eachrow(df)])
        blast[:,:well] = [read_name[x.qseqid][1] for x in eachrow(blast)]
        blast[:,:case] = [read_name[x.qseqid][2] for x in eachrow(blast)]

        # Count the number of matches
        clusters = combine(groupby(filter(x->x.mismatch .<= max_dist, blast), [:well, :case, :sseqid, :qseq, :mismatch]), :qseqid => length => :count)

        @info "Clusters after grouping by well, case, sseqid, qseq, mismatch: $(nrow(clusters)) rows"
        if verbose
            @info "Saving clusters to $(blast_tsv)-clusters.tsv"
            CSV.write(blast_tsv*"-clusters.tsv", clusters)
        end
        # Add the gene column
        transform!(clusters, :sseqid => ByRow(x -> split(x, "*")[1]) => :gene)

        # Compute the frequency of each gene in each well and case
        transform!(groupby(clusters, [:well, :case, :gene]), :count => (x -> x ./ maximum(x)) => :frequency)

        # Sort the results
        clusters_sorted = sort(clusters, [:well,:case, :sseqid], rev=false)

        # Filter the results
        filter!(x->x.count >= min_count, clusters_sorted)
        @info "Clusters after filtering by count: $(nrow(clusters_sorted)) rows"
        if verbose
            @info "Saving clusters after filtering by count to $(blast_tsv)-clusters-filtered-count.tsv"
            CSV.write(blast_tsv*"-clusters-filtered-count.tsv", clusters_sorted)
        end
        filter!(x->x.frequency >= min_frequency, clusters_sorted)
        @info "Clusters after filtering by frequency: $(nrow(clusters_sorted)) rows"
        if verbose
            @info "Saving clusters after filtering by frequency to $(blast_tsv)-clusters-filtered-freq.tsv"
            CSV.write(blast_tsv*"-clusters-filtered-freq.tsv", clusters_sorted)
        end
        filter!(x->length(x.qseq) >= min_length, clusters_sorted)
        @info "Clusters after filtering by minimum length: $(nrow(clusters_sorted)) rows"
        if verbose
            @info "Saving clusters after filtering by minimum length to $(blast_tsv)-clusters-filtered-minlen.tsv"
            CSV.write(blast_tsv*"-clusters-filtered-minlen.tsv", clusters_sorted)
        end
        # Save the results
        return clusters_sorted
    end
end
