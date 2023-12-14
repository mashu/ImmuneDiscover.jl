module data
    using FASTX
    using CodecZlib
    using ProgressMeter
    using CodecZlib
    using DataFrames
    using Statistics
    using UnicodePlots
    using CSV
    using MD5

    export load_fasta, plotgenes, unique_name, sequence_hash, load_demultiplex

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

    Wrapper funciton to update allele names
    """
    function unique_name(name, sequence; digits=4)
        "$(first(rsplit(name,"_")))_$(sequence_hash(sequence,digits=digits))"
    end

    """
        write_fastq(path, records)

    Save records with sequences and quality scores to fastq file
    """
    function write_fastq(path, records)
        FASTQ.Writer(open(path, "w")) do writer
            for record in records
                write(writer, record)
            end
        end
    end

    """
        write_gz_fastq(path, records)

    Save records with sequences and quality scores to fastq file
    """
    function write_gz_fastq(path, records)
        FASTQ.Writer(GzipCompressorStream(open(path, "w"))) do writer
            for record in records
                write(writer, record)
            end
        end
    end

    """
        write_fasta(path, records)

    Save records with sequences to fasta file
    """
    function write_fasta(path, records)
        FASTA.Writer(open(path, "w")) do writer
            for record in records
                write(writer, record)
            end
        end
    end

    """
        write_gz_fasta(path, records)

    Save records with sequences to fasta file
    """
    function write_gz_fasta(path, records)
        FASTA.Writer(GzipCompressorStream(open(path, "w"))) do writer
            for record in records
                write(writer, record)
            end
        end
    end

    """
        load_demultiplex(path)

    Load demultiplex file
    """
    function load_demultiplex(path; limit=nothing)
        table = CSV.File(path, delim='\t', types=Dict(:case => String), limit=limit) |> DataFrame
        @assert all([name in names(table) for name in ["well","case","name","genomic_sequence"]]) "File must contain following columns: well, case, name, genomic_sequence"
        return table
    end

    """
        validate_identifier(id::String)

    Validate that the identifier follows the specified format:
    - Contains a star (*) that separates it into two parts.
    - Each part contains at least one capital letter and one number.
    - If the second part does not contain 'S', it must be only numbers.
    """
    function validate_identifier(id)
        # Regular expression pattern
        # Part 1: [A-Z].*\d.* - At least one capital letter and one number before '*'
        # Part 2: (\*.*S.*\d.*|.*\*.*\d+$) - After '*', either contains 'S' and numbers or is only numbers
        pattern = r"[A-Z].*\d.*\*(.*S.*\d.*|.*\d+$)"
        letters = collect(id)
        underscore = sum(map(x->'_' == x, letters))
        star = sum(map(x->'*' == x, letters))
        # Check if the identifier matches the pattern
        return occursin(pattern, id) && (underscore <= 1) && (star <= 1)
    end

    """
        validate_sequence(seq::String)

    Check if the DNA sequence contains only 'A', 'T', 'G', and 'C'.
    """
    function validate_sequence(seq::String)
        # Define the allowed characters
        allowed_chars = Set(['A', 'T', 'G', 'C'])

        # Check each character in the sequence
        for char in seq
            if !(char in allowed_chars)
                return false
            end
        end
        return true
    end

    """
        load_fasta(path; validate=true)

    Function to load gz compressed or uncompressed FASTA file. If validate=true, 
    throws an error upon encountering a duplicated sequence or identifier.
    """
    function load_fasta(path; validate=false)
        records = Vector{Tuple{String, String}}()
        seen_descriptions = Set{String}()
        seen_sequences = Set{String}()

        open(FASTA.Reader, path) do reader
            for record in reader
                desc = string(FASTA.description(record))
                seq = string(FASTA.sequence(record))

                # By default everything is fine
                error = false

                if validate
                    # Check if sequence is valid
                    if length(seq) == 0
                        throw(ErrorException("Empty sequence found: $(desc)"))
                        error = true
                    end

                    # Check if identifier is valid
                    if !validate_identifier(desc)
                        throw(ErrorException("Invalid identifier found: $(desc)"))
                        error = true
                    end

                    if length(desc) == 0
                        throw(ErrorException("Empty identifier found: $(desc)"))
                        error = true
                    end

                    # Check for duplicate description
                    if desc in seen_descriptions
                        throw(ErrorException("Duplicate identifier found: $(desc)"))
                        error = true
                    end

                    # Check for duplicate sequence
                    if seq in seen_sequences
                        throw(ErrorException("Duplicate sequence found: $(seq)"))
                        error = true
                    end
                end

                # Exit on error
                if error
                    exit(1)
                end

                push!(seen_descriptions, desc)
                push!(seen_sequences, seq)
                push!(records, (desc, seq))
            end
        end
        return records
    end

    """
        load_fastq(path)

    Function to process either gz compressed or uncompressed FASTQ file with an `f` callback
    """
    function process_fastq(f, path)
        if endswith(path, ".gz")
            @info "Compressed file detected"
            reader = FASTQ.Reader(GzipDecompressorStream(open(path)))
        else
            @info "Uncompressed file detected"
            reader = FASTQ.Reader(open(path))
        end
        prog = ProgressUnknown(desc="Processing FASTQ reads")
        for record in reader
            f(record)
            ProgressMeter.next!(prog)
        end
        ProgressMeter.finish!(prog)
        close(reader)
    end

    """
        sort_by_median(d::Dict{String, Vector{Float64}})

    Sort the keys of a dictionary by the median of their associated Vector{Float64}.

    # Arguments
    - `d::Dict{String, Vector{Float64}}`: The dictionary to sort.

    # Returns
    - `Array{Tuple{String, Vector{Float64}}, 1}`: An array of (key, value) tuples, sorted by median of the value.
    """
    function sort_by_median(d::Dict{String, Vector{Float64}})
        # Calculate median for each key-value pair and store it as a tuple (key, value, median)
        median_tuples = [(k, v, median(v)) for (k, v) in d]

        # Sort the array of tuples based on the median value
        sort!(median_tuples, by = x -> x[3], rev=true)

        # Create an array of (key, value) tuples, ordered by their median
        sorted_by_median = [(t[1], t[2]) for t in median_tuples]

        return sorted_by_median
    end

    """
    plotgenes(df::DataFrame, output::String)

    Function to plotgenes results of hamming search
    """
    function plotgenes(df::DataFrame)
        # Compute log2 counts
        df[!, :log2_count] = log2.(df[:, :count])

        # Create a dictionary to hold arrays of log2_counts for each gene
        log2_counts_by_gene = Dict{String, Vector{Float64}}()
        for row in eachrow(df)
            gene = row[:gene]
            log2_count = row[:log2_count]
            if haskey(log2_counts_by_gene, gene)
                push!(log2_counts_by_gene[gene], log2_count)
            else
                log2_counts_by_gene[gene] = [log2_count]
            end
        end

        sorted_genes = sort_by_median(log2_counts_by_gene)
        labels = map(x->x[1], sorted_genes)
        boxes = map(x->x[2], sorted_genes);

        return boxplot(labels, boxes, title="Amplicon counts per gene", ylabel="log2(count)", xlabel="gene")
    end

end
