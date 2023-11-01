module data
    using FASTX
    using CodecZlib
    using ProgressMeter
    using CodecZlib
    using DataFrames
    using Statistics
    using UnicodePlots
    export load_fasta, plotgenes

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
        load_fasta(path)

    Function to load gz compressed or uncompressed FASTA file
    """
    function load_fasta(path)
        records = Vector{Tuple{String,String}}()
        open(FASTA.Reader, path) do reader
            for record in reader
                push!(records, (FASTA.description(record), string(FASTA.sequence(record))))
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
        sort!(median_tuples, by = x -> x[3])

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
