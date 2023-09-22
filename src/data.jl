module data
    using FASTX
    using CodecZlib
    using ProgressMeter
    using CodecZlib

    export load_fasta

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
        prog = ProgressUnknown("Processing FASTQ reads")
        for record in reader
            f(record)
            ProgressMeter.next!(prog)
        end
        ProgressMeter.finish!(prog)
        close(reader)
    end
end
