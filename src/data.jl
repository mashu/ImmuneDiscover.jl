module data
    using FASTX
    using TranscodingStreams, CodecZlib
    using ProgressMeter
    using CodecZlib

    """
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
    Save records with sequences and quality scores to fastq file
    """
    function write_gz_fastq(path, records)
        FASTQ.Writer(TranscodingStream(GzipCompressor(),open(path, "w"))) do writer
            for record in records
                write(writer, record)
            end
        end
    end

    """
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
    Save records with sequences to fasta file
    """
    function write_gz_fasta(path, records)
        FASTA.Writer(TranscodingStream(GzipCompressor(),open(path, "w"))) do writer
            for record in records
                write(writer, record)
            end
        end
    end

    """
    Function to load regular FASTA file
    """
    function load_fasta(path)
        records = []
        open(FASTA.Reader, path) do reader
            for record in reader
                push!(records, (FASTA.identifier(record), string(FASTA.sequence(record))))
            end
        end
        return records
    end

    """
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