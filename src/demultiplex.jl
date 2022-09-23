module demultiplex
    using FASTX
    using CSV
    using DataFrames
    using ProgressMeter
    using CodecZlib
    using Statistics

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