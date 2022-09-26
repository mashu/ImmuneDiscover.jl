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

    """
    Function to demultiplex plate with indices using double barcoding
    """
    function demux(fastq_path, indices_path)
        indices = CSV.File(indices_path) |> DataFrame
        @assert all(names(indices) .== ["forward_index","reverse_index","case"]) "File must contain following columns forward_index, reverse_index, case"
        @info "Loaded $(nrow(indices)) indices"
        data = []
        total = 0
        process_fastq(fastq_path) do record
            name, genomic_sequence = identifier(record), string(sequence(record))
            @inbounds for i in 1:nrow(indices)
                forward = indices[i,:forward_index]
                reverse = indices[i,:reverse_index]
                case = indices[i,:case]
                if startswith(genomic_sequence, forward) & endswith(genomic_sequence, reverse)
                    push!(data, (i, case, name, genomic_sequence))
                end
            end
            # Count number of reads
            total += 1
        end
        percent = round((length(data)/total),digits=4)*100
        @info "Demultiplexed $(length(data)) out of $total ($percent)% sequences from FASTQ"
        data_df = DataFrame(data)
        rename!(data_df, [:well, :case, :name, :genomic_sequence])
        # Display overall per well stats
        @info "Statistics"
        data_df[:,:length] = map(length, data_df[:,:genomic_sequence])
        stats = combine(groupby(data_df, [:well, :case])) do group
            (len_μ = mean(group.length), len_σ = std(group.length), len_q1=quantile(group.length,0.25), len_q2=quantile(group.length,0.5), len_q3=quantile(group.length,0.75), len_min = minimum(group.length), len_max=maximum(group.length), seq_count=nrow(group))
        end
        sort!(stats, :well)
        @info(stats)
        return data_df[:,[:case,:name,:genomic_sequence]]
    end
end