module demultiplex
    using CSV
    using FASTX
    using DataFrames
    using Statistics
    include("data.jl")
    using .data

    """
    Function to demultiplex plate with indices using double barcoding
    """
    function demux(fastq_path, indices_path)
        indices = CSV.File(indices_path) |> DataFrame
        @assert all(names(indices) .== ["forward_index","reverse_index","case"]) "File must contain following columns forward_index, reverse_index, case"
        @info "Loaded $(nrow(indices)) indices"
        records = []
        total = 0
        # Processing callback
        data.process_fastq(fastq_path) do record
            name, genomic_sequence = identifier(record), string(sequence(record))
            @inbounds for i in 1:nrow(indices)
                forward = indices[i,:forward_index]
                reverse = indices[i,:reverse_index]
                case = indices[i,:case]
                if startswith(genomic_sequence, forward) & endswith(genomic_sequence, reverse)
                    push!(records, (i, case, name, genomic_sequence))
                end
            end
            # Count number of reads
            total += 1
        end
        percent = round((length(records)/total),digits=4)*100
        @info "Demultiplexed $(length(records)) out of $total ($percent)% sequences from FASTQ"
        records_df = DataFrame(records)
        rename!(records_df, [:well, :case, :name, :genomic_sequence])
        # Display overall per well stats
        @info "Statistics"
        records_df[:,:length] = map(length, records_df[:,:genomic_sequence])
        stats = combine(groupby(records_df, [:well, :case])) do group
            (len_μ = mean(group.length), len_σ = std(group.length), len_q1=quantile(group.length,0.25), len_q2=quantile(group.length,0.5), len_q3=quantile(group.length,0.75), len_min = minimum(group.length), len_max=maximum(group.length), seq_count=nrow(group))
        end
        sort!(stats, :well)
        @info(stats)
        return records_df[:,[:case,:name,:genomic_sequence]]
    end
end