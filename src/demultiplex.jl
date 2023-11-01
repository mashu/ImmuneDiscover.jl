module demultiplex
    using CSV
    using FASTX
    using DataFrames
    using Statistics
    include("data.jl")
    using .data

    """
        demux(fastq_path, indices_path; min_length=0, save_fastq_files=false)

    Function to demultiplex plate with indices using double barcoding
    """
    function demux(fastq_path, indices_path; min_length=0, save_fastq_files=false)
        indices = CSV.File(indices_path)
        @assert all(string.(eachindex(first(indices))) .== ["forward_index","reverse_index","case"]) "File must contain following columns forward_index, reverse_index, case"
        @info "Loaded $(length(indices)) indices"

        records = []
        fastq_writers = Dict{String, FASTX.FASTQ.Writer}()

        total = 0
        short = 0

        # Processing callback
        data.process_fastq(fastq_path) do record
            name, genomic_sequence = identifier(record), string(sequence(record))

            @inbounds for i in eachindex(indices)
                forward = indices[i].forward_index
                reverse = indices[i].reverse_index
                case = indices[i].case

                if startswith(genomic_sequence, forward) && endswith(genomic_sequence, reverse)
                    if length(genomic_sequence) < min_length
                        short += 1
                        continue
                    end

                    push!(records, (i, case, name, genomic_sequence))

                    if save_fastq_files
                        dir_path = "$(fastq_path)_split"
                        isdir(dir_path) || mkpath(dir_path)

                        fastq_writer = get!(fastq_writers, case) do
                            io_stream = open(joinpath(dir_path, "$(case).fastq"), "w")
                            FASTX.FASTQ.Writer(io_stream)
                        end

                        fastq_record = FASTX.FASTQ.Record(identifier(record), sequence(record), quality(record))
                        FASTX.FASTQ.write(fastq_writer, fastq_record)
                    end
                end
            end

            total += 1
        end

        # Close all FASTQ writers
        if save_fastq_files
            for writer in values(fastq_writers)
                close(writer)  # `writer.io` should point to the original `io_stream`
            end
        end

        percent = round((length(records)/total), digits=4) * 100
        @info "Demultiplexed $(length(records)) out of $total ($percent)% sequences from FASTQ"
        @info "Dropped $short sequences shorter than $min_length"

        records_df = DataFrame(records)
        rename!(records_df, [:well, :case, :name, :genomic_sequence])

        precision(x) = round(x, digits=1)
        records_df[:, :length] = map(length, records_df[:, :genomic_sequence])
        stats = combine(groupby(records_df, [:well, :case])) do group
            (
                len_μ = precision(mean(group.length)),
                len_σ = precision(std(group.length)),
                len_q1 = precision(quantile(group.length, 0.25)),
                len_q2 = precision(quantile(group.length, 0.5)),
                len_q3 = precision(quantile(group.length, 0.75)),
                len_min = minimum(group.length),
                len_max = maximum(group.length),
                seq_count = nrow(group)
            )
        end
        sort!(stats, :well)

        return records_df[:, [:case, :name, :genomic_sequence]], stats
    end
end