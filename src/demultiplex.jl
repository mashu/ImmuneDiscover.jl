module demultiplex
    using CSV
    using FASTX
    using DataFrames
    using Statistics
    using UnicodePlots
    using Logging

    # Parent module provides shared IO and validation
    using ..data

    export demux, handle_demultiplex

    function demux(fastq_path, indices_path, array_index=""; min_length=0, save_fastq_files=false)
        indices = CSV.File(indices_path)
        @assert all(string.(eachindex(first(indices))) .== ["forward_index","reverse_index","case"]) "File must contain following columns forward_index, reverse_index, case"
        @info "Loaded $(length(indices)) indices"
        records = Vector{Tuple{Int, String, String, String}}()
        fastq_writers = Dict{String, FASTX.FASTQ.Writer}()

        total = 0
        short = 0

        pattern = if !isempty(array_index)
            [(forward = Regex(string("^.{2}", array_index, ".{24}", i.forward_index)),
              reverse = Regex(i.reverse_index),
              case = i.case,) for i in indices]
        else
            [(forward = Regex(i.forward_index),
              reverse = Regex(i.reverse_index),
              case = i.case,) for i in indices]
        end

        data.process_fastq(fastq_path) do record
            name, genomic_sequence = identifier(record), string(sequence(record))

            @inbounds for i in eachindex(pattern)
                if startswith(genomic_sequence, pattern[i].forward) && endswith(genomic_sequence, pattern[i].reverse)
                    if length(genomic_sequence) < min_length
                        short += 1
                        continue
                    end

                    push!(records, (i, String(pattern[i].case), String(name), String(genomic_sequence)))

                    if save_fastq_files
                        dir_path = "$(fastq_path)_split"
                        isdir(dir_path) || mkpath(dir_path)
                        fastq_writer = get!(fastq_writers, pattern[i].case) do
                            FASTQ.Writer(open(joinpath(dir_path, "$(pattern[i].case).fastq"), "w"))
                        end
                        write(fastq_writer, record)
                    end
                    break
                end
            end
            total += 1
        end

        for writer in values(fastq_writers)
            close(writer)
        end

        if isempty(records)
            @warn "No reads were demultiplexed. Check your indices and FASTQ file."
            return DataFrame(well=Int[], case=String[], name=String[], genomic_sequence=String[]),
                   DataFrame(well=Int[], total=Int[], short=Int[])
        end

        records_df = DataFrame(records)
        rename!(records_df, [:well, :case, :name, :genomic_sequence])

        stats = combine(groupby(records_df, :well),
                        nrow => :count,
                        :genomic_sequence => (x -> round(mean(length.(x)), digits=1)) => :mean_length,
                        :genomic_sequence => (x -> round(std(length.(x)), digits=1)) => :std_length)
        sort!(stats, :well)

        return records_df[:, [:well, :case, :name, :genomic_sequence]], stats
    end

    function handle_demultiplex(parsed_args, always_gz)
        @info "Demultiplexing"
        table, stats = demux(
            parsed_args["demultiplex"]["fastq"],
            parsed_args["demultiplex"]["indices"],
            parsed_args["demultiplex"]["forwardarrayindex"],
            min_length=parsed_args["demultiplex"]["length"],
            save_fastq_files=parsed_args["demultiplex"]["split"]
        )

        case_filter_regex = get(parsed_args["demultiplex"], "case-filter-regex", nothing)
        if case_filter_regex !== nothing
            @info "Filtering cases with regex: $case_filter_regex"
            original_count = nrow(table)
            compiled_regex = Regex(case_filter_regex)
            filter!(x -> startswith(x.case, compiled_regex), table)
            @info "Filtered from $original_count to $(nrow(table)) rows"
        end

        @info "Demultiplexing summary:"
        @info "  - Total sequences: $(nrow(table))"
        @info "  - Unique cases: $(length(unique(table.case)))"
        @info "  - Unique wells: $(length(unique(table.well)))"

        logfile = "$(parsed_args["demultiplex"]["output"]).log"
        CSV.write(logfile, stats, delim='\t')
        count_df = combine(groupby(table, :case), :case => length => :count)
        sort!(count_df, :count, rev=true)
        println(barplot(count_df.case, count_df.count))
        output = always_gz(parsed_args["demultiplex"]["output"])
        CSV.write(output, table, compress=true, delim='\t')
        @info "Demultiplexing data saved in compressed $output file"
    end
end
