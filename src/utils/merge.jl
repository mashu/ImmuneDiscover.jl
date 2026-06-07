module Merge
    using FASTX
    using Logging
    using ..Data: load_fasta, register_fasta_entry!

    export merge_fasta_files, handle_merge

    """
        merge_fasta_files(input_files, output_file; kwargs...)

    Merge multiple FASTA files into one. Every allele name and sequence must be unique across
    all inputs; conflicts abort with an error.

    # Arguments
    - `input_files::Vector{String}`: Paths to input FASTA files
    - `output_file::String`: Path to the output merged FASTA file

    # Keyword Arguments
    - `sort_by_name::Bool=true`: Sort sequences by name in output
    - `cleanup_pattern::Union{String,Nothing}=nothing`: Regex pattern to remove from sequence names
    - `add_source_prefix::Bool=false`: Add source filename as prefix to sequence names
    """
    function merge_fasta_files(input_files::Vector{String}, output_file::String;
                              sort_by_name::Bool=true,
                              cleanup_pattern::Union{String,Nothing}=nothing,
                              add_source_prefix::Bool=false)
        @info "Merging $(length(input_files)) FASTA files into $output_file"

        cleanup_regex = cleanup_pattern !== nothing ? Regex(cleanup_pattern) : nothing
        name_to_seq = Dict{String,String}()
        seq_to_name = Dict{String,String}()
        output_records = Vector{Tuple{String,String}}()

        for (file_idx, input_file) in enumerate(input_files)
            @info "Processing file $file_idx/$(length(input_files)): $input_file"
            isfile(input_file) || throw(ErrorException("File not found: $input_file"))

            source_name = splitext(basename(input_file))[1]
            file_sequences = 0

            for (name, sequence) in load_fasta(input_file)
                if cleanup_regex !== nothing && occursin(cleanup_regex, name)
                    name = strip(replace(name, cleanup_regex => ""))
                end
                if add_source_prefix
                    name = "$(source_name)_$name"
                end
                register_fasta_entry!(name, sequence, input_file, name_to_seq, seq_to_name)
                push!(output_records, (name, sequence))
                file_sequences += 1
            end

            @info "  Loaded $file_sequences sequences from $input_file"
        end

        sort_by_name && sort!(output_records, by=first)

        @info "Writing $(length(output_records)) sequences to $output_file"
        open(FASTA.Writer, output_file) do writer
            for (name, sequence) in output_records
                write(writer, FASTARecord(name, sequence))
            end
        end

        @info "Successfully merged FASTA files"
        return length(output_records)
    end

    function merge_fasta_files(input_file1::String, input_file2::String, output_file::String; kwargs...)
        return merge_fasta_files([input_file1, input_file2], output_file; kwargs...)
    end

    function handle_merge(parsed_args)
        merge_fasta_files(
            String.(parsed_args["fasta"]["merge"]["inputs"]),
            parsed_args["fasta"]["merge"]["output"];
            sort_by_name = !parsed_args["fasta"]["merge"]["no-sort"],
            cleanup_pattern = parsed_args["fasta"]["merge"]["cleanup"],
            add_source_prefix = parsed_args["fasta"]["merge"]["add-source-prefix"],
        )
    end

end
