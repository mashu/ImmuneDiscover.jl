module Merge
    using FASTX
    using Logging
    using ..Data: load_fasta, register_fasta_entry!
    using ..Option: Absent, Present, absent, optional

    export merge_fasta_files, handle_merge

    """
        merge_fasta_files(input_files, output_file; kwargs...)

    Merge multiple FASTA files into one. Every allele name and sequence must be unique across
    all inputs; conflicts abort with an error.
    """
    function merge_fasta_files(input_files::Vector{String}, output_file::String;
                              sort_by_name::Bool=true,
                              cleanup_pattern=absent,
                              add_source_prefix::Bool=false)
        @info "Merging $(length(input_files)) FASTA files into $output_file"
        cleanup = compiled_cleanup(optional(cleanup_pattern))
        name_to_seq = Dict{String,String}()
        seq_to_name = Dict{String,String}()
        output_records = Vector{Tuple{String,String}}()

        for (file_idx, input_file) in enumerate(input_files)
            @info "Processing file $file_idx/$(length(input_files)): $input_file"
            isfile(input_file) || throw(ErrorException("File not found: $input_file"))
            source_name = splitext(basename(input_file))[1]
            file_sequences = 0
            for (name, sequence) in load_fasta(input_file)
                name = apply_cleanup(name, cleanup)
                add_source_prefix && (name = "$(source_name)_$name")
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

    compiled_cleanup(::Absent) = absent
    compiled_cleanup(s::Present) = Present(Regex(s.value))

    apply_cleanup(name, ::Absent) = name
    function apply_cleanup(name, rx::Present{Regex})
        occursin(rx.value, name) || return name
        return strip(replace(name, rx.value => ""))
    end

    function merge_fasta_files(input_file1::String, input_file2::String, output_file::String; kwargs...)
        return merge_fasta_files([input_file1, input_file2], output_file; kwargs...)
    end

    function handle_merge(parsed_args)
        merge_fasta_files(
            String.(parsed_args["fasta"]["merge"]["inputs"]),
            parsed_args["fasta"]["merge"]["output"];
            sort_by_name = !parsed_args["fasta"]["merge"]["no-sort"],
            cleanup_pattern = optional(parsed_args["fasta"]["merge"]["cleanup"]),
            add_source_prefix = parsed_args["fasta"]["merge"]["add-source-prefix"],
        )
    end
end
