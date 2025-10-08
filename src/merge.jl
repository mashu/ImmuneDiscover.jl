module merge
    using FASTX
    using Logging
    using DataStructures

    export merge_fasta_files, handle_merge

    """
        merge_fasta_files(input_files, output_file; kwargs...)

    Merge multiple FASTA files, keeping only unique sequences.

    # Arguments
    - `input_files::Vector{String}`: Paths to input FASTA files
    - `output_file::String`: Path to the output merged FASTA file

    # Keyword Arguments
    - `sort_by_name::Bool=true`: Sort sequences by name in output
    - `cleanup_pattern::Union{String,Nothing}=nothing`: Regex pattern to remove from sequence names
    - `prefer_first::Bool=true`: When duplicate sequences have different names, prefer the first encountered
    - `add_source_prefix::Bool=false`: Add source filename as prefix to sequence names
    """
    function merge_fasta_files(input_files::Vector{String}, output_file::String;
                              sort_by_name::Bool=true,
                              cleanup_pattern::Union{String,Nothing}=nothing,
                              prefer_first::Bool=true,
                              add_source_prefix::Bool=false)
        
        @info "Merging $(length(input_files)) FASTA files into $output_file"
        
        # Dictionary to store unique sequences: sequence -> (name, source_file)
        unique_sequences = Dict{String, Tuple{String, String}}()
        # Track duplicates: sequence -> list of (name, source_file) that were duplicates
        duplicate_sequences = Dict{String, Vector{Tuple{String, String}}}()
        total_sequences = 0
        
        # Compile cleanup regex once before processing
        cleanup_regex = cleanup_pattern !== nothing ? Regex(cleanup_pattern) : nothing
        
        for (file_idx, input_file) in enumerate(input_files)
            @info "Processing file $file_idx/$(length(input_files)): $input_file"
            
            if !isfile(input_file)
                @warn "File not found: $input_file - skipping"
                continue
            end
            
            file_sequences = 0
            source_name = splitext(basename(input_file))[1]  # filename without extension
            
            open(FASTA.Reader, input_file) do reader
                for record in reader
                    identifier = FASTA.identifier(record)
                    description = FASTA.description(record)
                    # Only add description if it's different from identifier to avoid duplicates
                    if isempty(description) || description == identifier
                        name = identifier
                    else
                        name = "$identifier $description"
                    end
                    sequence = string(FASTA.sequence(record))
                    
                    # Clean up name if cleanup pattern is specified
                    if cleanup_regex !== nothing && occursin(cleanup_regex, name)
                        name = strip(replace(name, cleanup_regex => ""))
                    end
                    
                    # Add source prefix if requested
                    if add_source_prefix
                        name = "$(source_name)_$name"
                    end
                    
                    # Check if sequence already exists
                    if haskey(unique_sequences, sequence)
                        existing_name, existing_source = unique_sequences[sequence]
                        
                        # Track this duplicate
                        if !haskey(duplicate_sequences, sequence)
                            # First duplicate for this sequence - initialize with the original
                            duplicate_sequences[sequence] = [(existing_name, existing_source)]
                        end
                        push!(duplicate_sequences[sequence], (name, input_file))
                        
                        if !prefer_first
                            # Update with current name (prefer last)
                            unique_sequences[sequence] = (name, input_file)
                            @debug "Updated duplicate sequence: $existing_name → $name"
                        else
                            @debug "Skipping duplicate sequence: $name (keeping $existing_name)"
                        end
                    else
                        # New unique sequence
                        unique_sequences[sequence] = (name, input_file)
                        @debug "Added unique sequence: $name"
                    end
                    
                    file_sequences += 1
                    total_sequences += 1
                end
            end
            
            @info "  Loaded $file_sequences sequences from $input_file"
        end
        
        @info "Total sequences processed: $total_sequences"
        @info "Unique sequences found: $(length(unique_sequences))"
        
        # Prepare sequences for output
        output_records = [(name, sequence) for (sequence, (name, source)) in unique_sequences]
        
        # Sort by name if requested
        if sort_by_name
            @info "Sorting sequences by name"
            sort!(output_records, by=x->x[1])
        end
        
        # Write merged FASTA file
        @info "Writing $(length(output_records)) unique sequences to $output_file"
        open(FASTA.Writer, output_file) do writer
            for (name, sequence) in output_records
                record = FASTARecord(name, sequence)
                write(writer, record)
            end
        end
        
        duplicates_removed = total_sequences - length(unique_sequences)
        @info "Successfully merged FASTA files: $duplicates_removed duplicates removed"
        
        # Report detailed duplicate information
        if !isempty(duplicate_sequences)
            @info "Duplicate sequences found:"
            for (sequence, duplicate_list) in duplicate_sequences
                @info "  Sequence with $(length(duplicate_list)) duplicates:"
                for (dup_name, dup_source) in duplicate_list
                    @info "    - $dup_name (from $(basename(dup_source)))"
                end
            end
        end
        
        return length(unique_sequences)
    end

    """
        merge_fasta_files(input_file1, input_file2, output_file; kwargs...)

    Convenience method for merging exactly two FASTA files.
    """
    function merge_fasta_files(input_file1::String, input_file2::String, output_file::String; kwargs...)
        return merge_fasta_files([input_file1, input_file2], output_file; kwargs...)
    end

    """
        handle_merge(parsed_args)

    Handle merge command from CLI arguments
    """
    function handle_merge(parsed_args)
        merge_fasta_files(
            String.(parsed_args["fasta"]["merge"]["inputs"]),
            parsed_args["fasta"]["merge"]["output"];
            sort_by_name = !parsed_args["fasta"]["merge"]["no-sort"],
            cleanup_pattern = parsed_args["fasta"]["merge"]["cleanup"],
            prefer_first = !parsed_args["fasta"]["merge"]["prefer-last"],
            add_source_prefix = parsed_args["fasta"]["merge"]["add-source-prefix"]
        )
    end

end
