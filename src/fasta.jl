module fasta
    using CSV
    using DataFrames
    using FASTX
    using Logging

    export extract_sequences_to_fasta

    """
        extract_sequences_to_fasta(input_file, output_file; kwargs...)

    Extract sequences from a TSV file and save them to a FASTA file.

    # Arguments
    - `input_file::String`: Path to the input TSV file
    - `output_file::String`: Path to the output FASTA file

    # Keyword Arguments
    - `colname::String="allele_name"`: Name of the column with sequence names/IDs
    - `colseq::String="seq"`: Name of the column with sequences
    - `coldesc::Union{String,Nothing}=nothing`: Optional column with descriptions for FASTA headers
    - `filter_pattern::Union{String,Nothing}=nothing`: Optional regex pattern to filter colname column
    - `desc_filter_pattern::Union{String,Nothing}=nothing`: Optional regex pattern to filter coldesc column (use capture groups to include description in FASTA header)
    - `cleanup_pattern::Union{String,Nothing}=nothing`: Optional regex pattern to remove from sequence names (e.g., " Novel")
    - `sort_by_name::Bool=true`: Sort alleles by name before saving
    """
    function extract_sequences_to_fasta(input_file::String, output_file::String;
                                       colname::String="allele_name",
                                       colseq::String="seq", 
                                       coldesc::Union{String,Nothing}=nothing,
                                       filter_pattern::Union{String,Nothing}=nothing,
                                       desc_filter_pattern::Union{String,Nothing}=nothing,
                                       cleanup_pattern::Union{String,Nothing}=nothing,
                                       sort_by_name::Bool=true)
        
        @info "Extracting sequences to FASTA format"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        
        # Assert that the input file has the required columns
        @assert colname ∈ names(df) "Input file must have $colname column"
        @assert colseq ∈ names(df) "Input file must have $colseq column"
        
        # Check if description column is specified and exists
        if coldesc !== nothing
            @assert coldesc ∈ names(df) "Input file must have $coldesc column when coldesc is specified"
            @info "Using $coldesc column for FASTA descriptions"
        end
        
        @info "Loaded $(nrow(df)) rows from input file"
        
        # Apply regex filter if specified
        if filter_pattern !== nothing
            @info "Filtering rows where $colname matches regex pattern '$filter_pattern'"
            df = df[occursin.(Regex(filter_pattern), df[!, colname]), :]
            @info "After filtering: $(nrow(df)) rows remaining"
        end
        
        # Apply description regex filter if specified
        if desc_filter_pattern !== nothing
            if coldesc === nothing
                error("Description filter pattern specified but no coldesc column provided")
            end
            @info "Filtering rows where $coldesc matches regex pattern '$desc_filter_pattern'"
            df = df[occursin.(Regex(desc_filter_pattern), df[!, coldesc]), :]
            @info "After description filtering: $(nrow(df)) rows remaining"
        end
        
        # Remove duplicates based on both name and sequence columns
        before_unique = nrow(df)
        df = DataFrames.unique(df, [colname, colseq])
        after_unique = nrow(df)
        if before_unique != after_unique
            @info "Removed $(before_unique - after_unique) duplicate records (by $colname and $colseq)"
        end
        
        # Sort by name if requested
        if sort_by_name
            @info "Sorting alleles by name"
            df = sort(df, colname)
        end
        
        # Create FASTA records and write to file
        @info "Writing $(nrow(df)) sequences to $output_file"
        writer = FASTA.Writer(open(output_file, "w"))
        
        try
            for row in eachrow(df)
                name = string(row[colname])
                seq = string(row[colseq])
                
                # Clean up name using cleanup pattern (e.g., remove " Novel")
                if cleanup_pattern !== nothing && occursin(Regex(cleanup_pattern), name)
                    original_name = name
                    name = strip(replace(name, Regex(cleanup_pattern) => ""))
                    if name != original_name
                        @debug "Cleaned name: '$original_name' → '$name'"
                    end
                end
                
                # Also clean up name if it contains the filter pattern (legacy behavior)
                if filter_pattern !== nothing && occursin(Regex(filter_pattern), name)
                    name = rstrip(replace(name, Regex(filter_pattern) => ""))
                end
                
                # Add description if specified
                if coldesc !== nothing
                    desc = string(row[coldesc])
                    if !isempty(desc)
                        # If desc_filter_pattern has capture groups, use only the captured part
                        if desc_filter_pattern !== nothing
                            match_result = match(Regex(desc_filter_pattern), desc)
                            if match_result !== nothing && !isempty(match_result.captures) && match_result.captures[1] !== nothing
                                # Use first capture group if available
                                desc = match_result.captures[1]
                                name = "$name $desc"
                            elseif match_result !== nothing
                                # No capture groups, don't include description (just filter)
                                # name stays as is
                            else
                                # Pattern didn't match, skip this row (shouldn't happen due to filtering)
                                continue
                            end
                        else
                            # No filter pattern, include full description
                            name = "$name $desc"
                        end
                    end
                end
                
                # Create and write FASTA record
                record = FASTARecord(name, seq)
                write(writer, record)
            end
        finally
            close(writer)
        end
        
        @info "Successfully extracted $(nrow(df)) unique sequences to $output_file"
        return nothing
    end

end
