module Fasta
    using CSV
    using DataFrames
    using Logging

    export extract_sequences_to_fasta, handle_fasta_diff, handle_fasta_hash

    """
        extract_sequences_to_fasta(input_file, output_file; kwargs...)

    Extract sequences from a TSV file and save them to a FASTA file.

    # Arguments
    - `input_file::String`: Path to the input TSV file
    - `output_file::String`: Path to the output FASTA file

    # Keyword Arguments
    - `colname::String="allele_name"`: Name of the column with sequence names/IDs
    - `colseq::String="seq"`: Name of the column(s) with sequences. Can be a single column name or comma-separated list for concatenation
    - `coldesc::Union{String,Nothing}=nothing`: Optional column with descriptions for FASTA headers
    - `filter_pattern::Union{String,Nothing}=nothing`: Optional regex pattern to filter colname column
    - `desc_filter_pattern::Union{String,Nothing}=nothing`: Optional regex pattern to filter coldesc column (use capture groups to include description in FASTA header)
    - `cleanup_pattern::Union{String,Nothing}=nothing`: Optional regex pattern to remove from sequence names (e.g., " Novel")
    - `sort_by_name::Bool=true`: Sort alleles by name before saving
    - `mincase::Int=1`: Minimum number of donors (cases) that must have this allele to include it in output
    - `case_col::String="case"`: Name of the column containing donor/case identifiers
    - `unique_sequences::Bool=false`: Keep only unique sequences (ignore names, use first encountered name)
    """
    function extract_sequences_to_fasta(input_file::String, output_file::String;
                                       colname::String="allele_name",
                                       colseq::String="seq", 
                                       coldesc::Union{String,Nothing}=nothing,
                                       filter_pattern::Union{String,Nothing}=nothing,
                                       desc_filter_pattern::Union{String,Nothing}=nothing,
                                       cleanup_pattern::Union{String,Nothing}=nothing,
                                       sort_by_name::Bool=true,
                                       mincase::Int=1,
                                       case_col::String="case",
                                       unique_sequences::Bool=false)
        
        @info "Extracting sequences to FASTA format"
        df = CSV.File(input_file, delim='\t') |> DataFrame
        
        # Parse colseq as comma-separated list of columns
        colseq_list = strip.(split(colseq, ','))
        
        # Assert that the input file has the required columns
        @assert colname ∈ names(df) "Input file must have $colname column"
        for seq_col in colseq_list
            @assert seq_col ∈ names(df) "Input file must have $seq_col column"
        end
        
        if length(colseq_list) > 1
            @info "Using multiple sequence columns for concatenation: $(join(colseq_list, ", "))"
        end
        
        # Check if description column is specified and exists
        if coldesc !== nothing
            @assert coldesc ∈ names(df) "Input file must have $coldesc column when coldesc is specified"
            @info "Using $coldesc column for FASTA descriptions"
        end
        
        @info "Loaded $(nrow(df)) rows from input file"
        
        # Apply regex filter if specified
        if filter_pattern !== nothing
            @info "Filtering rows where $colname matches regex pattern '$filter_pattern'"
            filter_regex_prefilter = Regex(filter_pattern)
            df = df[occursin.(filter_regex_prefilter, df[!, colname]), :]
            @info "After filtering: $(nrow(df)) rows remaining"
        end
        
        # Apply description regex filter if specified
        if desc_filter_pattern !== nothing
            if coldesc === nothing
                error("Description filter pattern specified but no coldesc column provided")
            end
            @info "Filtering rows where $coldesc matches regex pattern '$desc_filter_pattern'"
            desc_regex_prefilter = Regex(desc_filter_pattern)
            df = df[occursin.(desc_regex_prefilter, df[!, coldesc]), :]
            @info "After description filtering: $(nrow(df)) rows remaining"
        end
        
        # Apply mincase filter if specified
        if mincase > 1
            @assert case_col ∈ names(df) "Input file must have $case_col column when mincase > 1"
            @info "Filtering alleles present in at least $mincase donors"
            
            # Count unique donors per allele (combination of colname and all colseq columns)
            grouping_cols = [colname; colseq_list]
            allele_donor_counts = combine(groupby(df, grouping_cols), 
                                        case_col => (x -> length(unique(x))) => :donor_count)
            
            # Filter alleles that meet the minimum donor count
            qualifying_alleles = allele_donor_counts[allele_donor_counts.donor_count .>= mincase, grouping_cols]
            
            # Join back to get all rows for qualifying alleles
            df = innerjoin(df, qualifying_alleles, on=grouping_cols)
            
            @info "After mincase filtering (≥$mincase donors): $(nrow(df)) rows remaining"
        end
        
        # Remove duplicates based on both name and sequence columns
        before_unique = nrow(df)
        unique_cols = [colname; colseq_list]
        df = DataFrames.unique(df, unique_cols)
        after_unique = nrow(df)
        if before_unique != after_unique
            @info "Removed $(before_unique - after_unique) duplicate records (by $colname and sequence columns)"
        end
        
        # If unique_sequences is true, keep only the first occurrence of each unique sequence
        if unique_sequences
            @info "Extracting unique sequences (ignoring sequence names)"
            before_seq_unique = nrow(df)
            df = DataFrames.unique(df, colseq_list)
            after_seq_unique = nrow(df)
            if before_seq_unique != after_seq_unique
                @info "Removed $(before_seq_unique - after_seq_unique) duplicate sequences with different names"
            end
        end
        
        # Sort by name if requested
        if sort_by_name
            @info "Sorting alleles by name"
            df = sort(df, colname)
        end
        
        # Create FASTA records and write to file (single-line format)
        @info "Writing $(nrow(df)) sequences to $output_file"
        
        # Compile regexes once before the loop
        cleanup_regex = cleanup_pattern !== nothing ? Regex(cleanup_pattern) : nothing
        filter_regex = filter_pattern !== nothing ? Regex(filter_pattern) : nothing
        desc_filter_regex = desc_filter_pattern !== nothing ? Regex(desc_filter_pattern) : nothing
        
        open(output_file, "w") do io
            for row in eachrow(df)
                name = string(row[colname])
                # Concatenate sequences from all columns in colseq_list
                seq = join([string(row[col]) for col in colseq_list], "")
                
                # Clean up name using cleanup pattern (e.g., remove " Novel")
                if cleanup_regex !== nothing && occursin(cleanup_regex, name)
                    original_name = name
                    name = strip(replace(name, cleanup_regex => ""))
                    if name != original_name
                        @debug "Cleaned name: '$original_name' → '$name'"
                    end
                end
                
                # Also clean up name if it contains the filter pattern (legacy behavior)
                if filter_regex !== nothing && occursin(filter_regex, name)
                    name = rstrip(replace(name, filter_regex => ""))
                end
                
                # Add description if specified
                if coldesc !== nothing
                    desc = string(row[coldesc])
                    if !isempty(desc)
                        # If desc_filter_pattern has capture groups, use only the captured part
                        if desc_filter_regex !== nothing
                            match_result = match(desc_filter_regex, desc)
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
                
                # Write FASTA record in single-line format
                println(io, ">$name")
                println(io, seq)
            end
        end
        
        @info "Successfully extracted $(nrow(df)) unique sequences to $output_file"
        return nothing
    end

    """
        handle_fasta_diff(parsed_args, immunediscover_module)

    Handle fasta diff command from CLI arguments
    """
    function handle_fasta_diff(parsed_args, immunediscover_module)
        fasta_paths = parsed_args["fasta"]["diff"]["fasta"]
        fasta_files = [(file=file, records=immunediscover_module.load_fasta.(file)) for file in fasta_paths]
        sets = [(file=x, set=immunediscover_module.KeyedSet(reverse.(y))) for (x,y) in fasta_files]
        for i in eachindex(sets)
            for j in i+1:length(sets)
                println("Comparing $(sets[i].file) vs $(sets[j].file)")
                println("Union: $(length(union(sets[i].set, sets[j].set)))")
                println("Intersection: $(length(intersect(sets[i].set, sets[j].set)))")
                ab = last.(collect(setdiff(sets[i].set, sets[j].set)))
                ba = last.(collect(setdiff(sets[j].set, sets[i].set)))
                println("Difference ($(length(ab))): \n$(join(ab, "\n"))")
                println("Difference ($(length(ba))): \n$(join(ba, "\n"))")
            end
        end
    end

    """
        handle_fasta_hash(parsed_args, immunediscover_module)

    Handle fasta hash command from CLI arguments
    """
    function handle_fasta_hash(parsed_args, immunediscover_module)
        @info "Hashing alleles"
        fastain = parsed_args["fasta"]["hash"]["fastain"]
        db = immunediscover_module.load_fasta(fastain)
        for (name, seq) in db
            newname = immunediscover_module.unique_name(name, seq)
            println(">$(newname)\n$(seq)")
        end
    end

end
