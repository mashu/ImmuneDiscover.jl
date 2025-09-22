module exact
    using CSV
    using DataFrames
    using ProgressMeter
    using Folds

    function filter_tuple_by_types(types_as_strings, mandatory_key; tuple_data)
        # Ensure mandatory_key is a valid option
        if mandatory_key ∉ ["prefix", "suffix"]
            error("mandatory_key must be either 'prefix' or 'suffix'")
        end

        # Define the full tuple structure as strings, including both prefix and suffix possibilities
        full_tuple_keys_as_strings = ["prefix", "sequence", "heptamer", "spacer", "nonamer", "suffix"]

        # Initialize selected data with mandatory 'sequence' and either 'prefix' or 'suffix'
        selected_data = Dict("sequence" => tuple_data.sequence)
        if mandatory_key in keys(tuple_data)
            selected_data[mandatory_key] = getfield(tuple_data, Symbol(mandatory_key))
        end

        # Convert the user-specified types to a set of strings for efficient look-up
        types_set = Set(types_as_strings)

        # Loop through each possible key and add it to the selected data if it's included in the types
        for key_str in full_tuple_keys_as_strings
            key_sym = Symbol(key_str) # Convert the string key to a symbol for field access
            if (key_str in types_as_strings || key_str in ["sequence", mandatory_key]) && key_sym in keys(tuple_data)
                # Dynamically select the field from the tuple and add it to the selected data
                selected_data[key_str] = getfield(tuple_data, key_sym)
            end
        end

        # Convert the dictionary back to a tuple for the result
        result_tuple = NamedTuple{Tuple(Symbol.(keys(selected_data)))}(values(selected_data))

        return result_tuple
    end

    """
        extract_flanking(genomic_sequence::String, range::Tuple{Int, Int}, gene_type::String, n::Int)

    Extract flanking sequences from genomic sequence
    """
    function extract_flanking(genomic_sequence::String, range::Tuple{Int, Int}, gene_type::String, n::Int, extension::Union{Int,Nothing}=nothing)
        # Ensure range is valid
        start_pos, end_pos = range
        if start_pos < 1 || end_pos > length(genomic_sequence) || start_pos > end_pos
            error("Invalid range")
        end

        sequence = genomic_sequence[start_pos:end_pos]

        if extension !== nothing
            if gene_type == "V"
                prefix = start_pos > n ? genomic_sequence[(start_pos - n):(start_pos - 1)] : genomic_sequence[1:(start_pos - 1)]
                extension_seq = genomic_sequence[end_pos + 1:min(end_pos + extension, length(genomic_sequence))]
                return (prefix=prefix, sequence=sequence, suffix=extension_seq)
            elseif gene_type == "J"
                extension_seq = genomic_sequence[max(1, start_pos - extension):(start_pos - 1)]
                suffix = genomic_sequence[end_pos + 1:min(end_pos + n, length(genomic_sequence))]
                return (prefix=extension_seq, sequence=sequence, suffix=suffix)
            elseif gene_type == "D"
                prefix = genomic_sequence[max(1, start_pos - extension):(start_pos - 1)]
                suffix = genomic_sequence[end_pos + 1:min(end_pos + extension, length(genomic_sequence))]
                return (prefix=prefix, sequence=sequence, suffix=suffix)
            else
                error("Invalid gene type")
            end
        end

        if gene_type == "V"
            prefix = start_pos > n ? genomic_sequence[(start_pos - n):(start_pos - 1)] : genomic_sequence[1:(start_pos - 1)]
            heptamer = genomic_sequence[end_pos + 1:min(end_pos + 7, length(genomic_sequence))]
            spacer = genomic_sequence[end_pos + 8:min(end_pos + 30, length(genomic_sequence))]
            nonamer = genomic_sequence[end_pos + 31:min(end_pos + 39, length(genomic_sequence))]
            return (prefix=prefix, sequence=sequence, heptamer=heptamer, spacer=spacer, nonamer=nonamer)
        elseif gene_type == "J"
            heptamer = genomic_sequence[max(1, start_pos - 7):(start_pos - 1)]
            spacer = genomic_sequence[max(1, start_pos - 30):max(1, start_pos - 8)]
            nonamer = genomic_sequence[max(1, start_pos - 39):max(1, start_pos - 31)]
            suffix = genomic_sequence[end_pos + 1:min(end_pos + n, length(genomic_sequence))]
            return (nonamer=nonamer, spacer=spacer, heptamer=heptamer, suffix=suffix, sequence=sequence)
        elseif gene_type == "D"
            pre_heptamer = genomic_sequence[max(1, start_pos - 7):(start_pos - 1)]
            pre_spacer = genomic_sequence[max(1, start_pos - 19):max(1, start_pos - 8)]
            pre_nonamer = genomic_sequence[max(1, start_pos - 28):max(1, start_pos - 20)]
            post_heptamer = genomic_sequence[end_pos + 1:min(end_pos + 7, length(genomic_sequence))]
            post_spacer = genomic_sequence[end_pos + 8:min(end_pos + 19, length(genomic_sequence))]
            post_nonamer = genomic_sequence[end_pos + 20:min(end_pos + 28, length(genomic_sequence))]
            return (pre_nonamer=pre_nonamer, pre_spacer=pre_spacer, pre_heptamer=pre_heptamer,sequence=sequence, post_heptamer=post_heptamer, post_spacer=post_spacer, post_nonamer=post_nonamer)
        else
            error("Invalid gene type")
        end
    end

    # Custom function to sort and select top N rows for each group
    function sort_and_select_top_n(df, N)
        sorted_df = sort(df, :flank_count, rev=true)
        return first(sorted_df, N)
    end

    # Function to skip filtering by ratio for certain genes
    function get_ratio(expect_dict, row, ratio)
        if row.db_name in keys(expect_dict)
            @info "Skipping allelic ratio filters for $(row.db_name) in case $(row.case)"
            return 0
        elseif row.gene in keys(expect_dict)
            @info "Skipping gene ratio filters for $(row.db_name) in case $(row.case)"
            return 0
        else
            return ratio # Return ratio if not in expect_dict
        end
    end

    """
        exact_search(table, query, gene; mincount=10, minratio=0.2, collapse=true)

    Exact search per case.
    Parameters:
        table::DataFrame: Table containing genomic sequences
        query::Dict{String, String}: Dictionary containing database name and sequence
        gene::String: Gene type (V, D, J)
        mincount::Int: Minimum number of reads for a sequence to be considered
        minratio::Float: Minimum ratio of reads for a sequence to be considered
        full_mincount::Int: Minimum number of reads for a full record to be considered
        full_minratio::Float: Minimum ratio of reads for a full record to be considered
        affix::Int: Number of bases to extract from the non-RSS side of the sequence
        rss::Vector{String}: Vector of RSS sequence names to filter by
        N::Int: Number of top records to return for each group
    """
    function exact_search(table, query, gene; mincount=10, minratio=0.01, affix=13, rss=["heptamer", "spacer", "nonamer"], extension=nothing, N=10, raw=nothing, expect_dict=Dict{String,Float64}())
        @info "Using mincount: $mincount, minratio: $minratio for genes: $gene"
        @assert all([name in names(table) for name in ["well","case","name","genomic_sequence"]]) "File must contain following columns: well, case, name, genomic_sequence"
        p = Progress(nrow(table))
        result = Folds.map(eachrow(table)) do row
            next!(p)
            case = row.case
            well = row.well
            matches = Vector{NamedTuple}()
            @inbounds for (name, seq) in query
                match = findfirst(seq, row.genomic_sequence)
                if match !== nothing
                    flanks = extract_flanking(row.genomic_sequence, (minimum(match), maximum(match)), gene, affix, extension)
                    local filtered_flanks
                    if extension !== nothing
                        filtered_flanks = flanks
                    else
                        filtered_flanks = gene == "V" ? filter_tuple_by_types(rss, "prefix", tuple_data = flanks) : gene == "J" ? filter_tuple_by_types(rss, "suffix", tuple_data = flanks) : flanks
                    end
                    meta = (well=string(well), case=string(case), db_name=string(name))
                    push!(matches, merge(meta, filtered_flanks))
                end
            end
            matches
        end
        valid_results = [r for r in result if r != []]
        if isempty(valid_results)
            result_df = DataFrame()
        else
            result_df = DataFrame(reduce(vcat, valid_results))
        end

        # Compute counts
        df = transform(groupby(result_df, names(result_df)), nrow => :full_count)
        transform!(groupby(df, [:well, :case, :db_name, :sequence]), nrow => :count)
        transform!(df, :db_name => ByRow(x -> first(split(x, '*'))) => :gene)

        # Save raw dawa
        if raw !== nothing
            CSV.write(raw*".gz", result_df, delim='\t', compress=true)
        end

        # Apply filters
        transform!(groupby(df, [:well, :case, :gene]), :full_count => (x->x./maximum(x)) => :full_ratio)
        transform!(groupby(df, [:well, :case, :gene]), :count => (x->x./maximum(x)) => :ratio)

        sort!(df, [:full_count, :count], rev=[true, true])
        # Collapse individual rows into unique rows before filtering
        udf = sort(unique(df),[:well, :case, :gene])
        # Filter by count and ratio
        # Function get_ratio is used to disable filtering by ratio for certain genes in expect_dict
        filter!(row -> (row.full_count >= mincount) & (row.full_ratio >= get_ratio(expect_dict, row, minratio)), udf)  # Order is important here to filter first on smaller numbers
        filter!(row -> (row.count >= mincount) && (row.ratio >= get_ratio(expect_dict, row, minratio)), udf)

        # Reorder columns
        priority_columns = ["well", "case", "gene", "db_name", "count", "full_count", "ratio", "full_ratio"]
        remaining_columns = setdiff(names(udf), priority_columns)
        udf[:, vcat(priority_columns, remaining_columns)]

        # Grouping by the specified columns
        gdf = groupby(udf, [:well, :case, :gene, :db_name, :sequence])

        # Adding flank_index which enumerates records within each group
        udf_indexed = transform(gdf, :well => (x -> 1:length(x)) => :flank_index)

        # Apply the custom function to each group
        limited_udf = filter(x->x.flank_index <= N, udf_indexed)

        return limited_udf
    end

    function transform_counts(group_df, name; count_col=:count)
        # Find the row where the name matches
        ref_row = filter(row -> startswith(row.db_name, name), group_df)

        # Check if the reference row is found
        ref_count = 1
        well, case = first(map(r->(r.well, r.case), eachrow(unique(group_df, [:well,:case]))))
        if isempty(ref_row)
            @warn "Reference name $name not found in well $well and case $case"
        else
            @info "Applying name $name to well $well and case $case"
            ref_count = ref_row.count
        end
        new_name = "$(count_col)_$(first(split(name,'*')))_ratio"

        # Calculate the ratio and add it as a new column
        group_df[!, new_name] = group_df[:, count_col] ./ ref_count

        return group_df
    end

    function grouped_ratios(counts_df, refgene; count_col=:count)
        transformed = []
        for group in groupby(counts_df, [:well, :case])
            if refgene != ""
                group = transform_counts(group, refgene, count_col=count_col)
            end
            push!(transformed, group)
        end
        transformed_df = reduce(vcat, transformed)
        return transformed_df
    end

    export grouped_ratios, transform_counts
end
