module exact
    using CSV
    using DataFrames
    using ProgressMeter
    using Folds
    using FASTX
    using UnicodePlots
    using Statistics
    
    # Typed rows for border statistics (row tables)
    const BORDER_ROW = NamedTuple{(:case, :gene, :matched_total, :accepted_total, :rejected_border, :rejected_ratio), Tuple{String, String, Int, Int, Int, Float64}}
    const BORDER_GENE_ROW = NamedTuple{(:gene, :mean_rejected_ratio, :matched_total, :rejected_border, :num_donors), Tuple{String, Float64, Int, Int, Int}}

    # Stores the most recent border-overlap statistics produced by exact_search
    const LAST_BORDER_STATS = Ref(Vector{BORDER_ROW}())
    const LAST_BORDER_GENE_STATS = Ref(Vector{BORDER_GENE_ROW}())
    
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

    # Merge helper for per-(gene) minimum integer values
    function merge_min!(dst::Dict{String,Int}, src::Dict{String,Int})
        for (k, v) in src
            if haskey(dst, k)
                dst[k] = min(dst[k], v)
            else
                dst[k] = v
            end
        end
        return dst
    end

    # Merge helper for per-(case,gene) integer counters
    function merge_counts!(dst::Dict{Tuple{String,String},Int}, src::Dict{Tuple{String,String},Int})
        for (k, v) in src
            dst[k] = get(dst, k, 0) + v
        end
        return dst
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

    # Overload supporting per-side extension overrides (prefix/suffix independently)
    function extract_flanking(genomic_sequence::String, range::Tuple{Int, Int}, gene_type::String, n::Int,
                              extension::Union{Int,Nothing}, prefix_ext::Union{Int,Nothing}, suffix_ext::Union{Int,Nothing})
        start_pos, end_pos = range
        if start_pos < 1 || end_pos > length(genomic_sequence) || start_pos > end_pos
            error("Invalid range")
        end

        sequence = genomic_sequence[start_pos:end_pos]

        if extension === nothing
            # Fall back to RSS extraction if extension not provided
            return extract_flanking(genomic_sequence, range, gene_type, n, extension)
        end

        # Effective per-side extensions (fallback to global extension if side-specific is nothing)
        left_ext = prefix_ext === nothing ? extension : prefix_ext
        right_ext = suffix_ext === nothing ? extension : suffix_ext

        if gene_type == "V"
            prefix = start_pos > n ? genomic_sequence[(start_pos - n):(start_pos - 1)] : genomic_sequence[1:(start_pos - 1)]
            extension_seq = genomic_sequence[end_pos + 1:min(end_pos + right_ext, length(genomic_sequence))]
            return (prefix=prefix, sequence=sequence, suffix=extension_seq)
        elseif gene_type == "J"
            extension_seq = genomic_sequence[max(1, start_pos - left_ext):(start_pos - 1)]
            suffix = genomic_sequence[end_pos + 1:min(end_pos + n, length(genomic_sequence))]
            return (prefix=extension_seq, sequence=sequence, suffix=suffix)
        elseif gene_type == "D"
            prefix = genomic_sequence[max(1, start_pos - left_ext):(start_pos - 1)]
            suffix = genomic_sequence[end_pos + 1:min(end_pos + right_ext, length(genomic_sequence))]
            return (prefix=prefix, sequence=sequence, suffix=suffix)
        else
            error("Invalid gene type")
        end
    end

    # Determine whether the requested extension overlaps read-end border regions
    function extension_overlaps_border(start_pos::Int, end_pos::Int, read_length::Int, gene_type::String, extension::Int, border::Int)
        if border <= 0 || extension <= 0
            return false
        end

        left_border_end = min(border, read_length)
        right_border_start = max(1, read_length - border + 1)

        if gene_type == "V"
            ext_start = end_pos + 1
            ext_end = min(end_pos + extension, read_length)
            return (ext_start <= ext_end) && (ext_end >= right_border_start)
        elseif gene_type == "J"
            ext_start = max(1, start_pos - extension)
            ext_end = start_pos - 1
            return (ext_start <= ext_end) && (ext_start <= left_border_end)
        elseif gene_type == "D"
            left_ext_start = max(1, start_pos - extension)
            left_ext_end = start_pos - 1
            right_ext_start = end_pos + 1
            right_ext_end = min(end_pos + extension, read_length)
            left_overlap = (left_ext_start <= left_ext_end) && (left_ext_start <= left_border_end)
            right_overlap = (right_ext_start <= right_ext_end) && (right_ext_end >= right_border_start)
            return left_overlap || right_overlap
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
    function exact_search(table, query, gene; mincount=10, minratio=0.01, affix=13, rss=["heptamer", "spacer", "nonamer"], extension=nothing, N=10, raw=nothing, expect_dict=Dict{String,Float64}(), sequence_lookup=nothing, border::Int=0, adjust_per_gene_extension::Bool=false, adjust_percent::Float64=1.0)
        @info "Using mincount: $mincount, minratio: $minratio for genes: $gene"
        @assert all([name in names(table) for name in ["well","case","name","genomic_sequence"]]) "File must contain following columns: well, case, name, genomic_sequence"
        # Optional pre-pass: calibrate per-gene per-side extension to avoid border overlaps
        per_gene_prefix = Dict{String,Int}()
        per_gene_suffix = Dict{String,Int}()
        if extension !== nothing && border > 0 && adjust_per_gene_extension
            @info "Calibrating per-gene extension (prefix/suffix) targeting ≥ $(Int(round(adjust_percent*100)))% safe (no-border) reads"
            # Collect allowed lengths per gene and side
            p_cal = Progress(nrow(table))
            tmp = Folds.map(eachrow(table)) do row
                next!(p_cal)
                local pre = Dict{String,Vector{Int}}()
                local suf = Dict{String,Vector{Int}}()
                read_length = length(row.genomic_sequence)
                left_border_end = min(border, read_length)
                right_border_start = max(1, read_length - border + 1)
                @inbounds for (name, seq) in query
                    match = findfirst(seq, row.genomic_sequence)
                    if match !== nothing
                        start_pos = minimum(match)
                        end_pos = maximum(match)
                        gene_base = first(split(string(name), '*'))
                        if gene == "V"
                            allowed_suf = max(0, right_border_start - 1 - end_pos)
                            vec = get!(suf, gene_base, Int[])
                            push!(vec, min(extension, allowed_suf))
                        elseif gene == "J"
                            allowed_pre = max(0, start_pos - 1 - left_border_end)
                            vec = get!(pre, gene_base, Int[])
                            push!(vec, min(extension, allowed_pre))
                        elseif gene == "D"
                            allowed_pre = max(0, start_pos - 1 - left_border_end)
                            allowed_suf = max(0, right_border_start - 1 - end_pos)
                            vecp = get!(pre, gene_base, Int[])
                            vecs = get!(suf, gene_base, Int[])
                            push!(vecp, min(extension, allowed_pre))
                            push!(vecs, min(extension, allowed_suf))
                        else
                            error("Invalid gene type")
                        end
                    end
                end
                (pre=pre, suf=suf)
            end
            # Merge per-row dictionaries concatenating vectors
            pre_vals = Dict{String,Vector{Int}}()
            suf_vals = Dict{String,Vector{Int}}()
            for r in tmp
                for (g, v) in r.pre
                    if haskey(pre_vals, g)
                        append!(pre_vals[g], v)
                    else
                        pre_vals[g] = copy(v)
                    end
                end
                for (g, v) in r.suf
                    if haskey(suf_vals, g)
                        append!(suf_vals[g], v)
                    else
                        suf_vals[g] = copy(v)
                    end
                end
            end
            # Helper to compute (1 - p) quantile via order statistic
            function safe_quantile(values::Vector{Int}, p::Float64)
                if isempty(values)
                    return extension
                end
                n = length(values)
                sorted = sort(values)
                # index for (1 - p) quantile, ensure within [1,n]
                idx = max(1, min(n, ceil(Int, (1 - p) * n)))
                return sorted[idx]
            end
            # Compute per-gene chosen extensions (cap by requested extension implicitly already applied)
            if gene == "V"
                for (g, arr) in suf_vals
                    per_gene_suffix[g] = safe_quantile(arr, adjust_percent)
                end
            elseif gene == "J"
                for (g, arr) in pre_vals
                    per_gene_prefix[g] = safe_quantile(arr, adjust_percent)
                end
            elseif gene == "D"
                for (g, arr) in pre_vals
                    per_gene_prefix[g] = safe_quantile(arr, adjust_percent)
                end
                for (g, arr) in suf_vals
                    per_gene_suffix[g] = safe_quantile(arr, adjust_percent)
                end
            else
                error("Invalid gene type")
            end
            @info "Per-gene extension calibrated for $(length(union(collect(keys(per_gene_prefix)), collect(keys(per_gene_suffix))))) genes"
        end

        p = Progress(nrow(table))
        result = Folds.map(eachrow(table)) do row
            next!(p)
            case = row.case
            well = row.well
            matches = Vector{NamedTuple}()
            totals = Dict{Tuple{String,String},Int}()  # (case,gene) => matched_total
            accepted = Dict{Tuple{String,String},Int}()  # (case,gene) => accepted_total
            @inbounds for (name, seq) in query
                match = findfirst(seq, row.genomic_sequence)
                if match !== nothing
                    start_pos = minimum(match)
                    end_pos = maximum(match)

                    # Track total potential matches only in extension mode with border logic
                    if extension !== nothing && border > 0
                        gene_base = first(split(string(name), '*'))
                        key = (string(case), gene_base)
                        totals[key] = get(totals, key, 0) + 1
                        # Apply border rejection: either global extension or per-gene side-specific calibrated extension
                        if adjust_per_gene_extension
                            read_length = length(row.genomic_sequence)
                            left_border_end = min(border, read_length)
                            right_border_start = max(1, read_length - border + 1)
                            if gene == "V"
                                suffE = get(per_gene_suffix, gene_base, extension)
                                right_ext_end = min(end_pos + suffE, read_length)
                                if right_ext_end >= right_border_start
                                    continue
                                end
                            elseif gene == "J"
                                prefE = get(per_gene_prefix, gene_base, extension)
                                left_ext_start = max(1, start_pos - prefE)
                                if left_ext_start <= left_border_end
                                    continue
                                end
                            elseif gene == "D"
                                prefE = get(per_gene_prefix, gene_base, extension)
                                suffE = get(per_gene_suffix, gene_base, extension)
                                left_ext_start = max(1, start_pos - prefE)
                                right_ext_end = min(end_pos + suffE, read_length)
                                if (left_ext_start <= left_border_end) || (right_ext_end >= right_border_start)
                                    continue
                                end
                            else
                                error("Invalid gene type")
                            end
                        else
                            if extension_overlaps_border(start_pos, end_pos, length(row.genomic_sequence), gene, extension, border)
                                # Reject this read due to border overlap
                                continue
                            end
                        end
                        accepted[key] = get(accepted, key, 0) + 1
                    end

                    # Use per-gene per-side extension if available
                    local flanks
                    if extension !== nothing && adjust_per_gene_extension && (border > 0)
                        gene_base = first(split(string(name), '*'))
                        prefE = get(per_gene_prefix, gene_base, extension)
                        suffE = get(per_gene_suffix, gene_base, extension)
                        flanks = extract_flanking(row.genomic_sequence, (start_pos, end_pos), gene, affix, extension, prefE, suffE)
                    else
                        flanks = extract_flanking(row.genomic_sequence, (start_pos, end_pos), gene, affix, extension)
                    end
                    local filtered_flanks
                    if extension !== nothing
                        filtered_flanks = flanks
                    else
                        filtered_flanks = gene == "V" ? filter_tuple_by_types(rss, "prefix", tuple_data = flanks) : gene == "J" ? filter_tuple_by_types(rss, "suffix", tuple_data = flanks) : flanks
                    end
                    meta = (well=string(well), case=string(case), db_name=string(name))
                    if extension !== nothing
                        # Include lengths of prefix and suffix for transparency in extension mode
                        local prefix_len = length(filtered_flanks.prefix)
                        local suffix_len = length(filtered_flanks.suffix)
                        push!(matches, merge(meta, merge(filtered_flanks, (prefix_len=prefix_len, suffix_len=suffix_len))))
                    else
                        push!(matches, merge(meta, filtered_flanks))
                    end
                end
            end
            (matches=matches, totals=totals, accepted=accepted)
        end
        # Gather accepted matches
        valid_match_lists = [r.matches for r in result if r.matches != []]
        if isempty(valid_match_lists)
            result_df = DataFrame()
        else
            result_df = DataFrame(reduce(vcat, valid_match_lists))
        end

        # Compute border-overlap statistics when applicable
        if extension !== nothing && border > 0
            # Combine per-row dictionaries
            totals_all = Dict{Tuple{String,String},Int}()
            accepted_all = Dict{Tuple{String,String},Int}()
            for r in result
                if !isempty(r.totals)
                    merge_counts!(totals_all, r.totals)
                end
                if !isempty(r.accepted)
                    merge_counts!(accepted_all, r.accepted)
                end
            end

            if isempty(totals_all)
                LAST_BORDER_STATS[] = BORDER_ROW[]
                LAST_BORDER_GENE_STATS[] = BORDER_GENE_ROW[]
            else
                # Build typed rows from dictionaries
                stats_rows = BORDER_ROW[]
                total_matched = 0
                total_rejected = 0
                for (c, g) in union(collect(keys(totals_all)), collect(keys(accepted_all)))
                    m = get(totals_all, (c, g), 0)
                    a = get(accepted_all, (c, g), 0)
                    r = m - a
                    ratio = m > 0 ? r / m : 0.0
                    push!(stats_rows, (case=c, gene=g, matched_total=m, accepted_total=a, rejected_border=r, rejected_ratio=ratio))
                    total_matched += m
                    total_rejected += r
                end

                total_pct = total_matched > 0 ? round(100 * total_rejected / total_matched; digits=2) : 0.0
                @info "Border filter rejected $(total_rejected) of $(total_matched) potential matches ($(total_pct)%)"

                LAST_BORDER_STATS[] = stats_rows

                # Per-gene summary averaged across donors
                # Aggregate via Dicts to avoid DataFrames
                sum_ratio = Dict{String,Float64}()
                sum_matched = Dict{String,Int}()
                sum_rejected = Dict{String,Int}()
                num_donors = Dict{String,Int}()
                for row in stats_rows
                    g = row.gene
                    sum_ratio[g] = get(sum_ratio, g, 0.0) + row.rejected_ratio
                    sum_matched[g] = get(sum_matched, g, 0) + row.matched_total
                    sum_rejected[g] = get(sum_rejected, g, 0) + row.rejected_border
                    num_donors[g] = get(num_donors, g, 0) + 1
                end
                gene_rows = BORDER_GENE_ROW[]
                for g in keys(sum_ratio)
                    mean_r = num_donors[g] > 0 ? sum_ratio[g] / num_donors[g] : 0.0
                    push!(gene_rows, (gene=g, mean_rejected_ratio=mean_r, matched_total=sum_matched[g], rejected_border=sum_rejected[g], num_donors=num_donors[g]))
                end
                LAST_BORDER_GENE_STATS[] = gene_rows
            end
        else
            LAST_BORDER_STATS[] = BORDER_ROW[]
            LAST_BORDER_GENE_STATS[] = BORDER_GENE_ROW[]
        end

        # Compute counts
        df = transform(groupby(result_df, names(result_df)), nrow => :full_count)
        transform!(groupby(df, [:well, :case, :db_name, :sequence]), nrow => :count)
        transform!(df, :db_name => ByRow(x -> first(split(x, '*'))) => :gene)
        
        # Add isin_db column if sequence_lookup is provided
        if sequence_lookup !== nothing
            @info "Adding isin_db column based on reference FASTA"
            df[!, :isin_db] = map(row -> get(sequence_lookup, row.sequence, false) ? "" : "Novel", eachrow(df))
        end

        # Save raw dawa
        if raw !== nothing
            CSV.write(raw*".gz", result_df, delim='\t', compress=true)
        end

        # Apply filters
        transform!(groupby(df, [:well, :case, :gene]), :full_count => (x->x./maximum(x)) => :full_ratio)
        transform!(groupby(df, [:well, :case, :gene]), :count => (x->x./maximum(x)) => :ratio)

        sort!(df, [:full_count, :count], rev=[true, true])
        # Collapse individual rows into unique rows before filtering
        # Keep unique combinations that preserve flanking sequence differences
        udf = sort(unique(df),[:well, :case, :gene, :db_name, :sequence])
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

    """
        build_sequence_lookup(ref_fasta_path::String) -> Dict{String, Bool}

    Build a dictionary mapping sequences to true if they exist in the reference FASTA file.
    This is used to add an isin_db column to exact search results where sequences in the 
    database get an empty string and novel sequences get "Novel".

    # Arguments
    - `ref_fasta_path::String`: Path to the reference FASTA file

    # Returns
    - `Dict{String, Bool}`: Dictionary mapping sequences to true
    """
    function build_sequence_lookup(ref_fasta_path::String)
        sequence_lookup = Dict{String, Bool}()
        
        @info "Building sequence lookup from reference FASTA: $ref_fasta_path"
        
        open(FASTA.Reader, ref_fasta_path) do reader
            for record in reader
                sequence = string(FASTA.sequence(record))
                sequence_lookup[sequence] = true
            end
        end
        
        @info "Loaded $(length(sequence_lookup)) sequences from reference FASTA"
        return sequence_lookup
    end

    """
        handle_exact(parsed_args, immunediscover_module, always_gz)

    Handle exact search command from CLI arguments
    """
    function handle_exact(parsed_args, immunediscover_module, always_gz)
        
        @info "Exact search"
        extension = parsed_args["search"]["exact"]["extension"]
        border = get(parsed_args["search"]["exact"], "border", 0)
        adjust_per_gene_extension = get(parsed_args["search"]["exact"], "adjust-per-gene-extension", false)
        adjust_percent = get(parsed_args["search"]["exact"], "adjust-percent", 1.0)
        limit = parsed_args["search"]["exact"]["limit"]
        refgenes = parsed_args["search"]["exact"]["refgene"]
        if length(refgenes) > 0
            @info "Using reference genes $refgenes"
        end
        
        table = immunediscover_module.load_demultiplex(parsed_args["search"]["exact"]["tsv"])
        if limit > 0
            @info "Limiting number of demultiplexed reads to $limit"
            table = table[1:limit,:]
        end
        
        db = immunediscover_module.load_fasta(parsed_args["search"]["exact"]["fasta"], validate=false)
        mincount = parsed_args["search"]["exact"]["mincount"]
        minratio = parsed_args["search"]["exact"]["minratio"]
        
        TRUST_MINCOUNT = 5
        if mincount < TRUST_MINCOUNT
            @warn "Decreasing mincount below $TRUST_MINCOUNT may lead to false positives due to sequencing errors - you've been warned!"
        end
        
        top = parsed_args["search"]["exact"]["top"]
        affix = parsed_args["search"]["exact"]["affix"]
        
        local rss
        if extension !== nothing
            @info "Using extension mode with length $extension instead of RSS elements"
            rss = String[]
        else
            rss = split(parsed_args["search"]["exact"]["rss"], ',')
            immunediscover_module.validate_types(rss)
            @info "Extract RSS: $(join(rss,','))"
        end
        
        if top != 1
            @info "Uncollapsed mode enabled; at most $top full records will be returned."
        end
        
        gene = parsed_args["search"]["exact"]["gene"]
        expect = parsed_args["search"]["exact"]["expect"]
        deletion = parsed_args["search"]["exact"]["deletion"]
        min_allele_case_medratio = parsed_args["search"]["exact"]["min-allele-mratio"]
        min_gene_case_medratio = parsed_args["search"]["exact"]["min-gene-mratio"]
        
        # Load expect if provided
        expect_df = DataFrame(name=[], ratio=[])
        if expect !== nothing
            expect_df = CSV.read(expect, DataFrame, delim='\t')
            @assert all([name in names(expect_df) for name in ["name","ratio"]]) "File must contain following columns: name, ratio"
        end
        if nrow(expect_df) > 0
            @info "Using expect file with $(nrow(expect_df)) entries"
        end
        expect_dict = Dict(zip(expect_df.name, expect_df.ratio))
        
        # Load deletions if provided
        deletion_df = DataFrame(name=[], ratio=[])
        if deletion !== nothing
            deletion_df = CSV.read(deletion, DataFrame, delim='\t')
            @assert all([name in names(deletion_df) for name in ["name","ratio"]]) "File must contain following columns: name, ratio"
        end
        if nrow(deletion_df) > 0
            @info "Using deletion_df file with $(nrow(expect_df)) entries"
        end
        deletion_dict = Dict(zip(deletion_df.name, deletion_df.ratio))
        
        raw = parsed_args["search"]["exact"]["raw"]
        locus = parsed_args["search"]["exact"]["locus"]
        ref_fasta = parsed_args["search"]["exact"]["ref-fasta"]
        
        # Build sequence lookup if reference FASTA is provided
        sequence_lookup = nothing
        if ref_fasta !== nothing
            sequence_lookup = build_sequence_lookup(ref_fasta)
        end
        
        counts_df = immunediscover_module.exact.exact_search(table, db, gene, 
                                                             mincount=mincount, minratio=minratio, 
                                                             expect_dict=expect_dict, affix=affix, 
                                                             rss=rss, extension=extension, N=top, 
                                                             raw=raw, sequence_lookup=sequence_lookup, border=border, adjust_per_gene_extension=adjust_per_gene_extension, adjust_percent=adjust_percent)
        sort!(counts_df, [:case, :db_name])
        
        if !parsed_args["search"]["exact"]["noplot"]
            if nrow(counts_df) > 0
                println(immunediscover_module.plotgenes(counts_df))
            else
                @warn "No exact matches to plot"
            end
        end
        
        # Add gene count frequency (aggregation with alleles starting with locus only)
        @info "Excluding genes not starting with $locus for correct allele and gene count frequency calculation"
        
        # Calculate gene counts per well/case/gene
        transform!(groupby(counts_df, [:well, :case, :gene])) do group
            filtered_group = filter(row -> startswith(row.db_name, locus), group)
            gene_count = isempty(filtered_group) ? 0 : sum(filtered_group.count)
            return DataFrame(gene_count = fill(gene_count, nrow(group)))
        end
        
        # Calculate total counts per well/case
        transform!(groupby(counts_df, [:well, :case])) do group
            filtered_group = filter(row -> startswith(row.db_name, locus), group)
            case_count = isempty(filtered_group) ? 0 : sum(filtered_group.count)
            return DataFrame(case_count = fill(case_count, nrow(group)))
        end
        
        # Calculate median count per gene across ALL cases
        transform!(groupby(counts_df, [:gene])) do group
            filtered_group = filter(row -> startswith(row.db_name, locus), group)
            cross_case_median_count = isempty(filtered_group) ? 0 : median(filtered_group.count)
            return DataFrame(cross_case_median_count = fill(cross_case_median_count, nrow(group)))
        end
        
        # Calculate median gene_count for each gene across all cases
        transform!(groupby(counts_df, [:gene])) do group
            filtered_group = filter(row -> startswith(row.db_name, locus), group)
            cross_case_median_gene_count = isempty(filtered_group) ? 0 : median(filtered_group.gene_count)
            return DataFrame(cross_case_median_gene_count = fill(cross_case_median_gene_count, nrow(group)))
        end
        
        # Calculate median allele count for each allele (db_name) across all cases
        transform!(groupby(counts_df, [:db_name])) do group
            filtered_group = filter(row -> startswith(row.db_name, locus), group)
            cross_case_median_allele_count = isempty(filtered_group) ? 0 : median(filtered_group.count)
            return DataFrame(cross_case_median_allele_count = fill(cross_case_median_allele_count, nrow(group)))
        end
        
        # Frequencies within a case
        counts_df[:,:allele_case_freq] = counts_df.count ./ counts_df.case_count
        counts_df[:,:gene_case_freq] = counts_df.gene_count ./ counts_df.case_count
        
        # Ratios compared to cross-case medians
        counts_df[:,:allele_to_cross_case_median_ratio] = counts_df.count ./ counts_df.cross_case_median_allele_count
        counts_df[:,:gene_to_cross_case_median_ratio] = counts_df.gene_count ./ counts_df.cross_case_median_gene_count
        
        # Calculate allele frequency within gene (per case)
        transform!(groupby(counts_df, [:well, :case, :gene]), :count => (x->x./sum(x)) => :allele_freq)
        
        # Apply original filters
        filter!(x -> (x.allele_freq >= immunediscover_module.get_ratio_threshold(expect_dict, x, type="allele_freq")), counts_df)
        filter!(x -> (x.gene_case_freq >= immunediscover_module.get_ratio_threshold(deletion_dict, x, type="gene_case_freq")), counts_df)
        
        # Use the existing min-allele-mratio and min-gene-mratio parameters for cross-case comparisons
        min_allele_cross_case_ratio = parsed_args["search"]["exact"]["min-allele-mratio"]
        min_gene_cross_case_ratio = parsed_args["search"]["exact"]["min-gene-mratio"]
        
        # Filter out alleles/genes that are much lower than expected across all cases
        filter!(x -> (x.allele_to_cross_case_median_ratio >= min_allele_cross_case_ratio), counts_df)
        filter!(x -> (x.gene_to_cross_case_median_ratio >= min_gene_cross_case_ratio), counts_df)
        
        output = always_gz(parsed_args["search"]["exact"]["output"])
        
        # Compute refgene ratios
        if length(refgenes) > 0
            for refgene in refgenes
                counts_df = grouped_ratios(counts_df, refgene, count_col=:count)
                transform!(groupby(counts_df, [:well, :case, :gene]), :count => sum => :ref_gene_count)
                counts_df = grouped_ratios(counts_df, refgene, count_col=:ref_gene_count)
            end
        end
        
        CSV.write(output, counts_df, compress=true, delim='\t')
        @info "Exact search data saved in compressed $output file"

        # No additional stats files are written; totals were reported via @info earlier
    end

    export grouped_ratios, transform_counts, build_sequence_lookup, handle_exact
end
