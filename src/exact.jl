module exact
    using CSV
    using DataFrames
    using ProgressMeter
    using Folds

    """
        extract_flanking(genomic_sequence::String, range::Tuple{Int, Int}, gene_type::String, n::Int)

    Extract flanking sequences from genomic sequence
    """
    function extract_flanking(genomic_sequence::String, range::Tuple{Int, Int}, gene_type::String, n::Int)
        # Ensure range is valid
        start_pos, end_pos = range
        if start_pos < 1 || end_pos > length(genomic_sequence) || start_pos > end_pos
            error("Invalid range")
        end

        sequence = genomic_sequence[start_pos:end_pos]
        if gene_type == "V"
            prefix = start_pos > n ? genomic_sequence[(start_pos - n):(start_pos - 1)] : genomic_sequence[1:(start_pos - 1)]
            heptamer = genomic_sequence[end_pos + 1:min(end_pos + 7, length(genomic_sequence))]
            spacer = genomic_sequence[end_pos + 8:min(end_pos + 30, length(genomic_sequence))]
            nonamer = genomic_sequence[end_pos + 31:min(end_pos + 39, length(genomic_sequence))]
            return (prefix=prefix, sequence=sequence, heptamer=heptamer, spacer=spacer, nonamer=nonamer)
        elseif gene_type == "J"
            nonamer = genomic_sequence[max(1, start_pos - 9):(start_pos - 1)]
            spacer = genomic_sequence[max(1, start_pos - 32):max(1, start_pos - 10)]
            heptamer = genomic_sequence[max(1, start_pos - 40):max(1, start_pos - 34)]
            suffix = genomic_sequence[end_pos + 1:min(end_pos + n, length(genomic_sequence))]
            return (nonamer=nonamer, spacer=spacer, heptamer=heptamer, suffix=suffix, sequence=sequence)
        elseif gene_type == "D"
            pre_nonamer = genomic_sequence[max(1, start_pos - 9):(start_pos - 1)]
            pre_spacer = genomic_sequence[max(1, start_pos - 21):max(1, start_pos - 10)]
            pre_heptamer = genomic_sequence[max(1, start_pos - 28):max(1, start_pos - 22)]
            post_heptamer = genomic_sequence[end_pos + 1:min(end_pos + 7, length(genomic_sequence))]
            post_spacer = genomic_sequence[end_pos + 8:min(end_pos + 19, length(genomic_sequence))]
            post_nonamer = genomic_sequence[end_pos + 20:min(end_pos + 28, length(genomic_sequence))]
            return (pre_nonamer=pre_nonamer, pre_spacer=pre_spacer, pre_heptamer=pre_heptamer,sequence=sequence, post_heptamer=post_heptamer, post_spacer=post_spacer, post_nonamer=post_nonamer)
        else
            error("Invalid gene type")
        end
    end

    """
        exact_search(table, query, gene; mincount=10, minfreq=0.2, collapse=true)

    Exact search per case.
    Parameters:
        table::DataFrame: Table containing genomic sequences
        query::Dict{String, String}: Dictionary containing database name and sequence
        gene::String: Gene type (V, D, J)
        mincount::Int: Minimum number of reads for a sequence to be considered
        minfreq::Float: Minimum frequency of reads for a sequence to be considered
        full_mincount::Int: Minimum number of reads for a full record to be considered
        full_minfreq::Float: Minimum frequency of reads for a full record to be considered
        collapse::Bool: If true, only return highest count full record, otherwise return all records
    """
    function exact_search(table, query, gene; mincount=10, minfreq=0.01, full_mincount=2, full_minfreq=0.01, collapse=true)
        @info "Using mincount: $mincount and minfreq: $minfreq for genes: $gene"
        @assert all(names(table) .== ["well","case","name","genomic_sequence"]) "File must contain following columns case, name, genomic_sequence"
        p = Progress(nrow(table))
        result = Folds.map(eachrow(table)) do row
            next!(p)
            case = row.case
            well = row.well
            matches = Vector{NamedTuple}()
            @inbounds for (name, seq) in query
                match = findfirst(seq, row.genomic_sequence)
                if match !== nothing
                    flanks = extract_flanking(row.genomic_sequence, (minimum(match), maximum(match)), gene, 10)
                    meta = (well=string(well), case=string(case), db_name=string(name))
                    push!(matches, merge(meta, flanks))
                end
            end
            matches
        end
        result_df = DataFrame(reduce(vcat,[r for r in result if r != []]))

        # Apply filters
        df = transform(groupby(result_df, names(result_df)), nrow => :full_count)
        transform!(groupby(df, [:well, :case, :db_name, :sequence]), nrow => :count)
        transform!(df, :db_name => ByRow(x -> first(split(x, '*'))) => :gene)
        transform!(groupby(df, [:well, :case, :gene]), :full_count => (x->x./maximum(x))=> :full_frequency)
        transform!(groupby(df, [:well, :case, :gene]), :count => (x->x./maximum(x))=> :frequency)
        sort!(df, [:full_count, :count], rev=[true,true])
        filter!(r->(r.full_count >= full_mincount) & (r.full_frequency >= full_minfreq), df)  # Order is important here to filter first on smaller numbers
        filter!(r->(r.count >= mincount) & (r.frequency >= minfreq) , df)
        udf = unique(df)

        # Reorder columns
        priority_columns = ["well", "case", "gene", "db_name", "count", "full_count", "frequency", "full_frequency"]
        remaining_columns = setdiff(names(udf), priority_columns)
        udf[:, vcat(priority_columns, remaining_columns)]
        # Return only first highest count full record, otherwise return all records
        if collapse
            return unique(udf, [:well, :case, :gene, :db_name, :sequence])
        end
        return udf
    end
end