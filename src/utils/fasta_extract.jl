"""
    extract_sequences_to_fasta(input_file, output_file; kwargs...)

Extract sequences from a TSV file and save them to a FASTA file.
Optional description / regex arguments are `absent` or `Present` (strings are converted).
"""
function extract_sequences_to_fasta(input_file::String, output_file::String;
                                   colname::String="allele_name",
                                   colseq::String="seq",
                                   coldesc=absent,
                                   filter_pattern=absent,
                                   desc_filter_pattern=absent,
                                   cleanup_pattern=absent,
                                   sort_by_name::Bool=true,
                                   mincase::Int=1,
                                   case_col::String="case",
                                   unique_sequences::Bool=false)
    @info "Extracting sequences to FASTA format"
    df = CSV.File(input_file, delim='\t') |> DataFrame
    colseq_list = strip.(split(colseq, ','))
    @assert colname ∈ names(df) "Input file must have $colname column"
    for seq_col in colseq_list
        @assert seq_col ∈ names(df) "Input file must have $seq_col column"
    end
    length(colseq_list) > 1 && @info "Using multiple sequence columns for concatenation: $(join(colseq_list, ", "))"

    desc = optional(coldesc)
    name_filter = optional(filter_pattern)
    desc_filter = optional(desc_filter_pattern)
    cleanup = optional(cleanup_pattern)
    require_desc_column!(df, desc)
    @info "Loaded $(nrow(df)) rows from input file"

    df = apply_name_filter(df, colname, name_filter)
    df = apply_desc_filter(df, desc, desc_filter)
    df = apply_mincase(df, colname, colseq_list, case_col, mincase)
    df = drop_name_seq_duplicates(df, colname, colseq_list)
    df = drop_duplicate_sequences(df, colseq_list, unique_sequences)
    sort_by_name && (@info "Sorting alleles by name"; df = sort(df, colname))

    @info "Writing $(nrow(df)) sequences to $output_file"
    write_fasta_records(output_file, df, colname, colseq_list, desc, name_filter, desc_filter, cleanup)
    @info "Successfully extracted $(nrow(df)) unique sequences to $output_file"
    return nothing
end

require_desc_column!(_, ::Absent) = nothing
function require_desc_column!(df, d::Present)
    @assert d.value ∈ names(df) "Input file must have $(d.value) column when coldesc is specified"
    @info "Using $(d.value) column for FASTA descriptions"
    return nothing
end

apply_name_filter(df, _, ::Absent) = df
function apply_name_filter(df, colname, p::Present)
    @info "Filtering rows where $colname matches regex pattern '$(p.value)'"
    df = df[occursin.(Regex(p.value), df[!, colname]), :]
    @info "After filtering: $(nrow(df)) rows remaining"
    return df
end

apply_desc_filter(df, ::Absent, ::Absent) = df
apply_desc_filter(_, ::Absent, ::Present) =
    error("Description filter pattern specified but no coldesc column provided")
apply_desc_filter(df, ::Present, ::Absent) = df
function apply_desc_filter(df, d::Present, p::Present)
    @info "Filtering rows where $(d.value) matches regex pattern '$(p.value)'"
    df = df[occursin.(Regex(p.value), df[!, d.value]), :]
    @info "After description filtering: $(nrow(df)) rows remaining"
    return df
end

function apply_mincase(df, colname, colseq_list, case_col, mincase)
    mincase <= 1 && return df
    @assert case_col ∈ names(df) "Input file must have $case_col column when mincase > 1"
    @info "Filtering alleles present in at least $mincase donors"
    grouping_cols = [colname; colseq_list]
    allele_donor_counts = combine(groupby(df, grouping_cols),
                                case_col => (x -> length(unique(x))) => :donor_count)
    qualifying_alleles = allele_donor_counts[allele_donor_counts.donor_count .>= mincase, grouping_cols]
    df = innerjoin(df, qualifying_alleles, on=grouping_cols)
    @info "After mincase filtering (≥$mincase donors): $(nrow(df)) rows remaining"
    return df
end

function drop_name_seq_duplicates(df, colname, colseq_list)
    before_unique = nrow(df)
    unique_cols = [colname; colseq_list]
    df = DataFrames.unique(df, unique_cols)
    after_unique = nrow(df)
    before_unique != after_unique &&
        @info "Removed $(before_unique - after_unique) duplicate records (by $colname and sequence columns)"
    return df
end

drop_duplicate_sequences(df, _, ::Val{false}) = df
drop_duplicate_sequences(df, colseq_list, keep_unique::Bool) =
    drop_duplicate_sequences(df, colseq_list, Val(keep_unique))
function drop_duplicate_sequences(df, colseq_list, ::Val{true})
    @info "Extracting unique sequences (ignoring sequence names)"
    before_seq_unique = nrow(df)
    df = DataFrames.unique(df, colseq_list)
    after_seq_unique = nrow(df)
    before_seq_unique != after_seq_unique &&
        @info "Removed $(before_seq_unique - after_seq_unique) duplicate sequences with different names"
    return df
end

compiled(::Absent) = absent
compiled(s::Present) = Present(Regex(s.value))

function write_fasta_records(output_file, df, colname, colseq_list, desc, name_filter, desc_filter, cleanup)
    cleanup_rx = compiled(cleanup)
    filter_rx = compiled(name_filter)
    desc_rx = compiled(desc_filter)
    open(output_file, "w") do io
        for row in eachrow(df)
            name = string(row[colname])
            seq = join([string(row[col]) for col in colseq_list], "")
            name = clean_name(name, cleanup_rx)
            name = strip_filter_match(name, filter_rx)
            name = append_description(name, row, desc, desc_rx)
            println(io, ">$name")
            println(io, seq)
        end
    end
end

clean_name(name, ::Absent) = name
function clean_name(name, rx::Present{Regex})
    occursin(rx.value, name) || return name
    return strip(replace(name, rx.value => ""))
end

strip_filter_match(name, ::Absent) = name
function strip_filter_match(name, rx::Present{Regex})
    occursin(rx.value, name) || return name
    return rstrip(replace(name, rx.value => ""))
end

append_description(name, _, ::Absent, _) = name
function append_description(name, row, col::Present, desc_rx)
    desc = string(row[col.value])
    isempty(desc) && return name
    return merge_description(name, desc, desc_rx)
end

merge_description(name, desc, ::Absent) = "$name $desc"
function merge_description(name, desc, rx::Present{Regex})
    match_result = match(rx.value, desc)
    match_result === nothing && return name
    return merge_capture(name, desc, match_result.captures)
end

merge_capture(name, desc, captures) = merge_capture(name, desc, captures, !isempty(captures) && captures[1] !== nothing)
merge_capture(name, _, captures, ::Val{true}) = "$name $(captures[1])"
merge_capture(name, _, _, ::Val{false}) = name
merge_capture(name, desc, captures, keep::Bool) = merge_capture(name, desc, captures, Val(keep))
