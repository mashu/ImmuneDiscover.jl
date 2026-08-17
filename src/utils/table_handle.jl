map_some_list(::Absent) = absent
map_some_list(s::Present) = Present(parse_csv_list(s.value))

function route_join(parsed_args, subcmd, how::JoinKind)
    block = parsed_args["table"][subcmd]
    left_keys = parse_csv_list(block["keys"])
    right_keys = or_default(map_some_list(optional(get(block, "right-keys", nothing))), left_keys)
    left_prefix = or_default(optional(get(block, "left-prefix", nothing)), "")
    right_prefix = or_default(optional(get(block, "right-prefix", nothing)), "")
    left_select = or_default(map_some_list(optional(get(block, "left-select", nothing))), String[])
    right_select = or_default(map_some_list(optional(get(block, "right-select", nothing))), String[])

    join_tsv(block["left"], block["right"], block["output"];
             how=how, left_keys=left_keys, right_keys=right_keys,
             left_prefix=left_prefix, right_prefix=right_prefix,
             left_select=left_select, right_select=right_select)
end

function handle_outerjoin(parsed_args, _, _)
    @info "Performing outer join"
    route_join(parsed_args, "outerjoin", OuterJoin())
end

function handle_leftjoin(parsed_args, _, _)
    @info "Performing left join"
    route_join(parsed_args, "leftjoin", LeftJoin())
end

function handle_transform(parsed_args, _, _)
    @info "Transforming TSV file"
    block = parsed_args["table"]["transform"]
    transform_tsv(block["input"], block["output"];
        column=block["column"], pattern=block["pattern"], replacement=block["replacement"],
        new_column=optional(get(block, "new-column", nothing)))
end

function handle_aggregate(parsed_args, _, _)
    @info "Aggregating TSV file"
    block = parsed_args["table"]["aggregate"]
    group_by = parse_csv_list(block["group-by"])
    keep_columns = optional(get(block, "keep-columns", nothing))
    count_column = get(block, "count-column", "count")
    aggregate_tsv(block["input"], block["output"];
        group_by=group_by, keep_columns=keep_columns, count_column=count_column)
end

function handle_unique(parsed_args, _, _)
    columns = parse_csv_list(parsed_args["table"]["unique"]["columns"])
    unique_tsv(parsed_args["table"]["unique"]["input"], parsed_args["table"]["unique"]["output"]; columns=columns)
end

function handle_sort(parsed_args, _, _)
    columns = parse_csv_list(parsed_args["table"]["sort"]["columns"])
    reverse = get(parsed_args["table"]["sort"], "reverse", false)
    sort_tsv(parsed_args["table"]["sort"]["input"], parsed_args["table"]["sort"]["output"]; columns=columns, reverse=reverse)
end

function handle_filter(parsed_args, _, _)
    block = parsed_args["table"]["filter"]
    filter_tsv(block["input"], block["output"], block["column"];
               pattern=optional(get(block, "pattern", nothing)),
               operator=optional(get(block, "operator", nothing)),
               threshold=optional(get(block, "threshold", nothing)))
end

function handle_select(parsed_args, _, _)
    columns = parse_csv_list(parsed_args["table"]["select"]["columns"])
    select_tsv(parsed_args["table"]["select"]["input"], parsed_args["table"]["select"]["output"]; columns=columns)
end

function handle_fasta_export(parsed_args, immunediscover_module, _)
    @info "Exporting TSV to FASTA"
    block = parsed_args["table"]["fasta"]
    immunediscover_module.Fasta.extract_sequences_to_fasta(
        block["input"], block["output"];
        colname = block["colname"],
        colseq = block["colseq"],
        coldesc = optional(block["coldesc"]),
        filter_pattern = optional(block["filter"]),
        desc_filter_pattern = optional(block["desc-filter"]),
        cleanup_pattern = optional(block["cleanup"]),
        sort_by_name = !block["no-sort"],
        mincase = block["mincase"],
        case_col = block["case-col"],
        unique_sequences = block["unique-sequences"]
    )
end

function handle_collect(parsed_args, _, _)
    @info "Collecting TSV files"
    pattern = parsed_args["table"]["collect"]["pattern"]
    output = parsed_args["table"]["collect"]["output"]
    files = glob(pattern)
    @info "Found $(length(files)) files matching pattern $pattern"
    collect_files(files, output)
end

function collect_files(files, output)
    isempty(files) && (@warn "No files found matching pattern"; return)
    collected = DataFrame[]
    first_df = load_and_tag(first(files))
    first_file_columns = names(first_df)
    push!(collected, first_df)
    for file in files[2:end]
        df = load_and_tag(file)
        @assert first_file_columns == names(df) "Column names in $file do not match first file"
        push!(collected, df)
    end
    collected_df = vcat(collected...)
    CSV.write(output, collected_df, delim='\t', compress=true)
    @info "Collected $(nrow(collected_df)) rows into $output"
end

function load_and_tag(file)
    df = load_tsv(file)
    df[:, :file] .= file
    return df
end

function handle_exclude(parsed_args, immunediscover_module, always_gz)
    @info "Exclude"
    db = immunediscover_module.load_fasta(parsed_args["table"]["exclude"]["fasta"])
    colname = parsed_args["table"]["exclude"]["colname"]
    colseq = parsed_args["table"]["exclude"]["colseq"]
    data_df = load_tsv(parsed_args["table"]["exclude"]["input"])
    discard_seqs = Set{String}()
    for row in eachrow(data_df)
        for (name, seq) in db
            if occursin(row[colseq], seq)
                @info "Allele $(row[colname]) is a substring of $(name)"
                push!(discard_seqs, row[colseq])
            end
            if occursin(seq, row[colseq])
                @info "Allele $(name) is a substring of $(row[colname])"
                push!(discard_seqs, row[colseq])
            end
        end
    end
    for seq in discard_seqs
        @info "Discarding sequences matching: $seq"
        filter!(x -> x[colseq] != seq, data_df)
    end
    output = always_gz(parsed_args["table"]["exclude"]["output"])
    CSV.write(output, data_df, compress=true, delim='\t')
end
