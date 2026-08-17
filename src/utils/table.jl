module Table
    using CSV
    using DataFrames
    using Logging
    using Glob
    using ..Option: Absent, Present, absent, optional, or_default

    export outerjoin_tsv, leftjoin_tsv, filter_tsv, transform_tsv, aggregate_tsv, unique_tsv, sort_tsv, select_tsv
    export OuterJoin, LeftJoin, JoinKind

    "Parse a comma-separated CLI value into trimmed column names."
    parse_csv_list(s::AbstractString) = String.(strip.(split(s, ',')))

    function write_tsv(output_file, df)
        output_gz = endswith(output_file, ".gz") ? output_file : output_file * ".gz"
        CSV.write(output_gz, df, compress=true, delim='\t')
        return output_gz
    end

    function require_column(df, column)
        column in names(df) && return nothing
        error("Column '$column' not found in input file. Available columns: $(join(names(df), ", "))")
    end

    function require_columns(df, columns)
        for col in columns
            require_column(df, col)
        end
        return nothing
    end

    load_tsv(path) = CSV.File(path, delim='\t') |> DataFrame

    include("table_join.jl")
    include("table_filter.jl")
    include("table_transform.jl")
    include("table_reshape.jl")
    include("table_handle.jl")
end
