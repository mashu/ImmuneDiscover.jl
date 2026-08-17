abstract type JoinKind end
struct OuterJoin <: JoinKind end
struct LeftJoin <: JoinKind end

join_apply(::OuterJoin, left, right; kwargs...) = outerjoin(left, right; kwargs...)
join_apply(::LeftJoin, left, right; kwargs...) = leftjoin(left, right; kwargs...)

prefix_columns!(df, keys, prefix::AbstractString) = prefix_columns!(df, keys, optional(prefix))
prefix_columns!(df, _, ::Absent) = df
function prefix_columns!(df, keys, p::Present)
    for col in names(df)
        col in keys && continue
        rename!(df, col => "$(p.value)_$(col)")
    end
    return df
end

"""
    join_tsv(left_file, right_file, output_file; how, left_keys, right_keys, left_prefix, right_prefix, left_select, right_select)

Perform a join (outer or left) on two TSV files by specified key columns.
"""
function join_tsv(left_file, right_file, output_file;
                  how::JoinKind,
                  left_keys,
                  right_keys=left_keys,
                  left_prefix="",
                  right_prefix="",
                  left_select=String[],
                  right_select=String[])

    @info "Loading left file: $left_file"
    left_df = load_tsv(left_file)
    @info "Loading right file: $right_file"
    right_df = load_tsv(right_file)

    left_df = select_join_columns(left_df, left_select, "left")
    right_df = select_join_columns(right_df, right_select, "right")

    @info "Left file: $(nrow(left_df)) rows, $(ncol(left_df)) columns"
    @info "Right file: $(nrow(right_df)) rows, $(ncol(right_df)) columns"
    @info "Join keys - Left: $left_keys, Right: $right_keys"

    left_df_prefixed = prefix_columns!(copy(left_df), left_keys, left_prefix)
    right_df_prefixed = prefix_columns!(copy(right_df), right_keys, right_prefix)

    @info "Performing join"
    result_df = join_apply(how, left_df_prefixed, right_df_prefixed, on=left_keys .=> right_keys, makeunique=true)
    @info "Result: $(nrow(result_df)) rows, $(ncol(result_df)) columns"
    output_gz = write_tsv(output_file, result_df)
    @info "Output saved to: $output_gz"
    return result_df
end

select_join_columns(df, cols::AbstractVector{<:AbstractString}, side) =
    isempty(cols) ? df : (@info "Selecting columns from $side file: $cols"; df[:, cols])
select_join_columns(df, ::Absent, _) = df
select_join_columns(df, s::Present{<:AbstractVector}, side) = select_join_columns(df, s.value, side)

function outerjoin_tsv(left_file, right_file, output_file; kwargs...)
    join_tsv(left_file, right_file, output_file; how=OuterJoin(), kwargs...)
end

function leftjoin_tsv(left_file, right_file, output_file; kwargs...)
    join_tsv(left_file, right_file, output_file; how=LeftJoin(), kwargs...)
end
