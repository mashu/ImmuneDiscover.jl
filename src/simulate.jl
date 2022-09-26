module simulate
    using FASTX
    using CSV
    using DataFrames

    """
    Generate fastq sequences
    """
    function generate_fastq(prefix::T ,core::T, suffix::T, n::Vector{Int}) where T<:Vector{String}
        @assert length(prefix) == length(core) == length(suffix) == length(n)
        records = []
        nread = 1
        for i in 1:length(n)
            for j in 1:n[i]
                seqlen = length(core[i])+length(prefix[i])+length(suffix[i])
                record = FASTQRecord("read$nread", prefix[i]*core[i]*suffix[i], fill(40, seqlen))
                push!(records, record)
                nread += 1  # Advance one read
            end
        end
        return records
    end

    function generate_indices(prefix::T, suffix::T,case::T) where T<:Vector{String}
        rows = []
        @assert length(case) == length(prefix) == length(suffix)
        for i in 1:length(case)
                push!(rows,(prefix[i],suffix[i],case[i]))
        end
        df = DataFrame(rows)
        rename!(df,[:forward_index,:reverse_index,:case])
        return df
    end
end