module simulate
    using FASTX
    using CSV
    using DataFrames

    """
        generate_fastq(prefix, core, suffix, n)

    Generate fastq sequences
    """
    function generate_fastq(prefix::T ,core::T, suffix::T, n::Vector{Int}) where T<:Vector{String}
        @assert length(prefix) == length(core) == length(suffix) == length(n)
        records = Vector{FASTQRecord}()
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

    """
        generate_fasta(prefix, core, suffix, n)

    Generate fasta sequences
    """
    function generate_fasta(prefix::T ,core::T, suffix::T, n::Vector{Int}) where T<:Vector{String}
        @assert length(prefix) == length(core) == length(suffix) == length(n)
        records = Vector{FASTARecord}()
        nread = 1
        for i in 1:length(n)
            for j in 1:n[i]
                seqlen = length(core[i])+length(prefix[i])+length(suffix[i])
                record = FASTARecord("read$nread", prefix[i]*core[i]*suffix[i])
                push!(records, record)
                nread += 1  # Advance one read
            end
        end
        return records
    end

    """
        generate_indices(prefix, suffix, case)

    Generate indices dataframe
    """
    function generate_indices(prefix::T, suffix::T,case::T) where T<:Vector{String}
        rows = []
        @assert length(case) == length(prefix) == length(suffix)
        for i in eachindex(case)
            push!(rows,(prefix[i],suffix[i],case[i]))
        end
        df = DataFrame(rows)
        rename!(df,[:forward_index,:reverse_index,:case])
        return df
    end
end