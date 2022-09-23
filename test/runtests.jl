using immunediscover
using Test
using FASTX

"""
Generate fastq for testing purpose
"""
function generate_fastq(prefix::T ,core::T, suffix::T, n::Vector{Int}; path="test.fastq") where T<:Vector{String}
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
    FASTQ.Writer(open(path, "w")) do writer
        for record in records
            write(writer, record)
        end
    end
end

generate_fastq(["AAA", "AAT"], ["CACAGTG", "CACAGTG"], ["TTT", "TTT"], [10, 5])

#
# Run testsets
#
@testset verbose = true "immunediscover" begin
    @testset "demultiplex.jl" begin
        n = 0
        immunediscover.demultiplex.process_fastq("test.fastq") do r
            n += 1
        end
        @test n == 15
    end

    @testset "cli.jl" begin
        empty!(ARGS)
        append!(ARGS, ["demultiplex", "test.fastq", "indices.tsv", "test.tsv"])
        parsed_args = immunediscover.cli.parse_commandline()
        @test parsed_args["%COMMAND%"] == "demultiplex"
        @test parsed_args["demultiplex"]["output"] == "test.tsv"
        @test parsed_args["demultiplex"]["fastq"] == "test.fastq"
        @test parsed_args["demultiplex"]["indices"] == "indices.tsv"    
    end
end