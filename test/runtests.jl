using immunediscover
using Test
using CSV
using DataFrames

#
# Run testsets
#
@testset verbose = true "immunediscover" begin
    @testset "simulate.jl" begin
        records = immunediscover.simulate.generate_fastq(["AAA", "AAT"], ["CACAGTG", "CACAGTG"], ["TTT", "TTT"], [10, 5])
        indices = immunediscover.simulate.generate_indices(["AAA", "AAT"], ["TTT", "TTT"], ["A1","A2"])
        immunediscover.data.write_fastq("test.fastq", records)
        CSV.write("test_indices.tsv", indices, delim='\t')
        @test length(records) == 15
        @test nrow(indices) == 2
    end

    @testset "demultiplex.jl" begin
        n = 0
        immunediscover.demultiplex.process_fastq("test.fastq") do r
            n += 1
        end
        @test n == 15
    end

    @testset "cli.jl" begin
        empty!(ARGS)
        append!(ARGS, ["demultiplex", "test.fastq", "test_indices.tsv", "test.tsv"])
        parsed_args = immunediscover.cli.parse_commandline()
        @test parsed_args["%COMMAND%"] == "demultiplex"
        @test parsed_args["demultiplex"]["output"] == "test.tsv"
        @test parsed_args["demultiplex"]["fastq"] == "test.fastq"
        @test parsed_args["demultiplex"]["indices"] == "test_indices.tsv"    
    end
    @testset "immunediscover.jl" begin
        @test immunediscover.julia_main() == 0
        @inferred immunediscover.real_main()
    end
end
