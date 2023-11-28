using immunediscover
using Test
using CSV
using DataFrames

include("../src/cli.jl")
include("../src/demultiplex.jl")
include("../src/simulate.jl")
include("../src/data.jl")
include("../src/profile.jl")
include("../src/trim.jl")
include("../src/exact.jl")

using .cli
using .demultiplex
using .simulate
using .data
using .profile
using .trim
using .exact

#
# Run testsets
#
@testset verbose = true "immunediscover" begin
    @testset "simulate.jl" begin
        fastq_records = simulate.generate_fastq(["A"^10, "A"^9*"T"], ["CACAGTGATGTAGCACAGTG", "CACAGTGATGTAGCACAGTG"], ["T"^10, "T"^10], [10, 5])
        fasta_records = simulate.generate_fasta(["CACAGTG", "CACAGTG"], ["ATGTAG", "ATGTAG"], ["CACAGTG", "CACAGTG"], [10, 5])
        indices = simulate.generate_indices(["A"^10, "A"^9*"T"], ["T"^10, "T"^10], ["A1", "A2"])
        data.write_fastq("test.fastq", fastq_records)
        data.write_fasta("test.fasta", fasta_records)
        data.write_gz_fastq("test.fastq.gz", fastq_records)
        data.write_gz_fasta("test.fasta.gz", fasta_records)
        db = data.load_fasta("test.fasta")
        @test length(db) == 15
        CSV.write("test_indices.tsv", indices, delim='\t')
        @test length(fastq_records) == 15
        @test length(fasta_records) == 15
        @test nrow(indices) == 2
        n = 0
        for file in ["test.fastq", "test.fastq.gz"]
            data.process_fastq(file) do r
                n += 1
            end
        end
        @test n == 30
    end

    @testset "demultiplex.jl" begin
        # CLI
        @test endswith(cli.always_gz("test.fasta.gz"), ".gz")
        @test endswith(cli.always_gz("test.fasta"), ".gz")
        empty!(ARGS)
        append!(ARGS, ["demultiplex", "test.fastq", "test_indices.tsv", "test.tsv"])
        parsed_args = cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "demultiplex"
        @test parsed_args["demultiplex"]["fastq"] == "test.fastq"
        @test parsed_args["demultiplex"]["indices"] == "test_indices.tsv"
        @test parsed_args["demultiplex"]["output"] == "test.tsv"
    end

    @testset "exact.jl" begin
        # CLI
        empty!(ARGS)
        append!(ARGS, ["exact", "test.tsv.gz", "test.fasta", "test_exact.tsv.gz"])
        parsed_args = cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "exact"
        @test parsed_args["exact"]["tsv"] == "test.tsv.gz"
        @test parsed_args["exact"]["fasta"] == "test.fasta"
        @test parsed_args["exact"]["output"] == "test_exact.tsv.gz"
        # Module
        table = CSV.File("test.tsv.gz", delim='\t') |> DataFrame
        db = data.load_fasta("test.fasta")
        mincount = 1
        minfreq = 0.01
        counts_df = exact.exact_search(table, db, mincount=mincount, minfreq=minfreq)
        sort!(counts_df, [:case, :db_name])
        @test nrow(counts_df) == 30
        @test counts_df[1, :count] == 10
    end

    motifs = ["CACAGTG", "CACAGTG", "CACATGT"]
    @testset "profile.jl" begin
        countsum = sum(profile.counts(motifs))
        @test countsum == 21
        profsum = sum(profile.motif_profile(motifs))
        @test profsum < 10
        @test profsum > 1
        prof = profile.motif_profile(motifs)
        @test size(prof) == (4, 7)
        prob = profile.motif_prob("CACATTT", prof)
        @test prob > 0.07
        prob = profile.motif_prob("-ACATTT", prof)
        @test prob == 0
    end

    @testset "trim.jl" begin
        # CLI
        empty!(ARGS)
        append!(ARGS, ["trim", "test.tsv.gz", "test.fasta", "test_trim.tsv","-l","7"])
        parsed_args = cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "trim"
        @test parsed_args["trim"]["input"] == "test.tsv.gz"
        @test parsed_args["trim"]["fasta"] == "test.fasta"
        @test parsed_args["trim"]["output"] == "test_trim.tsv"
        # Module
        read = "AAAAAAAAAACACAGTGATGTAGCACAGTGTTTTTTTTTT"
        motifs = ["CACAGTG", "CACAGTG", "CACATGT"]
        prof_start = profile.motif_profile(motifs)
        prof_stop = prof_start
        start = trim.find_maxprob_pos(read, prof_start) + size(prof_start, 2)
        start_trim = read[start:end]
        stop = trim.find_maxprob_pos(start_trim, prof_stop)-1
        stop_trim = start_trim[1:stop]
        p_start, p_stop = trim.trim_profiles(motifs, 3)
        tstart, tstop = trim.trim_read(read, p_start, p_stop)
        @test tstart > 0
        @test tstop > 0
        @test length(stop_trim) == 6
    end
end

