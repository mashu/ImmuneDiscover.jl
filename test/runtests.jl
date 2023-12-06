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
include("../src/hamming.jl")
include("../src/patterns.jl")
include("../src/heptamer.jl")

using .cli
using .demultiplex
using .simulate
using .data
using .profile
using .trim
using .exact
using .hamming
using .patterns
using .heptamer

# Initialize a dictionary to track test outcomes
test_outcomes = Dict(
    "simulate" => false,
    "demultiplex" => false
)

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
        CSV.write("test_indices.tsv", indices, delim='\t')
        db = data.load_fasta("test.fasta", validate=false)
        @test length(db) == 15
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
        test_outcomes["simulate"] = true
    end

    if test_outcomes["simulate"]
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
            table, stats = demultiplex.demux("test.fastq", "test_indices.tsv")
            CSV.write("test.tsv.gz", table, delim='\t')
        end
        test_outcomes["demultiplex"] = true
    else
        println("Skipping remaining tests due to failed dependencies")
    end

    if test_outcomes["demultiplex"]
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
            db = data.load_fasta("test.fasta", validate=false)
            mincount = 1
            minfreq = 0.01
            counts_df = exact.exact_search(table, db, mincount=mincount, minfreq=minfreq)
            sort!(counts_df, [:case, :db_name])
            @test nrow(counts_df) == 30
            @test counts_df[1, :count] == 10
        end

        @testset "hamming.jl" begin
            # CLI
            empty!(ARGS)
            append!(ARGS, ["hamming", "test.tsv.gz", "test.fasta", "test_hamming.tsv.gz", "--maxdist", "10", "-f", "genomic_sequence"])
            parsed_args = cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "hamming"
            @test parsed_args["hamming"]["tsv"] == "test.tsv.gz"
            @test parsed_args["hamming"]["fasta"] == "test.fasta"
            @test parsed_args["hamming"]["output"] == "test_hamming.tsv.gz"
            @test parsed_args["hamming"]["maxdist"] == 10
            @test parsed_args["hamming"]["column"] == "genomic_sequence"
            # Module
            table = CSV.File("test.tsv.gz", delim='\t') |> DataFrame
            db = data.load_fasta("test.fasta", validate=false)
            mincount = 1
            minfreq = 0.01
            counts = hamming.hamming_search(table, db, max_dist=10, column="genomic_sequence", check_bounds=true, umi=false)
            summary = hamming.summarize(counts, db, min_count=mincount, cluster_ratio=minfreq)
            @test nrow(counts) == 15
            @test nrow(summary) == 2
            @test summary[1, :count] == 10
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
            # Argument Parsing Tests
            @testset "parse_commandline" begin
                empty!(ARGS)
                append!(ARGS, ["trim", "test.tsv.gz", "test.fasta", "test_trim.tsv", "-l", "7"])
                parsed_args = cli.parse_commandline(ARGS)

                @test parsed_args["%COMMAND%"] == "trim"
                @test parsed_args["trim"]["input"] == "test.tsv.gz"
                @test parsed_args["trim"]["fasta"] == "test.fasta"
                @test parsed_args["trim"]["output"] == "test_trim.tsv"
                @test parsed_args["trim"]["length"] == 7
            end

            # Preparing test data
            read = "AAAAAAAAAACACAGTGATGTAGCACAGTGTTTTTTTTTT"
            motifs = ["CACAGTG", "CACAGTG", "CACATGT"]
            # Generating profiles for motifs
            prof_start = profile.motif_profile(motifs)
            prof_stop = prof_start

            # Find Max Probability Position Tests
            @testset "find_maxprob_pos" begin
                # Finding maximum probability positions for trimming
                start = trim.find_maxprob_pos(read, prof_start) + size(prof_start, 2)
                stop = trim.find_maxprob_pos(read[start:end], prof_stop) - 1

                # Applying trimming to the read
                start_trim = read[start:end]
                stop_trim = start_trim[1:stop]

                # Testing the trimming procedure
                @test start > 0
                @test stop > 0
                @test length(stop_trim) == 6
            end
            # Trim Profile Tests
            let prof
                @testset "trim_profiles" begin
                    # Test trimming process
                    prof = trim.trim_profiles(motifs, 3)
                    # Test profile positions
                    @test size(prof[1]) == (4, 3)
                    @test size(prof[2]) == (4, 3)
                end

                # Trim Read Tests
                @testset "trim_reads" begin
                    tstart, tstop = trim.trim_read(read, prof[1], prof[2])
                end
            end
        end

        @testset "pattern.jl" begin
            # Argument Parsing Tests
            @testset "parse_commandline" begin
                empty!(ARGS)
                append!(ARGS, ["pattern", "test.tsv.gz", "test.fasta", "test_pattern.tsv.gz"])
                parsed_args = cli.parse_commandline(ARGS)
                @test parsed_args["%COMMAND%"] == "pattern"
                @test parsed_args["pattern"]["input"] == "test.tsv.gz"
                @test parsed_args["pattern"]["fasta"] == "test.fasta"
                @test parsed_args["pattern"]["output"] == "test_pattern.tsv.gz"
            end

            # Preparing test data
            table = CSV.File("test.tsv.gz", delim='\t') |> DataFrame
            db = data.load_fasta("test.fasta", validate=false)
            db_df = DataFrame(db, [:id, :seq])
            db_df[:,:gene] = [split(x, "*")[1] for x in db_df[:,:id]]
            blacklist = DataFrame(id=["test1*01", "test2*02"], seq=["AAAAAAAA", "CCCCCCCC"], gene=["test1","test2"])

            # Search Patterns Tests
            @testset "search_patterns" begin
                # Testing the search_patterns function
                output = patterns.search_patterns(table, blacklist, db_df, min_count=1, min_freq=0.01, max_fragments=10, fragment_size=7, max_fragment_size=7)
                # Add specific assertions related to the output, if applicable
            end

            # Tests for split_kmers function
            @testset "split_kmers" begin
                @test split_kmers("ATCGATCG", 3) == ["CGA", "GAT"]
                @test isempty(split_kmers("", 3))
                @test isempty(split_kmers("ATCG", 5))
            end

            # Tests for find_unique_fragments function
            @testset "find_unique_fragments" begin
                group_df = DataFrame(seq = ["ATCGATCG", "ATCGATCGAAA", "ATCGATCGCCCC"])
                outgroup_df = DataFrame(seq = ["AATTCCGG", "TTAACCG"])

                @test find_unique_fragments(group_df, outgroup_df; fragment_size = 4, max_fragment_size = 4) == ["GATC", "CGAT"]
                @test isempty(find_unique_fragments(DataFrame(seq = String[]), DataFrame(seq = String[])))
            end
        end

        @testset "heptamer.jl" begin
            # Argument Parsing Tests
            @testset "parse_commandline" begin
                empty!(ARGS)
                append!(ARGS, ["heptamer", "test.tsv.gz", "test.fasta", "test_heptamer.tsv.gz", "test_summary.tsv"])
                parsed_args = cli.parse_commandline(ARGS)

                @test parsed_args["%COMMAND%"] == "heptamer"
                @test parsed_args["heptamer"]["tsv"] == "test.tsv.gz"
                @test parsed_args["heptamer"]["fasta"] == "test.fasta"
                @test parsed_args["heptamer"]["output"] == "test_heptamer.tsv.gz"
                @test parsed_args["heptamer"]["summary"] == "test_summary.tsv"
                # Assuming "chain" and other keys have default values if not provided in ARGS
                @test parsed_args["heptamer"]["chain"] == "IGHV"
            end

            # Preparing test data
            table = CSV.File("test.tsv.gz", delim='\t') |> DataFrame
            db = data.load_fasta("test.fasta", validate=false)
            parsed_args = cli.parse_commandline(["heptamer", "test.tsv.gz", "test.fasta", "test_heptamer.tsv.gz", "test_summary.tsv"])
            heptamers = heptamer.load_heptamers(parsed_args["heptamer"]["json"])

            # Extract Heptamers Tests
            let heptamer_df
                @testset "extract_heptamers" begin
                    # Extract heptamers and summarize data
                    heptamer_df = heptamer.extract_heptamers(table, db, heptamers[parsed_args["heptamer"]["chain"]]; max_dist=parsed_args["heptamer"]["maxdist"], b=parsed_args["heptamer"]["begin"]+1, e=parsed_args["heptamer"]["end"])

                    @test nrow(heptamer_df) == 225
                end

                # Summarize Heptamers Tests
                @testset "extract_heptamers" begin
                    summary_df = heptamer.summarize(heptamer_df, db, ratio=parsed_args["heptamer"]["ratio"], count=parsed_args["heptamer"]["mincount"])

                    @test nrow(summary_df) == 30
                    @test summary_df[1, :allele_count] == 10
                end
            end
        end
    else
        println("Skipping remaining tests due to failed dependencies")
    end
end

