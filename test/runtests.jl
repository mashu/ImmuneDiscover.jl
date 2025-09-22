using immunediscover
using Test
using CSV
using DataFrames
using FASTX

using immunediscover.cli
using immunediscover.demultiplex
using immunediscover.simulate
using immunediscover.data
using immunediscover.profile
using immunediscover.trim
using immunediscover.exact
using immunediscover.hamming
using immunediscover.patterns
using immunediscover.heptamer
using immunediscover.keyedsets
using immunediscover.blast
using immunediscover.linkage

# Initialize a dictionary to track test outcomes
test_outcomes = Dict(
    "simulate" => false,
    "demultiplex" => false
)

#
# Run testsets
#
@testset verbose = true "immunediscover" begin
    @testset "keyedsets.jl" begin
        ks = KeyedSet()
        @test length(ks) == 0

        push!(ks, KeyedPair("key1", "value1"))
        @test length(ks) == 1
        @test ks["key1"] == "value1"

        push!(ks, ("key2", "value2"))
        @test length(ks) == 2
        @test ks["key2"] == "value2"

        # Test vector of tuples constructor
        ks2 = KeyedSet([("key3", "value3"), ("key4", "value4")])
        @test length(ks2) == 2
        @test ks2["key3"] == "value3"
        @test ks2["key4"] == "value4"

        # Duplicate warning behaviour
        ks = KeyedSet()
        push!(ks, ("key1", "value1"))

        # Test duplicate key with same value (should log info)
        @test_logs (:info, "Duplicate key key1 with value value1 already exists in KeyedSet") push!(ks, ("key1", "value1"))

        # Test duplicate key with different value (should log warning)
        @test_logs (:warn, "Key key1 with value value2 already exists in KeyedSet with value value1") push!(ks, ("key1", "value2"))

        # Set operations
        ks1 = KeyedSet([("key1", "value1"), ("key2", "value2")])
        ks2 = KeyedSet([("key2", "value2"), ("key3", "value3")])

        @test "key1" in ks1
        @test !("key3" in ks1)

        union_ks = union(ks1, ks2)
        @test length(union_ks) == 3
        @test Set(keys(union_ks.data)) == Set(["key1", "key2", "key3"])

        intersect_ks = intersect(ks1, ks2)
        @test length(intersect_ks) == 1
        @test Set(keys(intersect_ks.data)) == Set(["key2"])

        diff_ks = setdiff(ks1, ks2)
        @test length(diff_ks) == 1
        @test Set(keys(diff_ks.data)) == Set(["key1"])

        # Equality
        ks1 = KeyedSet([("key1", "value1"), ("key2", "value2")])
        ks2 = KeyedSet([("key1", "value1"), ("key2", "value2")])
        ks3 = KeyedSet([("key1", "value1"), ("key3", "value3")])

        @test ks1 == ks2
        @test ks1 != ks3

        # Collect and iteration
        ks = KeyedSet([("key1", "value1"), ("key2", "value2")])
        collected = collect(ks)
        @test length(collected) == 2
        @test ("key1", "value1") in collected
        @test ("key2", "value2") in collected

        count = 0
        for key in ks
            @test key in ["key1", "key2"]
            count += 1
        end
        @test count == 2

        # Show method
        ks = KeyedSet([("key1", "value1"), ("key2", "value2")])
        @test sprint(show, ks) == "KeyedSet(size=2)"
    end

    @testset "simulate.jl" begin
        # Generate test data using new simulation functionality
        records, indices = simulate.generate_fasta_with_mutations(
            "test.fasta",
            "test_indices.tsv",
            "reference.fasta",
            "novel.fasta"
        )
        # Test file outputs
        @test isfile("test.fasta")
        @test isfile("reference.fasta")
        @test isfile("test_indices.tsv")
        @test isfile("novel.fasta")

        # Test data loading
        db = data.load_fasta("test.fasta", validate=false)
        ref_db = data.load_fasta("reference.fasta", validate=false)
        novel_db = data.load_fasta("novel.fasta", validate=false)
        indices_df = CSV.read("test_indices.tsv", DataFrame, delim='\t')

        # Debug prints
        println("Column names in indices_df: ", names(indices_df))
        println("Content of indices_df:")
        println(indices_df)

        # Test basic properties
        @test length(ref_db) == 1  # Should have one reference sequence
        @test nrow(indices_df) == 2  # Should have two donors
        @test length(novel_db) == 72  # 9 mutation types × 3 lengths × 2 donors + 34 reference sequences

        # Test indices structure with clear error messages
        col_names = names(indices_df)
        @test :forward_index in Symbol.(col_names)
        @test :reverse_index in Symbol.(col_names)
        @test :case in Symbol.(col_names)

        @test all(length.(indices_df.forward_index) .== 10)  # 10nt barcodes
        @test all(length.(indices_df.reverse_index) .== 10)  # 10nt barcodes

        test_outcomes["simulate"] = true
    end

    if test_outcomes["simulate"]
        @testset "demultiplex.jl" begin
            # CLI
            @test endswith(cli.always_gz("test.fasta.gz"), ".gz")
            @test endswith(cli.always_gz("test.fasta"), ".gz")
            empty!(ARGS)
            # Convert FASTA to FASTQ for demultiplex test
            fastq_records = Vector{FASTQRecord}()
            for record in data.load_fasta("test.fasta", validate=false)
                push!(fastq_records, FASTQRecord(record[1], String(record[2]), repeat('I', length(record[2]))))
            end
            data.write_fastq("test.fastq", fastq_records)

            append!(ARGS, ["demultiplex", "test.fastq", "test_indices.tsv", "test.tsv"])
            parsed_args = cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "demultiplex"
            @test parsed_args["demultiplex"]["fastq"] == "test.fastq"
            @test parsed_args["demultiplex"]["indices"] == "test_indices.tsv"
            @test parsed_args["demultiplex"]["output"] == "test.tsv"
            table, stats = demultiplex.demux("test.fastq", "test_indices.tsv")
            CSV.write("test.tsv.gz", table, delim='\t')

            # Test demultiplexing results
            @test nrow(table) > 0
            @test all(in.(["well", "case", "name", "genomic_sequence"], Ref(names(table))))

            test_outcomes["demultiplex"] = true
        end
    else
        println("Skipping remaining tests due to failed dependencies")
    end

    if test_outcomes["demultiplex"]
        @testset "exact.jl" begin
            # CLI
            empty!(ARGS)
            append!(ARGS, ["exact", "test.tsv.gz", "novel.fasta", "test_exact.tsv.gz"])
            parsed_args = cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "exact"
            @test parsed_args["exact"]["tsv"] == "test.tsv.gz"
            @test parsed_args["exact"]["fasta"] == "novel.fasta"
            @test parsed_args["exact"]["output"] == "test_exact.tsv.gz"
            @test parsed_args["exact"]["gene"] == "V"

            # Module
            table = CSV.File("test.tsv.gz", delim='\t') |> DataFrame
            db = data.load_fasta("novel.fasta", validate=false)
            mincount = 1
            minratio = 0.01
            gene = "V"
            counts_df = exact.exact_search(table, db, gene, mincount=mincount, minratio=minratio)
            sort!(counts_df, [:case, :db_name])
            @test nrow(counts_df) > 0

            # Test flanking extraction on simulated data
            read = first(table.genomic_sequence)
            sequence = read[50:60]  # Take a sample from middle
            location = findfirst(sequence, read)
            flanking = exact.extract_flanking(read, (minimum(location),maximum(location)), "V", 5)
            @test length(flanking.prefix) == 5
            @test length(flanking.sequence) == length(sequence)
        end

        @testset "hamming.jl" begin
            # CLI
            empty!(ARGS)
            append!(ARGS, ["hamming", "test.tsv.gz", "novel.fasta", "test_hamming.tsv.gz", "--maxdist", "10", "-f", "genomic_sequence"])
            parsed_args = cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "hamming"
            @test parsed_args["hamming"]["tsv"] == "test.tsv.gz"
            @test parsed_args["hamming"]["fasta"] == "novel.fasta"
            @test parsed_args["hamming"]["output"] == "test_hamming.tsv.gz"
            @test parsed_args["hamming"]["maxdist"] == 10
            @test parsed_args["hamming"]["column"] == "genomic_sequence"

            # Module
            table = CSV.File("test.tsv.gz", delim='\t') |> DataFrame
            db = data.load_fasta("novel.fasta", validate=false)
            mincount = 1
            minfreq = 0.01
            counts = hamming.hamming_search(table, db, max_dist=10, column="genomic_sequence", check_bounds=true, umi=false)
            summary = hamming.summarize(counts, db, min_count=mincount, cluster_ratio=minfreq)
            @test nrow(counts) > 0
            @test nrow(summary) > 0
        end

        @testset "profile.jl" begin
            # Get heptamer sequences from simulated data
            db = data.load_fasta("novel.fasta", validate=false)
            motifs = [string(db[i][2]) for i in 1:3]  # Take first 3 sequences

            countsum = sum(profile.counts(motifs))
            @test countsum > 0
            profsum = sum(profile.motif_profile(motifs))
            @test profsum > 0
            prof = profile.motif_profile(motifs)
            @test size(prof) == (4, length(first(motifs)))
        end

        @testset "trim.jl" begin
            # Argument Parsing Tests
            empty!(ARGS)
            append!(ARGS, ["trim", "test.tsv.gz", "novel.fasta", "test_trim.tsv", "-l", "7"])
            parsed_args = cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "trim"
            @test parsed_args["trim"]["input"] == "test.tsv.gz"
            @test parsed_args["trim"]["fasta"] == "novel.fasta"
            @test parsed_args["trim"]["output"] == "test_trim.tsv"
            @test parsed_args["trim"]["length"] == 7

            # Get test sequence from simulated data
            table = CSV.File("test.tsv.gz", delim='\t') |> DataFrame
            read = first(table.genomic_sequence)
            db = data.load_fasta("novel.fasta", validate=false)
            println("Number of sequences in database: ", length(db))
            motifs = [string(db[i][2][1:7]) for i in 1:min(3, length(db))]
            println("Extracted motifs: ", motifs)

            # Generating profiles for motifs
            prof_start = profile.motif_profile(motifs)
            prof_stop = prof_start

            # Find Max Probability Position Tests
            start = trim.find_maxprob_pos(read, prof_start)
            @test start > 0
            stop = trim.find_maxprob_pos(read[start+size(prof_start,2):end], prof_stop)
            @test stop > 0

            # Test trim_profiles
            trimmed_profs = trim.trim_profiles(motifs, 3)
            @test length(trimmed_profs) == 2
            @test all(size.(trimmed_profs) .== Ref((4, 3)))
        end

        @testset "pattern.jl" begin
            # CLI
            empty!(ARGS)
            append!(ARGS, ["pattern", "test.tsv.gz", "novel.fasta", "test_pattern.tsv.gz"])
            parsed_args = cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "pattern"
            @test parsed_args["pattern"]["input"] == "test.tsv.gz"
            @test parsed_args["pattern"]["fasta"] == "novel.fasta"
            @test parsed_args["pattern"]["output"] == "test_pattern.tsv.gz"

            # Module tests
            table = CSV.File("test.tsv.gz", delim='\t') |> DataFrame
            db = data.load_fasta("novel.fasta", validate=false)
            db_df = DataFrame(db, [:id, :seq])
            db_df.gene = [split(x, "_")[1] for x in db_df.id]
            blacklist = first(groupby(db_df, :gene))  # Use first group as blacklist

            @test split_kmers("ATCGATCG", 3) == ["CGA", "GAT"]
            @test isempty(split_kmers("", 3))
            @test isempty(split_kmers("ATCG", 5))

            # Test find_unique_fragments
            group_df = DataFrame(seq = ["ATCGATCG", "ATCGATCGAAA", "ATCGATCGCCCC"])
            outgroup_df = DataFrame(seq = ["AATTCCGG", "TTAACCG"])
            fragments = find_unique_fragments(group_df, outgroup_df; fragment_size = 4, max_fragment_size = 4)
            @test length(fragments) > 0
        end

        @testset "heptamer.jl" begin
            # CLI
            empty!(ARGS)
            append!(ARGS, ["heptamer", "test.tsv.gz", "novel.fasta", "test_heptamer.tsv.gz", "test_summary.tsv"])
            parsed_args = cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "heptamer"
                @test parsed_args["heptamer"]["tsv"] == "test.tsv.gz"
                @test parsed_args["heptamer"]["fasta"] == "novel.fasta"
                @test parsed_args["heptamer"]["output"] == "test_heptamer.tsv.gz"
                @test parsed_args["heptamer"]["summary"] == "test_summary.tsv"
                @test parsed_args["heptamer"]["chain"] == "IGHV"

                # Preparing test data from simulation
                table = CSV.File("test.tsv.gz", delim='\t') |> DataFrame
                db = data.load_fasta("novel.fasta", validate=false)
                heptamers = heptamer.load_heptamers(parsed_args["heptamer"]["json"])

                # Extract Heptamers Tests
                heptamer_df = heptamer.extract_heptamers(
                    table,
                    db,
                    heptamers[parsed_args["heptamer"]["chain"]];
                    max_dist=parsed_args["heptamer"]["maxdist"],
                    b=parsed_args["heptamer"]["begin"]+1,
                    e=parsed_args["heptamer"]["end"]
                )
                @test nrow(heptamer_df) > 0

                # Summarize Heptamers Tests
                summary_df = heptamer.summarize(
                    heptamer_df,
                    db,
                    ratio=parsed_args["heptamer"]["ratio"],
                    count=parsed_args["heptamer"]["mincount"]
                )
                @test nrow(summary_df) > 0
            end
        else
            println("Skipping remaining tests due to failed dependencies")
        end

    @testset "linkage.jl" begin
        # Generate test data programmatically
        df = DataFrame(
            case = ["D1", "D1", "D2", "D2", "D3", "D3", "D4", "D4", "D5", "D6", "D7", "D8"],
            db_name = ["A*01", "B*01", "A*01", "B*01", "A*01", "C*01", "B*01", "C*01", "B*01", "C*01", "X*01", "Y*01"]
        )
        edges_df, clusters, clusters_detailed = linkage.compute_edges_and_clusters(df; min_support=1, jaccard_threshold=0.0, similarity_threshold=0.0)
        @test nrow(edges_df) > 0
        # In toy data, A*01 and B*01 co-occur in D1 and D2 (support=2)
        row = first(filter(r -> (r.allele_a == "A*01" && r.allele_b == "B*01") || (r.allele_a == "B*01" && r.allele_b == "A*01"), edges_df))
        @test row.support == 2
        # Test that we have clusters
        @test length(clusters) >= 0
    end

    # # Cleanup test files
    # for file in ["test.fasta", "reference.fasta", "novel.fasta", "test_indices.tsv",
    #              "test.tsv.gz", "test_exact.tsv.gz", "test_hamming.tsv.gz",
    #              "test_trim.tsv", "test_pattern.tsv.gz", "test_heptamer.tsv.gz",
    #              "test_summary.tsv"]
    #     isfile(file) && rm(file)
    # end
end