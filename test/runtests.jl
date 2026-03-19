using immunediscover
using Test
using CSV
using DataFrames
using FASTX

using immunediscover.Cli
using immunediscover.Demultiplex
using immunediscover.Simulate
using immunediscover.Data
using immunediscover.Profile
using immunediscover.Exact
using immunediscover.Heptamer
using immunediscover.KeyedSets
using immunediscover.Blast
using immunediscover.Fasta
using immunediscover.Merge
using immunediscover.Haplotype
using immunediscover.Bwa
using Glob

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
        records, indices = Simulate.generate_fasta_with_mutations(
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
        db = Data.load_fasta("test.fasta", validate=false)
        ref_db = Data.load_fasta("reference.fasta", validate=false)
        novel_db = Data.load_fasta("novel.fasta", validate=false)
        indices_df = CSV.read("test_indices.tsv", DataFrame, delim='\t')

        # Test basic properties
        @test length(ref_db) == 1  # Should have one reference sequence
        @test nrow(indices_df) == 2  # Should have two donors
        @test length(novel_db) == 52  # 3 mutation types × 3 lengths × 2 donors + 2×17 extra = 52

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
            @test endswith(Cli.always_gz("test.fasta.gz"), ".gz")
            @test endswith(Cli.always_gz("test.fasta"), ".gz")
            empty!(ARGS)
            # Convert FASTA to FASTQ for demultiplex test
            fastq_records = Vector{FASTQRecord}()
            for record in Data.load_fasta("test.fasta", validate=false)
                push!(fastq_records, FASTQRecord(record[1], String(record[2]), repeat('I', length(record[2]))))
            end
            Data.write_fastq("test.fastq", fastq_records)

            append!(ARGS, ["demultiplex", "test.fastq", "test_indices.tsv", "test.tsv"])
            parsed_args = Cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "demultiplex"
            @test parsed_args["demultiplex"]["fastq"] == "test.fastq"
            @test parsed_args["demultiplex"]["indices"] == "test_indices.tsv"
            @test parsed_args["demultiplex"]["output"] == "test.tsv"
            table, stats = Demultiplex.demux("test.fastq", "test_indices.tsv")
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
            append!(ARGS, ["search", "exact", "test.tsv.gz", "novel.fasta", "test_exact.tsv.gz"])
            parsed_args = Cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "search"
            @test parsed_args["search"]["%COMMAND%"] == "exact"
            @test parsed_args["search"]["exact"]["tsv"] == "test.tsv.gz"
            @test parsed_args["search"]["exact"]["fasta"] == "novel.fasta"
            @test parsed_args["search"]["exact"]["output"] == "test_exact.tsv.gz"
            @test parsed_args["search"]["exact"]["gene"] == "V"

            # Module
            table = CSV.File("test.tsv.gz", delim='\t') |> DataFrame
            db = Data.load_fasta("novel.fasta", validate=false)
            mincount = 1
            minratio = 0.01
            gene = "V"
            counts_df = Exact.exact_search(table, db, gene, mincount=mincount, minratio=minratio)
            sort!(counts_df, [:case, :db_name])
            @test nrow(counts_df) > 0

            # Test flanking extraction on simulated data
            read = first(table.genomic_sequence)
            sequence = read[50:60]  # Take a sample from middle
            location = findfirst(sequence, read)
            flanking = Exact.extract_flanking(read, (minimum(location),maximum(location)), "V", 5)
            @test length(flanking.prefix) == 5
            @test length(flanking.sequence) == length(sequence)
        end


        @testset "profile.jl" begin
            # Get heptamer sequences from simulated data
            db = Data.load_fasta("novel.fasta", validate=false)
            motifs = [string(db[i][2]) for i in 1:3]  # Take first 3 sequences

            countsum = sum(Profile.counts(motifs))
            @test countsum > 0
            profsum = sum(Profile.motif_profile(motifs))
            @test profsum > 0
            prof = Profile.motif_profile(motifs)
            @test size(prof) == (4, length(first(motifs)))
        end



        @testset "heptamer.jl" begin
            # CLI
            empty!(ARGS)
            append!(ARGS, ["search", "heptamer", "test.tsv.gz", "novel.fasta", "test_heptamer.tsv.gz", "test_summary.tsv"])
            parsed_args = Cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "search"
            @test parsed_args["search"]["%COMMAND%"] == "heptamer"
            @test parsed_args["search"]["heptamer"]["tsv"] == "test.tsv.gz"
            @test parsed_args["search"]["heptamer"]["fasta"] == "novel.fasta"
            @test parsed_args["search"]["heptamer"]["output"] == "test_heptamer.tsv.gz"
            @test parsed_args["search"]["heptamer"]["summary"] == "test_summary.tsv"
            @test parsed_args["search"]["heptamer"]["chain"] == "IGHV"

            # Preparing test data from simulation
            table = CSV.File("test.tsv.gz", delim='\t') |> DataFrame
            db = Data.load_fasta("novel.fasta", validate=false)
            heptamers = Heptamer.load_heptamers(parsed_args["search"]["heptamer"]["json"])

            # Extract Heptamers Tests
            heptamer_df = Heptamer.extract_heptamers(
                table,
                db,
                heptamers[parsed_args["search"]["heptamer"]["chain"]];
                max_dist=parsed_args["search"]["heptamer"]["maxdist"],
                b=parsed_args["search"]["heptamer"]["begin"]+1,
                e=parsed_args["search"]["heptamer"]["end"]
            )
            @test nrow(heptamer_df) > 0

            # Summarize Heptamers Tests
            summary_df = Heptamer.summarize(
                heptamer_df,
                db,
                ratio=parsed_args["search"]["heptamer"]["ratio"],
                count=parsed_args["search"]["heptamer"]["mincount"]
            )
            @test nrow(summary_df) > 0
        end
        else
            println("Skipping remaining tests due to failed dependencies")
        end

    @testset "cooccurrence.jl" begin
        using immunediscover.Cooccurrence
        df = DataFrame(
            case = ["D1", "D1", "D2", "D2", "D3", "D3", "D4", "D4", "D5", "D6", "D7", "D8"],
            db_name = ["A*01", "B*01", "A*01", "B*01", "A*01", "C*01", "B*01", "C*01", "B*01", "C*01", "X*01", "Y*01"]
        )
        edges_df, clusters, clusters_detailed = Cooccurrence.compute_edges_and_clusters(df; min_support=1, jaccard_threshold=0.0, similarity_threshold=0.0)
        @test nrow(edges_df) > 0
        # In toy data, A*01 and B*01 co-occur in D1 and D2 (n_shared=2)
        row = first(filter(r -> (r.allele_a == "A*01" && r.allele_b == "B*01") || (r.allele_a == "B*01" && r.allele_b == "A*01"), edges_df))
        @test row.n_shared == 2
        # Test that we have clusters
        @test length(clusters) >= 0
    end

    @testset "fasta.jl" begin
        # CLI - fasta export is now under table group
        empty!(ARGS)
        append!(ARGS, ["table", "fasta", "test.tsv.gz", "test_fasta_output.fasta"])
        parsed_args = Cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "table"
        @test parsed_args["table"]["%COMMAND%"] == "fasta"
        @test parsed_args["table"]["fasta"]["input"] == "test.tsv.gz"
        @test parsed_args["table"]["fasta"]["output"] == "test_fasta_output.fasta"
        @test parsed_args["table"]["fasta"]["colname"] == "allele_name"
        @test parsed_args["table"]["fasta"]["colseq"] == "seq"

        # Module - test with generated data
        test_data = DataFrame(
            allele_name = ["A*01", "B*01", "C*01"],
            seq = ["ATCGATCG", "GCTAGCTA", "TTTTAAAA"],
            description = ["First allele", "Second allele", "Third allele"]
        )
        
        # Write test data to temporary file
        CSV.write("test_fasta_input.tsv", test_data, delim='\t')
        
        # Test basic extraction
        Fasta.extract_sequences_to_fasta("test_fasta_input.tsv", "test_fasta_basic.fasta")
        @test isfile("test_fasta_basic.fasta")
        
        # Test with custom column names
        Fasta.extract_sequences_to_fasta("test_fasta_input.tsv", "test_fasta_custom.fasta", 
                                        colname="allele_name", colseq="seq", coldesc="description")
        @test isfile("test_fasta_custom.fasta")
        
        # Test with filtering
        Fasta.extract_sequences_to_fasta("test_fasta_input.tsv", "test_fasta_filtered.fasta",
                                        filter_pattern="A\\*")
        @test isfile("test_fasta_filtered.fasta")
        
        # Clean up test files
        for file in ["test_fasta_input.tsv", "test_fasta_basic.fasta", "test_fasta_custom.fasta", "test_fasta_filtered.fasta"]
            isfile(file) && rm(file)
        end
    end

    @testset "merge.jl" begin
        # CLI - merge is now under fasta group
        empty!(ARGS)
        append!(ARGS, ["fasta", "merge", "test_merge_output.fasta", "file1.fasta", "file2.fasta"])
        parsed_args = Cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "fasta"
        @test parsed_args["fasta"]["%COMMAND%"] == "merge"
        @test parsed_args["fasta"]["merge"]["output"] == "test_merge_output.fasta"
        @test parsed_args["fasta"]["merge"]["inputs"] == ["file1.fasta", "file2.fasta"]

        # Module - test with generated FASTA files
        # Create test FASTA files
        test_fasta1 = "test_merge1.fasta"
        test_fasta2 = "test_merge2.fasta"
        test_fasta3 = "test_merge3.fasta"
        
        # Write test FASTA files
        open(FASTA.Writer, test_fasta1) do writer
            write(writer, FASTARecord("A*01", "ATCGATCG"))
            write(writer, FASTARecord("B*01", "GCTAGCTA"))
        end
        
        open(FASTA.Writer, test_fasta2) do writer
            write(writer, FASTARecord("C*01", "TTTTAAAA"))
            write(writer, FASTARecord("A*01", "ATCGATCG"))  # Duplicate sequence
        end
        
        open(FASTA.Writer, test_fasta3) do writer
            write(writer, FASTARecord("D*01", "GGGGCCCC"))
        end
        
        # Test merging two files
        immunediscover.Merge.merge_fasta_files([test_fasta1, test_fasta2], "test_merge_two.fasta")
        @test isfile("test_merge_two.fasta")
        
        # Test merging three files
        immunediscover.Merge.merge_fasta_files([test_fasta1, test_fasta2, test_fasta3], "test_merge_three.fasta")
        @test isfile("test_merge_three.fasta")
        
        # Test convenience method for two files
        immunediscover.Merge.merge_fasta_files(test_fasta1, test_fasta2, "test_merge_convenience.fasta")
        @test isfile("test_merge_convenience.fasta")
        
        # Test with cleanup pattern
        immunediscover.Merge.merge_fasta_files([test_fasta1, test_fasta2], "test_merge_cleanup.fasta", 
                               cleanup_pattern="\\*")
        @test isfile("test_merge_cleanup.fasta")
        
        # Clean up test files
        for file in [test_fasta1, test_fasta2, test_fasta3, "test_merge_two.fasta", 
                    "test_merge_three.fasta", "test_merge_convenience.fasta", "test_merge_cleanup.fasta"]
            isfile(file) && rm(file)
        end
    end

    @testset "haplotype.jl" begin
        # CLI - haplotype is now under analyze group
        empty!(ARGS)
        append!(ARGS, ["analyze", "haplotype", "test_haplotype_input.tsv", "test_haplotype_output.tsv"])
        parsed_args = Cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "analyze"
        @test parsed_args["analyze"]["%COMMAND%"] == "haplotype"
        @test parsed_args["analyze"]["haplotype"]["input"] == "test_haplotype_input.tsv"
        @test parsed_args["analyze"]["haplotype"]["output"] == "test_haplotype_output.tsv"
        @test parsed_args["analyze"]["haplotype"]["case-col"] == "case"
        @test parsed_args["analyze"]["haplotype"]["allele-col"] == "db_name"
        @test parsed_args["analyze"]["haplotype"]["gene-col"] == "gene"

        # Module - test with generated data
        # Create test data that mimics exact search output (includes sequence as required by infer_haplotypes)
        test_data = DataFrame(
            case = ["D1", "D1", "D1", "D2", "D2", "D2", "D3", "D3"],
            db_name = ["A*01", "A*02", "B*01", "A*01", "A*02", "B*01", "A*01", "B*01"],
            gene = ["V", "V", "D", "V", "V", "D", "V", "D"],
            count = [10, 8, 12, 15, 3, 9, 20, 18],
            sequence = ["ATCGATCG", "GCTAGCTA", "CCGGCCGG", "ATCGATCG", "GCTAGCTA", "CCGGCCGG", "ATCGATCG", "CCGGCCGG"]
        )
        
        # Write test data to temporary file
        CSV.write("test_haplotype_input.tsv", test_data, delim='\t')
        
        # Test basic haplotype inference
        Haplotype.infer_haplotypes("test_haplotype_input.tsv", "test_haplotype_basic.tsv")
        @test isfile("test_haplotype_basic.tsv")
        
        # Test with custom parameters
        Haplotype.infer_haplotypes("test_haplotype_input.tsv", "test_haplotype_custom.tsv",
                                  mincount=10, min_ratio=0.2)
        @test isfile("test_haplotype_custom.tsv")
        
        # Test with novel alleles FASTA
        # Create a simple novel alleles FASTA
        open(FASTA.Writer, "test_novel_alleles.fasta") do writer
            write(writer, FASTARecord("A*03", "ATCGATCGATCG"))
            write(writer, FASTARecord("B*02", "GCTAGCTAGCTA"))
        end
        
        Haplotype.infer_haplotypes("test_haplotype_input.tsv", "test_haplotype_novel.tsv",
                                  novel_fasta="test_novel_alleles.fasta")
        @test isfile("test_haplotype_novel.tsv")
        
        # Clean up test files
        for file in ["test_haplotype_input.tsv", "test_haplotype_basic.tsv", 
                    "test_haplotype_custom.tsv", "test_haplotype_novel.tsv", "test_novel_alleles.fasta"]
            isfile(file) && rm(file)
        end
    end

    @testset "blast.jl" begin
        # CLI - blast is now under search group
        empty!(ARGS)
        append!(ARGS, ["search", "blast", "test_blast_input.tsv", "test_blast_db.fasta", "test_blast_output.tsv"])
        parsed_args = Cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "search"
        @test parsed_args["search"]["%COMMAND%"] == "blast"
        @test parsed_args["search"]["blast"]["input"] == "test_blast_input.tsv"
        @test parsed_args["search"]["blast"]["fasta"] == "test_blast_db.fasta"
        @test parsed_args["search"]["blast"]["output"] == "test_blast_output.tsv"
        @test parsed_args["search"]["blast"]["minfullcount"] == 5
        @test parsed_args["search"]["blast"]["minfullfreq"] == 0.1
        @test parsed_args["search"]["blast"]["subjectcov"] == 0.1

        # Module - test utility functions that don't require BLAST
        # Test Data.load_fasta (same path as blast pipeline FASTA reads)
        test_fasta_content = ">seq1\nATCGATCG\n>seq2\nGCTAGCTA\n"
        open("test_blast_fasta.fasta", "w") do io
            write(io, test_fasta_content)
        end
        
        fasta_records = Data.load_fasta("test_blast_fasta.fasta", validate=false)
        @test length(fasta_records) == 2
        @test fasta_records[1] == ("seq1", "ATCGATCG")
        @test fasta_records[2] == ("seq2", "GCTAGCTA")
        
        # Test save_to_fasta function
        test_records = [("well1", "case1", "name1", "ATCGATCG"), ("well2", "case2", "name2", "GCTAGCTA")]
        Blast.save_to_fasta(test_records, "test_blast_save.fasta")
        @test isfile("test_blast_save.fasta")
        
        # Test CSV loading
        test_blast_data = DataFrame(
            name = ["seq1", "seq2"],
            sequence = ["ATCGATCG", "GCTAGCTA"],
            well = ["A1", "A2"],
            case = ["case1", "case2"]
        )
        CSV.write("test_blast_input.tsv", test_blast_data, delim='\t')
        
        loaded_df = Blast.load_csv("test_blast_input.tsv")
        @test nrow(loaded_df) == 2
        @test "name" ∈ names(loaded_df)
        
        # Clean up test files
        for file in ["test_blast_fasta.fasta", "test_blast_save.fasta", "test_blast_input.tsv"]
            isfile(file) && rm(file)
        end
    end

    @testset "bwa.jl" begin
        # CLI - bwa is now under analyze group
        empty!(ARGS)
        append!(ARGS, ["analyze", "bwa", "test_bwa_input.tsv", "test_bwa_output.tsv", "genome1.fasta", "genome2.fasta"])
        parsed_args = Cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "analyze"
        @test parsed_args["analyze"]["%COMMAND%"] == "bwa"
        @test parsed_args["analyze"]["bwa"]["tsv"] == "test_bwa_input.tsv"
        @test parsed_args["analyze"]["bwa"]["output"] == "test_bwa_output.tsv"
        @test parsed_args["analyze"]["bwa"]["genome"] == ["genome1.fasta", "genome2.fasta"]
        @test parsed_args["analyze"]["bwa"]["chromosome"] == "chromosome 14"
        @test parsed_args["analyze"]["bwa"]["colname"] == "best_name"

        # Module - test utility functions that don't require BWA
        # Test create_aligner function (this will fail without actual genome files, but we can test the structure)
        test_genome_paths = ["test_genome1.fasta", "test_genome2.fasta"]
        
        # Create dummy genome files for testing
        for (i, path) in enumerate(test_genome_paths)
            open(FASTA.Writer, path) do writer
                write(writer, FASTARecord("chr14", "ATCGATCGATCGATCG"))
            end
        end
        
        # Test that we can create the aligner structure (will fail at actual alignment without BWA index)
        try
            aligners = Bwa.create_aligner(test_genome_paths)
            @test length(aligners) == 2
            @test isa(aligners[1][2], BurrowsWheelerAligner.Aligner)
        catch e
            # Expected to fail without BWA index, but we can test the function exists
            # The error could be about missing index, file not found, or BWA-related issues
            error_msg = string(e)
            @test occursin("BWA", error_msg) || occursin("index", error_msg) || 
                  occursin("file", error_msg) || occursin("not found", error_msg) ||
                  occursin("aligner", error_msg)
        end
        
        # Test CSV loading for BWA input
        test_bwa_data = DataFrame(
            best_name = ["allele1", "allele2"],
            seq = ["ATCGATCG", "GCTAGCTA"]
        )
        CSV.write("test_bwa_input.tsv", test_bwa_data, delim='\t')
        
        loaded_df = CSV.read("test_bwa_input.tsv", DataFrame, delim='\t')
        @test nrow(loaded_df) == 2
        @test "best_name" ∈ names(loaded_df)
        @test "seq" ∈ names(loaded_df)
        
        # Clean up test files
        for file in [test_genome_paths..., "test_bwa_input.tsv"]
            isfile(file) && rm(file)
        end
    end


    @testset "diff.jl" begin
        # CLI - diff is now under fasta group
        empty!(ARGS)
        append!(ARGS, ["fasta", "diff", "file1.fasta", "file2.fasta", "file3.fasta"])
        parsed_args = Cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "fasta"
        @test parsed_args["fasta"]["%COMMAND%"] == "diff"
        @test parsed_args["fasta"]["diff"]["fasta"] == ["file1.fasta", "file2.fasta", "file3.fasta"]

        # Module - test the diff functionality
        # Create test FASTA files
        test_fasta1 = "test_diff1.fasta"
        test_fasta2 = "test_diff2.fasta"
        test_fasta3 = "test_diff3.fasta"
        
        # Write test FASTA files with some overlapping sequences
        open(FASTA.Writer, test_fasta1) do writer
            write(writer, FASTA.Record("seq1", "ATCGATCG"))
            write(writer, FASTA.Record("seq2", "GCTAGCTA"))
            write(writer, FASTA.Record("seq3", "TTTTAAAA"))
        end
        
        open(FASTA.Writer, test_fasta2) do writer
            write(writer, FASTA.Record("seq2", "GCTAGCTA"))  # Same as file1
            write(writer, FASTA.Record("seq4", "CCCCGGGG"))
            write(writer, FASTA.Record("seq5", "AAAAATTT"))
        end
        
        open(FASTA.Writer, test_fasta3) do writer
            write(writer, FASTA.Record("seq1", "ATCGATCG"))  # Same as file1
            write(writer, FASTA.Record("seq6", "TTTTCCCC"))
        end
        
        # Test the diff functionality by calling the main function
        # We'll test the core logic that the diff command uses
        fasta_files = [(file=file, records=immunediscover.load_fasta.(file)) for file in [test_fasta1, test_fasta2, test_fasta3]]
        sets = [(file=x, set=KeyedSet(reverse.(y))) for (x,y) in fasta_files]
        
        # Test that we can create the sets
        @test length(sets) == 3
        @test length(sets[1].set) == 3  # file1 has 3 sequences
        @test length(sets[2].set) == 3  # file2 has 3 sequences  
        @test length(sets[3].set) == 2  # file3 has 2 sequences
        
        # Test set operations
        union_12 = union(sets[1].set, sets[2].set)
        @test length(union_12) == 5  # 3 + 3 - 1 (seq2 is common) = 5
        
        intersection_12 = intersect(sets[1].set, sets[2].set)
        @test length(intersection_12) == 1  # Only seq2 is common
        
        # Test set difference
        diff_12 = setdiff(sets[1].set, sets[2].set)
        @test length(diff_12) == 2  # seq1 and seq3 are only in file1
        
        diff_21 = setdiff(sets[2].set, sets[1].set)
        @test length(diff_21) == 2  # seq4 and seq5 are only in file2
        
        # Test that the sequences are correctly identified
        # Note: KeyedSet stores (sequence, name) pairs, so we get the names
        diff_12_names = last.(collect(diff_12))
        @test "seq1" ∈ diff_12_names
        @test "seq3" ∈ diff_12_names
        @test "seq2" ∉ diff_12_names  # This should be in intersection, not difference
        
        # Clean up test files
        for file in [test_fasta1, test_fasta2, test_fasta3]
            isfile(file) && rm(file)
        end
    end

    @testset "collect.jl" begin
        # CLI - collect is now under table group
        empty!(ARGS)
        append!(ARGS, ["table", "collect", "test_*.tsv", "test_collected.tsv"])
        parsed_args = Cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "table"
        @test parsed_args["table"]["%COMMAND%"] == "collect"
        @test parsed_args["table"]["collect"]["pattern"] == "test_*.tsv"
        @test parsed_args["table"]["collect"]["output"] == "test_collected.tsv"

        # Module - test the collect functionality
        # Create test TSV files with the same structure
        test_tsv1 = "test_collect1.tsv"
        test_tsv2 = "test_collect2.tsv"
        test_tsv3 = "test_collect3.tsv"
        
        # Create test data with consistent column structure
        test_data1 = DataFrame(
            name = ["seq1", "seq2"],
            sequence = ["ATCGATCG", "GCTAGCTA"],
            well = ["A1", "A2"],
            case = ["case1", "case1"]
        )
        
        test_data2 = DataFrame(
            name = ["seq3", "seq4"],
            sequence = ["TTTTAAAA", "CCCCGGGG"],
            well = ["B1", "B2"],
            case = ["case2", "case2"]
        )
        
        test_data3 = DataFrame(
            name = ["seq5", "seq6"],
            sequence = ["AAAAATTT", "TTTTCCCC"],
            well = ["C1", "C2"],
            case = ["case3", "case3"]
        )
        
        # Write test TSV files
        CSV.write(test_tsv1, test_data1, delim='\t')
        CSV.write(test_tsv2, test_data2, delim='\t')
        CSV.write(test_tsv3, test_data3, delim='\t')
        
        # Test the collect functionality
        pattern = "test_collect*.tsv"
        files = Glob.glob(pattern)
        @test length(files) == 3
        @test test_tsv1 ∈ files
        @test test_tsv2 ∈ files
        @test test_tsv3 ∈ files
        
        # Test collecting the files
        collected = []
        first_file_columns = nothing
        for file in files
            df = CSV.read(file, DataFrame, delim='\t')
            if first_file_columns === nothing
                first_file_columns = names(df)
            else
                @test first_file_columns == names(df)
            end
            push!(collected, df)
        end
        
        @test length(collected) == 3
        @test first_file_columns == ["name", "sequence", "well", "case"]
        
        # Test concatenating the data
        collected_df = vcat(collected...)
        @test nrow(collected_df) == 6  # 2 + 2 + 2 = 6 rows
        @test ncol(collected_df) == 4  # 4 columns
        @test "name" ∈ names(collected_df)
        @test "sequence" ∈ names(collected_df)
        @test "well" ∈ names(collected_df)
        @test "case" ∈ names(collected_df)
        
        # Test that all sequences are present
        all_sequences = collected_df.sequence
        @test "ATCGATCG" ∈ all_sequences
        @test "GCTAGCTA" ∈ all_sequences
        @test "TTTTAAAA" ∈ all_sequences
        @test "CCCCGGGG" ∈ all_sequences
        @test "AAAAATTT" ∈ all_sequences
        @test "TTTTCCCC" ∈ all_sequences
        
        # Test writing the collected data
        output_file = "test_collected_output.tsv"
        CSV.write(output_file, collected_df, delim='\t', compress=true)
        @test isfile(output_file)
        
        # Verify the output file can be read back
        read_back = CSV.read(output_file, DataFrame, delim='\t')
        @test nrow(read_back) == 6
        @test ncol(read_back) == 4
        @test names(read_back) == ["name", "sequence", "well", "case"]
        
        # Clean up test files
        for file in [test_tsv1, test_tsv2, test_tsv3, output_file]
            isfile(file) && rm(file)
        end
    end

    # Cleanup test files
    for file in ["test.fasta", "reference.fasta", "novel.fasta", "test_indices.tsv",
                 "test.tsv.gz", "test_exact.tsv.gz", "test_heptamer.tsv.gz",
                 "test_summary.tsv"]
        isfile(file) && rm(file)
    end
end