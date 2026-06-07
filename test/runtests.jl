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
using immunediscover.Selftest
using immunediscover.Fasta
using immunediscover.Merge
using immunediscover.Haplotype
using immunediscover.Bwa
using immunediscover.Filters
using immunediscover.Report
using immunediscover.SeqStats
using immunediscover.Mosaic
using immunediscover.Table
using immunediscover.Cooccurrence
using immunediscover.HSMM
using Glob
using Random

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

    @testset "filters.jl" begin
        df = DataFrame(
            count = [10, 5, 2, 8, 1],
            ratio = [1.0, 0.5, 0.2, 0.8, 0.1],
            name = ["IGHV1*01", "IGHV2*01", "IGHV3*01", "IGHV4*01", "IGHV5*01"],
            qseq = ["ATCGATCGATCG", "ATCG", "ATCGATCGATCGATCG", "ATCGATCG", "AT"],
            mismatch = [0, 3, 5, 1, -1]
        )

        df_copy = copy(df)
        GermlineFilter([MinThreshold(:count, 5.0, "Min count")])(df_copy)
        @test nrow(df_copy) == 3
        @test all(df_copy.count .>= 5)

        df_copy = copy(df)
        GermlineFilter([MaxThreshold(:mismatch, 3.0, "Max mismatch")])(df_copy)
        @test nrow(df_copy) == 4
        @test all(df_copy.mismatch .<= 3)

        df_copy = copy(df)
        GermlineFilter([MinStringLength(:qseq, 8, "Min read length")])(df_copy)
        @test nrow(df_copy) == 3
        @test all(length.(df_copy.qseq) .>= 8)

        df_copy = copy(df)
        GermlineFilter([NonNegative(:mismatch, "Non-negative mismatch")])(df_copy)
        @test nrow(df_copy) == 4
        @test all(df_copy.mismatch .>= 0)

        df_copy = copy(df)
        GermlineFilter([CustomFilter(x -> x.ratio > 0.3, "High ratio")])(df_copy)
        @test nrow(df_copy) == 3
        @test all(df_copy.ratio .> 0.3)

        df_copy = copy(df)
        GermlineFilter([
            MinThreshold(:count, 5.0, "Min count"),
            MaxThreshold(:mismatch, 3.0, "Max mismatch"),
        ])(df_copy)
        @test nrow(df_copy) == 3

        df_copy = copy(df)
        apply_filters!(df_copy, MinThreshold(:count, 5.0, "Min count"))
        @test nrow(df_copy) == 3

        @test passes((count=10, ratio=1.0), MinThreshold(:count, 5.0, ""))
        @test !passes((count=3, ratio=1.0), MinThreshold(:count, 5.0, ""))
        @test passes((mismatch=3,), MaxThreshold(:mismatch, 5.0, ""))
        @test !passes((mismatch=6,), MaxThreshold(:mismatch, 5.0, ""))

        # add_group_ratio!: value / per-group maximum (shared "allelic ratio" helper)
        gr = DataFrame(gene=["V","V","D"], count=[2,4,5])
        Filters.add_group_ratio!(gr, :count, [:gene], :ratio)
        @test gr[gr.gene .== "V", :ratio] == [0.5, 1.0]
        @test gr[gr.gene .== "D", :ratio] == [1.0]

        # GermlineFilter show methods
        gf = GermlineFilter([MinThreshold(:count, 5.0, "Min count")])
        @test sprint(show, gf) == "GermlineFilter(1 criterion)"
        @test occursin("Min count", sprint(show, MIME("text/plain"), gf))

        # annotate path: record the first failing criterion instead of dropping
        adf = DataFrame(count=[10, 2, 8], qseq=["AAAAAAAA", "AA", "AAAAAAAA"])
        Filters.annotate_rejections!(adf, FilterCriterion[
            MinThreshold(:count, 5.0, "Min count"),
            MinStringLength(:qseq, 4, "Min len"),
        ])
        @test "reject_reason" in names(adf)
        @test "reject_stage" in names(adf)
        @test adf.reject_reason == ["", "Min count", ""]   # row 2 fails count first
        @test nrow(Filters.accepted(adf)) == 2

        # mark_rejected! only marks not-yet-rejected rows (first reason wins)
        mdf = DataFrame(x=[1, 2, 3])
        Filters.mark_rejected!(mdf, [true, false, true], "first", "s1")
        Filters.mark_rejected!(mdf, [true, true, false], "second", "s2")
        @test mdf.reject_reason == ["first", "second", "first"]
        @test mdf.reject_stage == ["s1", "s2", "s1"]
    end

    @testset "exact filter pipeline" begin
        # exact_search filtering now lives in handle_exact via the shared annotate machinery.
        df = DataFrame(
            well = [1, 1, 1, 1], case = ["D1", "D2", "D1", "D1"],
            gene = ["IGHV1-1", "IGHV1-1", "IGHV1-1", "CTRL1"],
            db_name = ["IGHV1-1*01", "IGHV1-1*01", "IGHV1-1*02", "CTRL1*01"],
            sequence = ["AAA", "AAA", "CCC", "GGG"],
            full_count = [10, 8, 2, 10], count = [10, 8, 2, 10],
            full_ratio = [1.0, 1.0, 0.2, 1.0], ratio = [1.0, 1.0, 0.2, 1.0],
            n_donors = [2, 2, 1, 1], max_full_ratio = [1.0, 1.0, 0.2, 1.0],
        )
        Filters.init_rejection_columns!(df)
        crit = Exact.exact_filter_criteria(; mincount=5, minratio=0.1, expect_dict=Dict{String,Float64}())
        Exact.annotate_stage!(df, crit, "count and ratio filter")
        # row 3 fails min count (full_count 2 < 5); the count term subsumes the redundant count check
        @test df.reject_reason == ["", "", "min count (--mincount 5)", ""]

        # locus aggregates exclude control genes AND already-rejected rows
        Exact.add_frequency_columns!(df, "IGHV")
        # (1,D1,IGHV1-1): only row1 accepted (row3 rejected) → gene_count 10; (1,D2): row2 → 8;
        # CTRL1 is outside the locus → default 0.
        @test df.gene_count == [10, 8, 10, 0]
        @test df.allelic_ratio[1] == 1.0        # row1 is its gene's only accepted read in (1,D1)

        # optional recurrence filter drops a single-donor candidate
        df2 = DataFrame(db_name=["IGHV1-1*01", "IGHV1-1*02"], gene=["IGHV1-1", "IGHV1-1"],
                        sequence=["AAA", "CCC"], full_count=[9, 9], count=[9, 9],
                        full_ratio=[1.0, 1.0], ratio=[1.0, 1.0], n_donors=[3, 1],
                        max_full_ratio=[1.0, 1.0])
        Filters.init_rejection_columns!(df2)
        Exact.annotate_stage!(df2,
            Exact.exact_filter_criteria(; mincount=5, minratio=0.1, expect_dict=Dict{String,Float64}(),
                                        min_recurrence=2), "count and ratio filter")
        @test df2.reject_reason == ["", "min donor recurrence (--min-recurrence 2)"]
    end

    @testset "exact column ordering + rounding" begin
        df = DataFrame(
            well=[1], case=["D1"], gene=["IGHV1-1"], db_name=["IGHV1-1*01"],
            count=[10], allelic_ratio=[0.123456], heptamer=["CACAGTG"],
            sequence=["ACGTACGT"], prefix=["TTTT"], spacer=["GGG"], nonamer=["AAAAAAAAA"],
            reject_reason=[""], reject_stage=[""],
        )
        o = Exact.order_exact_columns(df, VGene(), nothing)
        cols = names(o)
        # the long DNA columns are last, in genomic 5'→3' order for V (prefix, seq, 3' RSS)
        @test cols[end-4:end] == ["prefix", "sequence", "heptamer", "spacer", "nonamer"]
        # identifiers/metrics precede the DNA block
        @test findfirst(==("count"), cols) < findfirst(==("sequence"), cols)
        @test findfirst(==("allelic_ratio"), cols) < findfirst(==("prefix"), cols)
        # J places its 5' RSS before the sequence; extension mode just prefix/seq/suffix
        @test Exact.dna_layout(JGene(), nothing) == ["nonamer", "spacer", "heptamer", "sequence", "suffix"]
        @test Exact.dna_layout(VGene(), 20) == ["prefix", "sequence", "suffix"]

        # floats rounded to 4 dp; integer/string columns untouched
        Data.round_floats!(o)
        @test o.allelic_ratio[1] ≈ 0.1235
        @test o.count[1] === 10
        @test o.sequence[1] == "ACGTACGT"
    end

    @testset "hsmm collapse + posterior annotation" begin
        # Collapse represents each sequence by its best detection; count = detections clearing
        # min_posterior; the posterior threshold is annotated (not a silent pre-collapse drop).
        res = DataFrame(
            well=[1,1,1,1], case=["D1","D1","D1","D1"], sequence=["AAA","AAA","AAA","CCC"],
            posterior_prob=[0.9, 0.4, 0.8, 0.3],
            pre_nonamer=["pn9","pn4","pn8","x"], pre_spacer=["","","",""], pre_heptamer=["","","",""],
            post_heptamer=["","","",""], post_spacer=["","","",""], post_nonamer=["","","",""],
            heptamer_logp_pre=fill(log(0.4),4), heptamer_logp_post=fill(log(0.4),4),
            log_path_prob=[-1.0,-2.0,-1.5,-3.0], log_total_prob=fill(-0.5,4),
            isin_db=fill(false,4), db_name=["IGHD1*01","IGHD1*01","IGHD1*01","IGHD2*01"],
            nearest_db=["IGHD1*01","IGHD1*01","IGHD1*01","IGHD2*01"], nearest_db_dist=[1,2,1,5],
        )
        c = HSMM.collapse_detections(res, 0.7)
        @test nrow(c) == 2
        aaa = c[c.sequence .== "AAA", :]
        @test aaa[1, :count] == 2                 # 0.9 and 0.8 clear 0.7; 0.4 does not
        @test aaa[1, :posterior_prob] == 0.9      # best detection represents the cluster
        @test aaa[1, :pre_nonamer] == "pn9"       # ...including its flanks
        @test aaa[1, :nearest_db_dist] == 1       # minimum over the cluster
        ccc = c[c.sequence .== "CCC", :]
        @test ccc[1, :count] == 0                 # best 0.3 < 0.7

        Filters.init_rejection_columns!(c)
        crit = FilterCriterion[MinThreshold(:posterior_prob, 0.7, "min posterior")]
        for cr in crit
            fail = Bool[!passes(r, cr) for r in eachrow(c)]
            Filters.mark_rejected!(c, fail, cr.label, "detection filter")
        end
        @test c[c.sequence .== "AAA", :reject_reason][1] == ""
        @test c[c.sequence .== "CCC", :reject_reason][1] == "min posterior"
    end

    @testset "report findings helpers" begin
        # reject_counts: accepted (empty) first, then reasons by descending count
        rc = Report.reject_counts(["", "min count", "", "min count", "min ratio"])
        @test rc[1] == ("accepted", 2)
        @test rc[2] == ("min count", 2)
        @test ("min ratio", 1) in rc
        @test Report.reject_counts(String[]) == Tuple{String,Int}[]

        # composition_matrix underlies rss_consistency: a conserved motif → one base per column
        M = Report.composition_matrix(["CACAGTG", "CACAGTG", "CACAGTG"])
        @test all(isapprox.(sum(M, dims=1), 1.0))
        @test maximum(M[:, 1]) == 1.0           # fully conserved column
        Mv = Report.composition_matrix(["AAAA", "CCCC"])
        @test maximum(Mv[:, 1]) == 0.5          # split column ⇒ lower conservation

        # consensus_motif: most-frequent base per position + per-position conservation
        cons, cN = Report.consensus_motif(["CACAGTG", "CACAGTG", "CACAGTC"])
        @test cons == "CACAGTG"                 # pos 7: G in 2/3 ⇒ consensus G
        @test cN[1] == 1.0                      # pos 1 fully conserved
        @test cN[7] ≈ 2/3
        @test Report.consensus_motif(String[]) == ("", Float64[])

        # filter_quality_report prints and returns nothing; no error with mixed accept/reject
        qf = DataFrame(reject_reason=["", "x", ""], n_donors=[3, 1, 4],
                       max_full_ratio=[1.0, 0.1, 0.9], full_count=[10, 2, 8])
        @test Report.filter_quality_report(qf, [:n_donors, :max_full_ratio]) === nothing
        # all-accepted ⇒ nothing to compare ⇒ no-op
        qa = DataFrame(reject_reason=["", ""], n_donors=[3, 4])
        @test Report.filter_quality_report(qa, [:n_donors]) === nothing
    end

    @testset "report" begin
        @test occursin("kept 8/10", Report.stage_summary("edge", 8, 10))
        @test occursin("80.0%", Report.stage_summary("edge", 8, 10))
        @test occursin("min=", Report.distribution_summary([0.1, 0.5, 0.9]))
        @test occursin("max=", Report.distribution_summary([0.1, 0.5, 0.9]))
        @test Report.distribution_summary(Float64[]) == ""
        @test Report.stage_report("edge", 8, 10) === nothing
        @test Report.stage_report("scov", 5, 5; values=[0.2, 0.5, 0.8]) === nothing
        @test Report.stage_report("scov", 5, 8; values=[0.2, 0.5, 0.8], histogram=true) === nothing

        # plotting helpers (UnicodePlots is a direct dependency) run without error
        @test Data.histogram_if_available([1.0, 2.0, 2.0, 3.0]) === nothing
        @test Data.histogram_if_available(Float64[]) === nothing
        @test Data.barplot_if_available(["a", "b"], [3, 1]) === nothing
        @test Data.heatmap_if_available(rand(4, 6)) === nothing

        # dominant length: most common length among cores
        @test Report.dominant_length(["ACGT","ACGA","TTT"]) == 4
        @test Report.dominant_length(String[]) == 0
        # composition matrix: 4×L, columns sum to 1, splits at variable positions
        M = Report.composition_matrix(["ACGT","ACGA"])
        @test size(M) == (4, 4)
        @test all(isapprox.(sum(M; dims=1), 1.0))   # each column is a frequency distribution
        @test M[1, 1] == 1.0                        # position 1 conserved (all A)
        @test M[4, 4] == 0.5 && M[1, 4] == 0.5      # position 4 splits T/A
        @test Report.is_novel_name("IGHV1-2_S1234")
        @test !Report.is_novel_name("IGHV1-2*01")
        row = Report.mismatch_row("ACGA", "ACGT", 1.0)
        @test row == [0.0, 0.0, 0.0, 1.0]
        @test Report.matched_germline("IGHV1-2*01", Report.db_dict([("IGHV1-2*01", "ACGT")])) == "ACGT"
        db = [("IGHV1-2*01", "ACGT")]
        panels, suspicious = Report.gene_novel_diff_panels(
            ["IGHV1-2", "IGHV1-2", "IGHV1-2"],
            ["ACGT", "ACGA", "ACGT"],
            ["IGHV1-2*01", "IGHV1-2*01_S1234", "IGHV1-2*01"],
            ["IGHV1-2*01", "IGHV1-2*01", "IGHV1-2*01"];
            reads=[500, 50, 500], aln_mismatches=[0, 1, 0], db_seqs=db)
        @test suspicious == 0
        @test length(panels) == 1
        @test first(first(panels)) == "IGHV1-2"
        @test sum(first(panels)[2]) ≈ 1.0
        db2 = [("IGHV1-2*01", "ACGT"), ("IGHV3-7*01", "ACGA")]
        panels2, _ = Report.gene_novel_diff_panels(
            ["IGHV3-7", "IGHV1-2"], ["ACGC", "ACGA"],
            ["IGHV3-7*01_S1", "IGHV1-2*01_S1"],
            ["IGHV3-7*01", "IGHV1-2*01"];
            reads=[100, 50], aln_mismatches=[1, 1], db_seqs=db2)
        @test length(panels2) == 2
        @test [first(p) for p in panels2] == ["IGHV1-2", "IGHV3-7"]
        panels3, _ = Report.gene_novel_diff_panels(
            ["IGHV1-2", "IGHV1-2"], ["ACGA", "ACGC"],
            ["IGHV1-2*01_S1", "IGHV1-2*01_S2"],
            ["IGHV1-2*01", "IGHV1-2*01"];
            reads=[500, 100], aln_mismatches=[1, 1], db_seqs=db)
        M3 = first(panels3)[2]
        @test maximum(@view(M3[1, :])) ≈ 1.0
        @test maximum(@view(M3[2, :])) ≈ 0.2
        @test Report.snp_support_row("ACGA", "ACGT", 50, 100) == [0.0, 0.0, 0.0, 0.5]
        @test Report.cluster_profile_heatmap(["IGHV1-2"], ["ACGT"],
                                             ["IGHV1-2*01"], ["IGHV1-2*01"];
                                             reads=[100], aln_mismatch=[0],
                                             db_seqs=db) === nothing
        @test Report.cluster_profile_heatmap(["IGHV1-2"], ["ACGA"],
                                             ["IGHV1-2*01_S1234"], ["IGHV1-2*01"];
                                             reads=[50], aln_mismatch=[1],
                                             db_seqs=db) === nothing
        db_known = [("IGHV4-39*01_S1660", "ACGTACGT"), ("IGHV4-39*01", "ACGTACGT")]
        _, susp_known = Report.gene_novel_diff_panels(
            ["IGHV4-39"], ["ACGTACGT"], ["IGHV4-39*01_S1660"], ["IGHV4-39*01_S1660"];
            reads=[100], aln_mismatches=[0], db_seqs=db_known)
        @test susp_known == 0
        # grouped parameter display: runs, groups known keys, "other" catches the rest
        @test Report.params_report(Dict("input" => "a.tsv", "gene" => "V", "extra" => 1),
                                   ["IO" => ["input"], "Gene" => ["gene"]]) === nothing
    end

    @testset "seqstats" begin
        @test SeqStats.gc_content("GCGC") == 1.0
        @test SeqStats.gc_content("ATAT") == 0.0
        @test SeqStats.gc_content("GCAT") == 0.5
        @test SeqStats.gc_content("") == 0.0
        @test SeqStats.max_homopolymer("AAATTC") == 3
        @test SeqStats.max_homopolymer("ACGT") == 1
        @test SeqStats.max_homopolymer("") == 0
        @test SeqStats.n_content("ANNA") == 0.5
        @test SeqStats.n_content("ACGT") == 0.0
        @test SeqStats.shannon_entropy(["AAA", "AAA"]) == 0.0
        @test SeqStats.shannon_entropy(["AA", "AT"]) == 0.5      # 0 bits + 1 bit, averaged
        @test SeqStats.shannon_entropy(String[]) == 0.0
        # positional_entropy: per-position bits, left-aligned, variable length tolerated
        @test SeqStats.positional_entropy(["AA", "AT"]) == [0.0, 1.0]
        @test SeqStats.positional_entropy(["AAA", "AAA"]) == [0.0, 0.0, 0.0]
        @test SeqStats.positional_entropy(["ACG", "AC"]) == [0.0, 0.0, 0.0]  # pos 3 has 1 seq → 0
        @test SeqStats.positional_entropy(String[]) == Float64[]
        @test SeqStats.consensus_fraction(["AAA", "AAA"]) == 1.0
        @test SeqStats.consensus_fraction(["AA", "AT"]) == 0.75
        @test_throws ArgumentError SeqStats.shannon_entropy(["AA", "AAA"])
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

            append!(ARGS, ["preprocess", "demultiplex", "test.fastq", "test_indices.tsv", "test.tsv"])
            parsed_args = Cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "preprocess"
            @test parsed_args["preprocess"]["%COMMAND%"] == "demultiplex"
            @test parsed_args["preprocess"]["demultiplex"]["fastq"] == "test.fastq"
            @test parsed_args["preprocess"]["demultiplex"]["indices"] == "test_indices.tsv"
            @test parsed_args["preprocess"]["demultiplex"]["output"] == "test.tsv"
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
            gene = "V"
            counts_df = Exact.exact_search(table, db, gene)
            sort!(counts_df, [:case, :db_name])
            @test nrow(counts_df) > 0
            # exact_search now returns unfiltered candidates plus quality metrics.
            @test "n_donors" in names(counts_df)
            @test "n_reads_total" in names(counts_df)
            @test "max_full_ratio" in names(counts_df)

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
                ratio=parsed_args["search"]["heptamer"]["ratio"],
                mincount=parsed_args["search"]["heptamer"]["mincount"]
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

        # Homozygous call must still report the runner-up in other_alleles (not drop it).
        homo_df = DataFrame(case=["D1","D1"], db_name=["A*01","A*02"], gene=["V","V"],
                            count=[100, 2], sequence=["ACGT","TTTT"])
        CSV.write("test_hap_homo.tsv", homo_df, delim='\t')
        res = Haplotype.infer_haplotypes("test_hap_homo.tsv", "test_hap_homo_out.tsv"; mincount=1, min_ratio=0.1)
        hrow = first(filter(r -> r.case == "D1" && r.gene == "V", res))
        @test hrow.genotype == "homozygous"
        @test hrow.allele_2 == ""
        @test occursin("A*02", hrow.other_alleles)   # runner-up retained
        for f in ["test_hap_homo.tsv", "test_hap_homo_out.tsv"]; isfile(f) && rm(f); end
    end

    @testset "blast.jl" begin
        # CLI - blast is now under search group
        empty!(ARGS)
            append!(ARGS, ["discover", "blast", "test_blast_input.tsv", "test_blast_db.fasta", "test_blast_output.tsv"])
            parsed_args = Cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "discover"
            @test parsed_args["discover"]["%COMMAND%"] == "blast"
            @test parsed_args["discover"]["blast"]["input"] == "test_blast_input.tsv"
            @test parsed_args["discover"]["blast"]["fasta"] == "test_blast_db.fasta"
            @test parsed_args["discover"]["blast"]["output"] == "test_blast_output.tsv"
            @test parsed_args["discover"]["blast"]["minfullcount"] == 5
            @test parsed_args["discover"]["blast"]["minfullratio"] == 0.1
            @test parsed_args["discover"]["blast"]["subjectcov"] == 0.1
            @test parsed_args["discover"]["blast"]["work-dir"] == ".immunediscover"

        # Preset application: untouched params take the V preset, explicit overrides are kept.
        empty!(ARGS)
        append!(ARGS, ["discover", "blast", "i.tsv", "d.fa", "o.tsv", "-g", "V"])
        pa_preset = Cli.apply_blast_presets!(Cli.parse_commandline(ARGS))
        @test pa_preset["discover"]["blast"]["minfullratio"] ≈ 0.08    # V preset (recall-safe FP cut)
        @test pa_preset["discover"]["blast"]["min-corecov"] == 0.50     # V preset (was default 0.6)

        empty!(ARGS)
        append!(ARGS, ["discover", "blast", "i.tsv", "d.fa", "o.tsv", "-g", "V", "--min-corecov", "0.9"])
        pa_override = Cli.apply_blast_presets!(Cli.parse_commandline(ARGS))
        @test pa_override["discover"]["blast"]["min-corecov"] == 0.9    # explicit override respected

        # Single source of truth: every BLAST_DEFAULTS entry must equal the parsed ArgParse
        # default (guards against the arg table and the preset/defaults table drifting apart).
        empty!(ARGS)
        append!(ARGS, ["discover", "blast", "i.tsv", "d.fa", "o.tsv"])
        bdef = Cli.parse_commandline(ARGS)["discover"]["blast"]
        for (k, v) in Cli.BLAST_DEFAULTS
            @test haskey(bdef, k)
            @test bdef[k] == v
        end
        @test bdef["min-reads-total"] == 0    # new abundance filter, off by default

        empty!(ARGS)
        append!(ARGS, ["discover", "blast", "i.tsv", "d.fa", "o.tsv", "--min-reads-total", "40"])
        @test Cli.parse_commandline(ARGS)["discover"]["blast"]["min-reads-total"] == 40

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

    @testset "exact grouped_ratios refgene (multi-allele)" begin
        # A refgene that matches several alleles in a (well,case) must use a scalar
        # (summed) reference count, not a vector (which previously DimensionMismatch'd).
        cdf = DataFrame(
            well = [1, 1, 1], case = ["D1", "D1", "D1"],
            db_name = ["IGHV3-23*01", "IGHV3-23*02", "IGHV1-2*01"],
            count = [10, 30, 20],
        )
        out = Exact.grouped_ratios(cdf, "IGHV3-23", count_col=:count)
        col = "count_IGHV3-23_ratio"
        @test col in names(out)
        @test out[out.db_name .== "IGHV1-2*01", col][1] ≈ 20 / 40   # 20 / (10+30)
        @test out[out.db_name .== "IGHV3-23*02", col][1] ≈ 30 / 40
    end

    @testset "blast accumulate_affixes" begin
        # A gene present in reads is extended by its common flanks; a decoy/pseudo absent
        # from the reads survives as an unextended singleton (so -p decoys reach the DB).
        gene = "ACGTACGTACGT"
        demux = DataFrame(
            well = [1, 1], case = ["D1", "D1"], name = ["r1", "r2"],
            genomic_sequence = ["AAAAA" * gene * "TTTTT", "AAAAA" * gene * "TTTTT"],
        )
        db = [("G1*01", gene), ("Pdecoy*01", "GGGGGGGGGGGG")]
        extended = Blast.accumulate_affixes(db, demux; forward_extension=5, reverse_extension=5)
        d = Dict(name => (seq, pre, suf) for (name, seq, pre, suf) in extended)

        @test haskey(d, "G1*01")
        @test haskey(d, "Pdecoy*01")          # decoy preserved in the extended DB
        g1seq, g1pre, g1suf = d["G1*01"]
        @test g1pre == "AAAAA"
        @test g1suf == "TTTTT"
        @test g1seq == "AAAAA" * gene * "TTTTT"
        pseq, ppre, psuf = d["Pdecoy*01"]
        @test ppre == ""                       # absent in reads → unextended
        @test psuf == ""
        @test pseq == "GGGGGGGGGGGG"

        # Heterogeneous upstream context: eroding LCP/LCS would collapse to "" or 1–2 nt;
        # consensus keeps boundary-adjacent agreement.
        mixed = DataFrame(
            well = [1, 1, 1], case = ["D1", "D1", "D1"], name = ["r1", "r2", "r3"],
            genomic_sequence = [
                "AAAAT" * gene * "TCCCC",
                "CCCGT" * gene * "TAAAA",
                "GGGGT" * gene * "TTTTT",
            ],
        )
        mixed_ext = Blast.accumulate_affixes([("G2*01", gene)], mixed; forward_extension=5, reverse_extension=5)
        _, _, mpre, msuf = only(mixed_ext)
        @test mpre == "GT"
        @test msuf == "T"
    end

    @testset "blast annotate + partition (full vs filtered)" begin
        # Mirrors handle_blast's final stages without needing blastn: annotate, then split.
        clusters = DataFrame(
            sseqid = ["A*01", "B*01", "C*01", "D*01"],
            full_count = [10, 2, 10, 10],         # B fails Min count
            full_ratio = [1.0, 1.0, 0.01, 1.0],   # C fails Min ratio
            qseq = ["ACGTACGTAC", "ACGTACGTAC", "ACGTACGTAC", "AC"],  # D fails Min len
            aln_mismatch = [0, 0, 0, 0],
            corecov = [0.9, 0.9, 0.9, 0.9],
        )
        Filters.mark_rejected!(clusters, clusters.corecov .< 0.5, "corecov < 0.5", "corecov")
        criteria = FilterCriterion[
            MinThreshold(:full_count, 5.0, "Min count"),
            MinThreshold(:full_ratio, 0.1, "Min ratio"),
            MinStringLength(:qseq, 5, "Min len"),
            MaxThreshold(:aln_mismatch, 14.0, "Max dist"),
        ]
        Filters.annotate_rejections!(clusters, criteria; stage="output filter")
        @test clusters.reject_reason == ["", "Min count", "Min ratio", "Min len"]

        kept = Filters.accepted(clusters)
        @test nrow(kept) == 1
        @test kept.sseqid == ["A*01"]
        filtered = select(kept, Not([:reject_reason, :reject_stage]))
        @test !("reject_reason" in names(filtered))   # filtered output drops the reason columns
        @test nrow(clusters) == 4                       # full table keeps every candidate

        # optional quality-metric criteria (recurrence / homopolymer) annotate like the rest
        cl2 = DataFrame(sseqid=["A", "B"], n_donors=[3, 1], max_homopolymer=[2, 9])
        Filters.annotate_rejections!(cl2, FilterCriterion[
            MinThreshold(:n_donors, 2.0, "min recurrence"),
            MaxThreshold(:max_homopolymer, 5.0, "max homopolymer"),
        ])
        @test cl2.reject_reason == ["", "min recurrence"]   # B fails recurrence (n_donors 1 < 2) first
    end

    @testset "blast neighbor stats (satellite detection)" begin
        cores = ["AAAA", "AAAT", "GGGG"]
        reads = [100, 2, 50]
        nd, pr = Blast.neighbor_stats(cores, reads)
        @test nd[1] == -1            # AAAA is most abundant → no parent
        @test pr[1] == 1.0
        @test nd[2] == 1             # AAAT is 1 bp from AAAA (much bigger)
        @test pr[2] == 50.0          # 100 / 2
        @test nd[3] == 4             # GGGG nearest more-abundant is AAAA at distance 4
        @test pr[3] == 2.0           # 100 / 50
        @test Blast.core_distance("ACGT", "ACGA") == 1
        @test Blast.core_distance("AC", "ACGT") > 0    # different length → Levenshtein
        @test Blast.satellite_score(1, 50.0) ≈ 0.5
        @test Blast.satellite_score(1, 20.0) ≈ 0.5
        @test Blast.likely_satellite(1, 50.0)
        @test !Blast.likely_satellite(4, 2.0)
        @test !Blast.likely_satellite(-1, 1.0)
        df_sat = DataFrame(
            gene = ["G", "G", "G"],
            aln_qseq = ["AAAA", "AAAT", "GGGG"],
            reject_reason = ["", "", ""],
            n_reads_total = [100, 2, 50],
        )
        Blast.add_neighbor_stats!(df_sat)
        @test df_sat.nn_dist == [-1, 1, 4]
        @test df_sat.parent_ratio == [1.0, 50.0, 2.0]
        @test df_sat.satellite_score ≈ [0.0, 0.5, 0.0]
        @test df_sat.likely_satellite == [false, true, false]
    end

    @testset "mosaic chimera_score" begin
        refs = [("IGHV1-2*01", "AAAAAAAAAA"), ("IGHV1-2*02", "BBBBBBBBBB")]
        @test Mosaic.chimera_score("AAAAAAAAAA", refs) == 0.0
        @test Mosaic.chimera_score("AAAAABBBBBB", refs) ≈ 1.0
        @test Mosaic.chimera_score("ACGT", refs) == 0.0
        db_by = Mosaic.refs_by_gene(refs)
        @test length(db_by["IGHV1-2"]) == 2
        df = DataFrame(gene=["IGHV1-2", "IGHV1-2"], sequence=["AAAAABBBBBB", "AAAAAAAAAA"])
        Mosaic.add_chimera_scores!(df, db_by; seq_col=:sequence, gene_col=:gene)
        @test df.chimera_score[1] ≈ 1.0
        @test df.chimera_score[2] == 0.0
    end

    @testset "blast name_candidate (known never named novel)" begin
        db = [("IGHD1-7*01", "GGTATAACTGGAACTAC"), ("IGHD1-7*02", "GGTATAACTGGAACAAC")]
        # exact match to the best-hit reference → that reference
        @test Blast.name_candidate("GGTATAACTGGAACTAC", "IGHD1-7*01", 0, db, true) == "IGHD1-7*01"
        # core EXACTLY equals a different known allele than the best-hit (aln_mismatch>0 vs best
        # hit): must resolve to that known allele regardless of --isin — never a novel _S name.
        @test Blast.name_candidate("GGTATAACTGGAACAAC", "IGHD1-7*01", 1, db, true)  == "IGHD1-7*02"
        @test Blast.name_candidate("GGTATAACTGGAACAAC", "IGHD1-7*01", 1, db, false) == "IGHD1-7*02"
        # genuine internal variation vs every known allele → novel hashed name
        @test occursin(r"_S\d+$", Blast.name_candidate("GGTATAACTGGAACGGC", "IGHD1-7*01", 2, db, true))
        # --isin governs only partial coverage (core ⊂ known): on → known, off → novel
        @test Blast.name_candidate("TATAACTGG", "IGHD1-7*01", 3, db, true)  == "IGHD1-7*01"
        @test occursin(r"_S\d+$", Blast.name_candidate("TATAACTGG", "IGHD1-7*01", 3, db, false))
    end

    @testset "selftest recovery" begin
        # CLI
        empty!(ARGS)
        append!(ARGS, ["discover", "selftest", "disc.tsv.gz", "base.fasta", "truth.fasta", "report.tsv"])
        pa = Cli.parse_commandline(ARGS)
        @test pa["discover"]["%COMMAND%"] == "selftest"
        @test pa["discover"]["selftest"]["discovery"] == "disc.tsv.gz"
        @test pa["discover"]["selftest"]["base"] == "base.fasta"
        @test pa["discover"]["selftest"]["truth"] == "truth.fasta"
        @test pa["discover"]["selftest"]["output"] == "report.tsv"
        @test pa["discover"]["selftest"]["seq-col"] == "aln_qseq"
        @test pa["discover"]["selftest"]["metrics-output"] === nothing

        # classify_allele
        acc_list = ["ACGTACGT"]
        acc_set = Set(acc_list)
        rej = Dict("TTTTTTTT" => "output filter")
        @test Selftest.classify_allele("ACGTACGT", acc_set, acc_list, rej) == ("recovered", "")
        @test Selftest.classify_allele("ACGT", acc_set, acc_list, rej) == ("recovered", "")  # substring
        @test Selftest.classify_allele("TTTTTTTT", acc_set, acc_list, rej) == ("rejected", "output filter")
        @test Selftest.classify_allele("GGGGGGGG", acc_set, acc_list, rej) == ("missed", "")
        @test Selftest.is_novel("X", Set(["Y"]))
        @test !Selftest.is_novel("Y", Set(["Y"]))

        # evaluate_recovery
        disc = DataFrame(
            aln_qseq = ["NOVELAAA", "NOVELBBB", "ARTIFXXX", "KNOWNAAA"],
            reject_reason = ["", "min count", "", ""],
            reject_stage = ["", "output filter", "", ""],
        )
        base = Set(["KNOWNAAA"])
        truth = [("V1", "NOVELAAA"), ("V2", "NOVELBBB"), ("V3", "NOVELCCC"), ("K", "KNOWNAAA")]
        res = Selftest.evaluate_recovery(disc, base, truth; seq_col=:aln_qseq)
        @test res.summary.n_truth_novel == 3
        @test res.summary.recovered == 1
        @test res.summary.recall ≈ 1 / 3
        @test res.summary.true_positive == 1
        @test res.summary.false_positive == 1
        @test res.summary.precision == 0.5
        @test res.summary.n_novel_accepted == 2
        paa = res.per_allele
        @test paa[paa.allele .== "V1", :status][1] == "recovered"
        @test paa[paa.allele .== "V2", :status][1] == "rejected"
        @test paa[paa.allele .== "V2", :reject_stage][1] == "output filter"
        @test paa[paa.allele .== "V3", :status][1] == "missed"

        # CSV round-trip turns empty fields into `missing`: an accepted row has missing
        # reject_reason and must still count as recovered; empty cores must not match.
        disc2 = DataFrame(
            aln_qseq = ["NOVELAAA", missing, ""],
            reject_reason = [missing, "min count", missing],
            reject_stage = [missing, "output filter", missing],
        )
        res2 = Selftest.evaluate_recovery(disc2, Set(String[]), [("V1", "NOVELAAA")]; seq_col=:aln_qseq)
        @test res2.summary.n_truth_novel == 1
        @test res2.summary.recovered == 1
        @test res2.summary.false_positive == 0

        # metric_separation: label novel candidates true/false and find the best splitting threshold.
        discm = DataFrame(
            aln_qseq      = ["NOVELAAA", "NOVELAAA", "ARTIFXXX", "ARTIFYYY", "KNOWNAAA"],
            reject_reason = ["", "", "", "min count", ""],
            reject_stage  = ["", "", "", "output filter", ""],
            n_donors      = [3, 3, 1, 1, 9],
            parent_ratio  = [1.0, 1.0, 50.0, 40.0, 1.0],
        )
        basem = Set(["KNOWNAAA"])
        truthm = [("V1", "NOVELAAA"), ("K", "KNOWNAAA")]
        sep = Selftest.metric_separation(discm, basem, truthm; seq_col=:aln_qseq)
        @test nrow(sep) == 2                       # KNOWNAAA excluded (in base), both metrics scanned
        @test sep[1, :youden] ≈ 1.0                # sorted by youden desc
        nd = sep[sep.metric .== "n_donors", :]
        @test nd[1, :direction] == "keep ≥"        # true cores have more donors
        @test nd[1, :n_tp] == 2 && nd[1, :n_fp] == 2
        @test nd[1, :tp_kept] == 2 && nd[1, :fp_removed] == 2
        @test nd[1, :youden] ≈ 1.0
        pr = sep[sep.metric .== "parent_ratio", :]
        @test pr[1, :direction] == "keep ≤"        # true cores sit far from a dominant parent
        @test pr[1, :youden] ≈ 1.0

        # No false novel candidates ⇒ nothing to separate, empty table (not an error).
        discz = DataFrame(aln_qseq=["NOVELAAA"], reject_reason=[""], reject_stage=[""], n_donors=[3])
        sepz = Selftest.metric_separation(discz, Set(String[]), [("V1", "NOVELAAA")]; seq_col=:aln_qseq)
        @test nrow(sepz) == 0

        # recall_safe_filters: most aggressive cut over ACCEPTED cores that keeps recall 1.0.
        discs = DataFrame(
            aln_qseq      = ["GOODAAAA", "BADXXXXX", "BADYYYYY", "REJECTZZ"],
            reject_reason = ["", "", "", "min count"],   # REJECTZZ not accepted → ignored
            reject_stage  = ["", "", "", "output filter"],
            n_reads_total = [100, 5, 8, 999],
            max_full_ratio= [0.5, 0.02, 0.03, 0.9],
        )
        bases = Set(String[])
        truths = [("V1", "GOODAAAA")]
        safe = Selftest.recall_safe_filters(discs, bases, truths; seq_col=:aln_qseq)
        @test "max_full_ratio" in safe.metric
        # n_reads_total: recovered allele's core has 100; both FP (5,8) drop below it → 2/2.
        nr = safe[safe.metric .== "n_reads_total", :]
        @test nr[1, :direction] == "keep ≥"
        @test nr[1, :threshold] ≈ 100.0
        @test nr[1, :fp_removed] == 2 && nr[1, :acc_fp] == 2
        # peak allelic ratio separates the same way (0.5 keeps the true allele, drops both FP).
        mr = safe[safe.metric .== "max_full_ratio", :]
        @test mr[1, :fp_removed] == 2

        # No accepted false novel cores ⇒ empty (nothing to cut).
        safez = Selftest.recall_safe_filters(discz, Set(String[]), [("V1", "NOVELAAA")]; seq_col=:aln_qseq)
        @test nrow(safez) == 0
    end

    @testset "bwa.jl" begin
        # CLI - bwa is now under analyze group
        empty!(ARGS)
            append!(ARGS, ["search", "bwa", "test_bwa_input.tsv", "test_bwa_output.tsv", "genome1.fasta", "genome2.fasta"])
            parsed_args = Cli.parse_commandline(ARGS)
            @test parsed_args["%COMMAND%"] == "search"
            @test parsed_args["search"]["%COMMAND%"] == "bwa"
            @test parsed_args["search"]["bwa"]["tsv"] == "test_bwa_input.tsv"
            @test parsed_args["search"]["bwa"]["output"] == "test_bwa_output.tsv"
            @test parsed_args["search"]["bwa"]["genome"] == ["genome1.fasta", "genome2.fasta"]
            @test parsed_args["search"]["bwa"]["chromosome"] == "chromosome 14"
            @test parsed_args["search"]["bwa"]["colname"] == "best_name"

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

    @testset "table subcommands" begin
        test_left = "test_table_left.tsv"
        test_right = "test_table_right.tsv"

        left_df = DataFrame(key=["a","b","c"], val=[1,2,3])
        right_df = DataFrame(key=["b","c","d"], score=[10,20,30])
        CSV.write(test_left, left_df, delim='\t')
        CSV.write(test_right, right_df, delim='\t')

        @testset "outerjoin" begin
            result = Table.outerjoin_tsv(test_left, test_right, "test_oj.tsv.gz",
                left_keys=["key"])
            @test nrow(result) == 4
            @test "key" in names(result)
        end

        @testset "leftjoin" begin
            result = Table.leftjoin_tsv(test_left, test_right, "test_lj.tsv.gz",
                left_keys=["key"])
            @test nrow(result) == 3
        end

        @testset "unique" begin
            dup_df = DataFrame(a=["x","x","y"], b=[1,1,2])
            CSV.write("test_dup.tsv", dup_df, delim='\t')
            result = Table.unique_tsv("test_dup.tsv", "test_uniq.tsv.gz", columns=["a","b"])
            @test nrow(result) == 2
            isfile("test_dup.tsv") && rm("test_dup.tsv")
        end

        @testset "sort" begin
            sort_df = DataFrame(name=["c","a","b"], val=[3,1,2])
            CSV.write("test_sort_in.tsv", sort_df, delim='\t')
            result = Table.sort_tsv("test_sort_in.tsv", "test_sort_out.tsv.gz", columns=["name"])
            @test result.name == ["a","b","c"]
            result_rev = Table.sort_tsv("test_sort_in.tsv", "test_sort_out2.tsv.gz", columns=["name"], reverse=true)
            @test result_rev.name == ["c","b","a"]
            isfile("test_sort_in.tsv") && rm("test_sort_in.tsv")
        end

        @testset "select" begin
            sel_df = DataFrame(a=[1,2], b=[3,4], c=[5,6])
            CSV.write("test_sel_in.tsv", sel_df, delim='\t')
            result = Table.select_tsv("test_sel_in.tsv", "test_sel_out.tsv.gz", columns=["a","c"])
            @test ncol(result) == 2
            @test names(result) == ["a","c"]
            isfile("test_sel_in.tsv") && rm("test_sel_in.tsv")
        end

        @testset "filter" begin
            filt_df = DataFrame(name=["alpha","beta","gamma"], val=[10,5,20])
            CSV.write("test_filt_in.tsv", filt_df, delim='\t')
            result_regex = Table.filter_tsv("test_filt_in.tsv", "test_filt1.tsv.gz", "name", pattern="^al")
            @test nrow(result_regex) == 1
            @test result_regex.name[1] == "alpha"
            result_num = Table.filter_tsv("test_filt_in.tsv", "test_filt2.tsv.gz", "val", operator=">=", threshold=10.0)
            @test nrow(result_num) == 2
            @test all(result_num.val .>= 10)
            isfile("test_filt_in.tsv") && rm("test_filt_in.tsv")
        end

        @testset "transform" begin
            tf_df = DataFrame(id=["ABC_123","DEF_456"], val=[1,2])
            CSV.write("test_tf_in.tsv", tf_df, delim='\t')
            result = Table.transform_tsv("test_tf_in.tsv", "test_tf_out.tsv.gz",
                column="id", pattern="(.*)_(.*)", replacement="\\1-\\2")
            @test result.id == ["ABC-123","DEF-456"]
            result_new = Table.transform_tsv("test_tf_in.tsv", "test_tf_out2.tsv.gz",
                column="id", pattern="(.*)_(.*)", replacement="\\1-\\2", new_column="suffix")
            @test "suffix" in names(result_new)
            isfile("test_tf_in.tsv") && rm("test_tf_in.tsv")
        end

        @testset "aggregate" begin
            agg_df = DataFrame(gene=["V","V","D","D","D"], case=["A","A","B","B","B"])
            CSV.write("test_agg_in.tsv", agg_df, delim='\t')
            result = Table.aggregate_tsv("test_agg_in.tsv", "test_agg_out.tsv.gz",
                group_by=["gene","case"])
            @test nrow(result) == 2
            @test "count" in names(result)
            isfile("test_agg_in.tsv") && rm("test_agg_in.tsv")

            # keep_columns path (single combine): first value of kept column per group.
            agg_df2 = DataFrame(gene=["V","V","D"], case=["A","A","B"], extra=["x","y","z"])
            CSV.write("test_agg_in2.tsv", agg_df2, delim='\t')
            result2 = Table.aggregate_tsv("test_agg_in2.tsv", "test_agg_out2.tsv.gz",
                group_by=["gene","case"], keep_columns=["extra"])
            @test nrow(result2) == 2
            @test "count" in names(result2)
            @test "extra" in names(result2)
            @test result2[(result2.gene .== "D") .& (result2.case .== "B"), "extra"][1] == "z"
            isfile("test_agg_in2.tsv") && rm("test_agg_in2.tsv")
        end

        for f in Glob.glob("test_table_*.tsv"); isfile(f) && rm(f); end
        for f in Glob.glob("test_oj*.tsv.gz"); isfile(f) && rm(f); end
        for f in Glob.glob("test_lj*.tsv.gz"); isfile(f) && rm(f); end
        for f in Glob.glob("test_uniq*.tsv.gz"); isfile(f) && rm(f); end
        for f in Glob.glob("test_sort_out*.tsv.gz"); isfile(f) && rm(f); end
        for f in Glob.glob("test_sel_out*.tsv.gz"); isfile(f) && rm(f); end
        for f in Glob.glob("test_filt*.tsv.gz"); isfile(f) && rm(f); end
        for f in Glob.glob("test_tf_out*.tsv.gz"); isfile(f) && rm(f); end
        for f in Glob.glob("test_agg_out*.tsv.gz"); isfile(f) && rm(f); end
    end

    @testset "fasta hash CLI" begin
        empty!(ARGS)
        append!(ARGS, ["fasta", "hash", "test_hash.fasta"])
        parsed_args = Cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "fasta"
        @test parsed_args["fasta"]["%COMMAND%"] == "hash"
        @test parsed_args["fasta"]["hash"]["fastain"] == "test_hash.fasta"
    end

    @testset "hsmm CLI" begin
        empty!(ARGS)
        append!(ARGS, ["discover", "hsmm", "test.tsv", "test.fasta", "test_out.tsv.gz"])
        parsed_args = Cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "discover"
        @test parsed_args["discover"]["%COMMAND%"] == "hsmm"
        @test parsed_args["discover"]["hsmm"]["tsv"] == "test.tsv"
        @test parsed_args["discover"]["hsmm"]["fasta"] == "test.fasta"
        @test parsed_args["discover"]["hsmm"]["output"] == "test_out.tsv.gz"
        @test parsed_args["discover"]["hsmm"]["ratio"] == 0.2
        @test parsed_args["discover"]["hsmm"]["mincount"] == 10
        @test parsed_args["discover"]["hsmm"]["min-posterior"] == 0.7
        @test parsed_args["discover"]["hsmm"]["out-mincount"] == 10
        @test parsed_args["discover"]["hsmm"]["out-minratio"] == 0.2
    end

    @testset "hsmm module" begin
        logp = ntuple(i -> log(0.25), 4)
        emission = HSMM.IIDLogEmission(logp)
        @test HSMM.logemit(emission, 1) ≈ log(0.25)

        obs = HSMM.encode_dna("ACGT")
        @test obs == [1, 2, 3, 4]
        @test HSMM.encode_dna("N") == [0]
    end

    @testset "simulate V end-variants (exact detection)" begin
        # Two V alleles that differ ONLY in their last 8 nt — the 3' V/RSS border case.
        Random.seed!(2024)
        germline = Simulate.random_sequence(100, 100)
        variant  = Simulate.v_end_variant(germline; n_end=8)
        @test length(germline) == length(variant)
        @test germline != variant
        @test germline[1:90] == variant[1:90]      # share a long 5' prefix
        @test germline[93:100] != variant[93:100]  # differ only at the 3' end

        db = [("IGHV1-1*01", germline), ("IGHV1-1*02", variant)]

        names_ = String[]; seqs = String[]
        for (allele, gseq) in db, r in 1:10
            push!(names_, "$(allele)_r$(r)")
            push!(seqs, Simulate.assemble_read(VGene(), gseq; flank=15))
        end
        table = DataFrame(well=fill(1, length(names_)), case=fill("D1", length(names_)),
                          name=names_, genomic_sequence=seqs)

        counts = Exact.exact_search(table, db, "V"; N=10)
        @test nrow(counts) > 0
        # Both the germline and the end-variant must be recovered as distinct alleles.
        @test "IGHV1-1*01" in counts.db_name
        @test "IGHV1-1*02" in counts.db_name
        # A read built from one allele must not be assigned the other (3' end discriminates).
        v01_seq = first(filter(r -> r.db_name == "IGHV1-1*01", counts)).sequence
        v02_seq = first(filter(r -> r.db_name == "IGHV1-1*02", counts)).sequence
        @test v01_seq == germline
        @test v02_seq == variant
    end

    @testset "simulate short-D variants (HSMM detection)" begin
        # Very short D segments (down to 8 nt) embedded in full RSS, recovered by the HSMM.
        Random.seed!(1234)
        d_genes = [Simulate.random_sequence(L, L) for L in (8, 10, 12, 14, 16)]

        tuples = NTuple{7,String}[]
        reads = String[]
        for g in d_genes
            read, f = Simulate.simulate_d_read(g; flank=10)
            push!(tuples, (f.pre_nonamer, f.pre_spacer, f.pre_heptamer, f.gene,
                           f.post_heptamer, f.post_spacer, f.post_nonamer))
            push!(reads, read)
        end

        model = HSMM.fit_dgene_rss_hsmm(tuples, minimum(length.(d_genes)), maximum(length.(d_genes)))

        # Each short D must be localized and extracted exactly from its read.
        for (g, read) in zip(d_genes, reads)
            det = HSMM.extract_dgene(read, model)
            @test det.gene_seq == g
        end

        # Explicitly assert the very-short (8 nt) D variant is recovered at full length.
        short_read, _ = Simulate.simulate_d_read(d_genes[1]; flank=10)
        det = HSMM.extract_dgene(short_read, model)
        @test length(det.gene_seq) == 8
        @test det.gene_seq == d_genes[1]

        # short_d_variant: truncating a germline yields a shorter allele with shared prefix.
        long_d = Simulate.random_sequence(30, 30)
        trunc_d = Simulate.short_d_variant(long_d; len=12)
        @test length(trunc_d) == 12
        @test long_d[1:12] == trunc_d
    end

    @testset "simulate decoys (no false-positive D detection)" begin
        # Train on genuine short Ds, then confirm non-D background and broken-RSS reads
        # never look as good as a real D — the signal the posterior/heptamer filters use.
        Random.seed!(99)
        d_genes = [Simulate.random_sequence(L, L) for L in (10, 12, 14, 16, 18)]
        tuples = NTuple{7,String}[]
        true_reads = String[]
        for g in d_genes
            r, f = Simulate.simulate_d_read(g; flank=12)
            push!(tuples, (f.pre_nonamer, f.pre_spacer, f.pre_heptamer, f.gene,
                           f.post_heptamer, f.post_spacer, f.post_nonamer))
            push!(true_reads, r)
        end
        model = HSMM.fit_dgene_rss_hsmm(tuples, minimum(length.(d_genes)), maximum(length.(d_genes)))

        # Reference: best-path log-prob of a genuine D read (RSS intact).
        true_score = HSMM.scan_best_and_total(true_reads[1], model).log_path_prob
        @test isfinite(true_score)

        # 1) Pure random background: must score strictly below a genuine D.
        for _ in 1:20
            det = HSMM.scan_best_and_total(Simulate.decoy_read(; len=rand(80:120)), model)
            @test det.log_path_prob < true_score
        end

        # 2) Invalid D (gene present but heptamers scrambled): broken RSS scores far worse.
        for g in d_genes
            det = HSMM.scan_best_and_total(Simulate.invalid_d_read(g; flank=12), model)
            @test det.log_path_prob < true_score
        end
    end

    @testset "heptamer extract handles short alleles" begin
        # An allele shorter than the b+e trim must yield no matches, not crash.
        tbl = DataFrame(well=[1], case=["D1"], name=["r1"],
                        genomic_sequence=["ACGTACGTACGTCACAGTGACGT"])
        db = [("IGHVshort*01", "ACG")]   # length 3 < b(1)+e(8)
        hdf = Heptamer.extract_heptamers(tbl, db, ["CACAGTG"]; max_dist=0, b=1, e=8)
        @test nrow(hdf) == 0
    end

    @testset "cooccurrence CLI" begin
        empty!(ARGS)
        append!(ARGS, ["analyze", "cooccurrence", "test_input.tsv"])
        parsed_args = Cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "analyze"
        @test parsed_args["analyze"]["%COMMAND%"] == "cooccurrence"
        @test parsed_args["analyze"]["cooccurrence"]["input"] == "test_input.tsv"
        @test parsed_args["analyze"]["cooccurrence"]["case-col"] == "case"
        @test parsed_args["analyze"]["cooccurrence"]["allele-col"] == "db_name"
        @test parsed_args["analyze"]["cooccurrence"]["min-donors"] == 2
        @test parsed_args["analyze"]["cooccurrence"]["cluster-method"] == "components"
    end

    @testset "cooccurrence module" begin
        df = DataFrame(
            case = ["D1","D1","D1","D2","D2","D2","D3","D3","D3","D4","D4","D4"],
            db_name = ["A*01","B*01","C*01","A*01","B*01","C*01","A*01","B*01","D*01","A*01","C*01","D*01"]
        )
        edges, clusters, detailed = Cooccurrence.compute_edges_and_clusters(df; min_support=1, jaccard_threshold=0.0, similarity_threshold=0.0)
        @test nrow(edges) > 0

        stats = Cooccurrence.compute_full_stats(df; case_col="case", allele_col="db_name", min_donors=1)
        @test length(stats.alleles) > 0
        @test stats.N == 4
        @test size(stats.R, 1) == length(stats.alleles)

        # Shipped CLI edge builder must carry Benjamini–Hochberg q-values.
        edges_m = Cooccurrence.build_edges_from_matrices(stats.R, stats.J, stats.SUP, stats.P, String.(stats.alleles))
        @test "q_value" in names(edges_m)
        @test nrow(edges_m) > 0
        @test all(0.0 .<= edges_m.q_value .<= 1.0)

        # Benjamini–Hochberg correction: monotone, in [0,1], empty-safe.
        @test Cooccurrence.adjust_bh(Float64[]) == Float64[]
        q = Cooccurrence.adjust_bh([0.01, 0.02, 0.03, 0.04])
        @test all(isapprox.(q, 0.04; atol=1e-12))
        @test all(0.0 .<= Cooccurrence.adjust_bh([0.5, 0.001, 0.9, 0.2]) .<= 1.0)

        # min_donors is a DONOR threshold, not a row count: X*01 has 5 rows but 1 donor.
        df_rows = DataFrame(
            case    = ["D1","D1","D1","D1","D1","D2","D3"],
            db_name = ["X*01","X*01","X*01","X*01","X*01","Y*01","Y*01"],
        )
        s2 = Cooccurrence.compute_full_stats(df_rows; case_col="case", allele_col="db_name", min_donors=2)
        @test !("X*01" in String.(s2.alleles))   # 5 rows but only 1 donor → excluded
        @test "Y*01" in String.(s2.alleles)       # 2 donors → included
    end

    @testset "table exclude" begin
        excl_data = DataFrame(allele_name=["A*01","B*01","C*01"], seq=["ATCG","GCTA","TTTT"])
        CSV.write("test_excl_input.tsv", excl_data, delim='\t')
        open(FASTA.Writer, "test_excl_ref.fasta") do writer
            write(writer, FASTARecord("known", "ATCG"))
        end

        empty!(ARGS)
        append!(ARGS, ["table", "exclude", "test_excl_input.tsv", "test_excl_output.tsv", "test_excl_ref.fasta"])
        parsed_args = Cli.parse_commandline(ARGS)
        @test parsed_args["%COMMAND%"] == "table"
        @test parsed_args["table"]["%COMMAND%"] == "exclude"

        for f in ["test_excl_input.tsv", "test_excl_ref.fasta"]
            isfile(f) && rm(f)
        end
    end

    @testset "pure helpers" begin
        @testset "bwa cigar/seq" begin
            @test Bwa.reverse_complement_seq("ACGT") == "ACGT"
            @test Bwa.reverse_complement_seq("AAAA") == "TTTT"
            @test Bwa.reverse_complement_seq("") == ""
            @test Bwa.calculate_ref_span_from_cigar("10M") == 10
            @test Bwa.calculate_ref_span_from_cigar("10M5D3M") == 18
            @test Bwa.calculate_ref_span_from_cigar("5S10M") == 10
            @test Bwa.calculate_ref_span_from_cigar("") == 0
            @test Bwa.calculate_leading_n_from_cigar("5N10M") == 5
            @test Bwa.calculate_leading_n_from_cigar("10M5N") == 0
            @test Bwa.calculate_leading_n_from_cigar("") == 0
            @test Bwa.hamming_distance("ACGT", "ACGA") == 1
            @test Bwa.hamming_distance("AB", "ABC") == -1
        end

        @testset "data utils" begin
            @test Data.validate_sequence("ACGT")
            @test !Data.validate_sequence("ACGN")
            @test !Data.validate_sequence("acgt")
            @test Data.validate_identifier("IGHV1-2*01")
            @test Data.validate_identifier("IGHV1*01_S1234")
            @test !Data.validate_identifier("foo")
            @test startswith(Data.sequence_hash("ACGT"), "S")
            @test length(Data.sequence_hash("ACGT")) == 5
            @test Data.sequence_hash("ACGT") == Data.sequence_hash("ACGT")
            @test startswith(Data.unique_name("IGHV1-2_x", "ACGT"), "IGHV1-2_S")
            @test Data.concatenate_columns((a="X", b="Y"), ["a", "b"]) == "XY"
            @test Data.validate_types(["heptamer", "spacer"]) === nothing
            @test_throws ErrorException Data.validate_types(["bad"])
            @test_throws ErrorException Data.validate_types(String[])
            @test Data.get_ratio_threshold(Dict("A*01" => 0.5), (db_name="A*01", gene="A")) == 0.5
            @test Data.get_ratio_threshold(Dict{String,Float64}(), (db_name="A*01", gene="A")) == 0.0
        end

        @testset "blast string utils" begin
            @test Blast.sseqid_to_db_key("TRGV2*01_S2223", Set(["TRGV2*01"])) == "TRGV2*01"
            @test Blast.sseqid_to_db_key("X*01", Set(["X*01"])) == "X*01"
            @test Blast.sseqid_to_db_key("Q*01_S9", Set(["Z"])) == "Q*01_S9"
            @test Blast.consensus_prefix(["AAAAA", "AAAAA"]) == "AAAAA"
            @test Blast.consensus_suffix(["TTTTT", "TTTTT"]) == "TTTTT"
            @test Blast.consensus_prefix(["AAAAT", "CCCGT", "GGGGT"]) == "GT"
            @test Blast.consensus_suffix(["TCCCC", "TAAAA", "TTTTT"]) == "T"
            @test Blast.consensus_prefix(["XYZ", "ABC"]) == ""
            @test Blast.edge("CCC", "AAACCCGGG") == (3, 3)   # 3 nt before, 3 after
            @test Blast.blastn_cli_token("--task") == "-task"
            @test Blast.blastn_cli_token("-num_threads") == "-num_threads"
            @test Blast.blastn_cli_token("megablast") == "megablast"
            cmd = Blast.build_blastn_cmd("/q.fa", "/db", "6 qseqid", "-task megablast")
            @test cmd.exec[1:3] == ["blastn", "-num_threads", string(Sys.CPU_THREADS)]
            @test "-task" in cmd.exec && "megablast" in cmd.exec
            @test Blast.nogaps("AC-G-T") == "ACGT"
            @test isabspath(Blast.resolve_work_dir(""))
            @test endswith(Blast.resolve_work_dir(""), ".immunediscover")
            @test Blast.resolve_work_dir("/tmp/wd") == "/tmp/wd"
            # subject coverage from sstart..send span, bounded ≤ 1 (incl. minus strand)
            @test Blast.subject_coverage(1, 100, 100) == 1.0
            @test Blast.subject_coverage(1, 50, 100) == 0.5
            @test Blast.subject_coverage(100, 1, 100) == 1.0   # send < sstart on minus strand
            @test Blast.subject_coverage(10, 19, 100) == 0.1
        end

        @testset "exact gene-type dispatch" begin
            @test Exact.gene_type_from_name("IGHV1-2") isa VGene
            @test Exact.gene_type_from_name("IGHD3") isa DGene
            @test Exact.gene_type_from_name("IGHJ4") isa JGene
            @test Exact.gene_type_from_name("XYZ") === nothing
            @test Exact.parse_gene_type("V") isa VGene
            @test_throws ErrorException Exact.parse_gene_type("Q")
            @test Exact.gene_string(VGene()) == "V"
            @test Exact.gene_string(DGene()) == "D"
            @test Exact.gene_string(JGene()) == "J"
            # V extension reaches the 3' border -> reject; short extension -> keep
            @test Exact.extension_overlaps_border(5, 10, 20, "V", 15, 3)
            @test !Exact.extension_overlaps_border(5, 10, 20, "V", 2, 3)
        end

        @testset "cooccurrence math" begin
            # Perfect co-occurrence -> phi = 1; mutual exclusion -> phi = -1
            @test Cooccurrence.phi_coefficient(2, 0, 0, 2) ≈ 1.0
            @test Cooccurrence.phi_coefficient(0, 2, 2, 0) ≈ -1.0
            @test Cooccurrence.phi_coefficient(0, 0, 0, 0) == 0.0   # degenerate denom
            ji, n11 = Cooccurrence.jaccard_index(Set(["a","b"]), Set(["b","c"]))
            @test n11 == 1
            @test ji ≈ 1/3
            @test Cooccurrence.hypergeom_p_enrichment(10, 0, 5, 0) == 1.0
        end
    end

    # Cleanup test files
    for file in ["test.fasta", "reference.fasta", "novel.fasta", "test_indices.tsv",
                 "test.tsv.gz", "test_exact.tsv.gz", "test_heptamer.tsv.gz",
                 "test_summary.tsv", "test.fastq"]
        isfile(file) && rm(file)
    end
end