module immunediscover
    include("cli.jl")
    include("demultiplex.jl")
    include("simulate.jl")
    include("data.jl")
    include("profile.jl")
    include("trim.jl")
    include("exact.jl")
    include("heptamer.jl")
    include("hamming.jl")
    include("patterns.jl")
    include("regex.jl")
    include("bwa.jl")
    include("nwpattern.jl")
    include("blast.jl")
    include("keyedsets.jl")
    include("linkage.jl")

    using .cli
    using .demultiplex
    using .simulate
    using .data
    using .profile
    using .trim
    using .exact
    using .heptamer
    using .hamming
    using .nwpattern
    using .patterns
    using .regex
    using .bwa
    using .blast
    using .keyedsets
    using .linkage

    using CSV
    using DataFrames
    using UnicodePlots
    using Glob
    using Statistics
    using DataStructures

    export load_fasta, blast_discover

    const TRUST_MINCOUNT = 5

    mutable struct AlignmentStats
        total_attempts::Int
        prefix_failures::Int
        suffix_failures::Int
        successful_trims::DefaultDict{String, Vector{String}}
        failed_genes::DefaultDict{String, Vector{String}}
    end

    AlignmentStats() = AlignmentStats(0, 0, 0, DefaultDict{String, Vector{String}}(Vector{String}), DefaultDict{String, Vector{String}}(Vector{String}))

    """
        concatenate_columns(row, col_names)

    Helper function to concatenate content from columns
    """
    function concatenate_columns(row, col_names)
        concatenated_string = ""
        for col_name in col_names
            concatenated_string *= getproperty(row, Symbol(col_name))
        end
        concatenated_string
    end

    function validate_types(types)
        allowed_types = ["heptamer", "spacer", "nonamer"]
        for type in types
            if !(type in allowed_types)
                error("Invalid type: $type. Allowed types are: heptamer, spacer, nonamer.")
            end
        end

        if isempty(types)
            error("At least one type must be specified.")
        end
    end

    # Function to get threshold for new gene
    function get_ratio_threshold(expect_dict, row; type="allele_ratio")
        if row.db_name in keys(expect_dict)
            @info "Applying $type >= $(expect_dict[row.db_name]) for $(row.db_name)"
            return expect_dict[row.db_name]
        elseif row.gene in keys(expect_dict)
            @info "Applying $type >= $(expect_dict[row.gene]) for $(row.db_name)"
            return expect_dict[row.gene]
        else
            return 0  # Return ratio if not in expect_dict
        end
    end

    """
        real_main(args=[])

    Main function for Julia
    """
    function real_main(args=[])
        parsed_args = parse_commandline(args)
        if parsed_args !== nothing
            if get(parsed_args,"%COMMAND%","") == "demultiplex"
                @info "Demultiplexing"
                table, stats = immunediscover.demultiplex.demux(parsed_args["demultiplex"]["fastq"],
                                                         parsed_args["demultiplex"]["indices"],
                                                         parsed_args["demultiplex"]["forwardarrayindex"],
                                                         min_length=parsed_args["demultiplex"]["length"],
                                                         save_fastq_files=parsed_args["demultiplex"]["split"])
                
                # Apply case filtering if regex is provided
                case_filter_regex = get(parsed_args["demultiplex"], "case-filter-regex", nothing)
                if case_filter_regex !== nothing
                    @info "Filtering cases with regex: $case_filter_regex"
                    original_count = nrow(table)
                    filter!(x -> startswith(x.case, Regex(case_filter_regex)), table)
                    filtered_count = nrow(table)
                    @info "Filtered from $original_count to $filtered_count rows ($(original_count - filtered_count) rows removed)"
                end
                
                logfile = "$(parsed_args["demultiplex"]["output"]).log"
                CSV.write(logfile, stats, delim='\t')
                count_df = combine(groupby(table, :case), :case => length => :count)
                sort!(count_df, :count, rev=true)
                @info "Demultiplexing statistics"
                println(barplot(count_df.case, count_df.count))
                output = cli.always_gz(parsed_args["demultiplex"]["output"])
                CSV.write(output, table, compress=true, delim='\t')
                @info "Demultiplexing data saved in compressed $output file"
                @info "Demultiplexing detailed statistics saved in $logfile file"
            end

            if get(parsed_args,"%COMMAND%","") == "heptamer"
                @info "Extracting heptamers"
                heptamers = heptamer.load_heptamers(parsed_args["heptamer"]["json"])
                chain = parsed_args["heptamer"]["chain"]
                @info "Parsing for $(chain)"
                @info "Using heptamers $(join(heptamers[chain],','))"
                db = load_fasta(parsed_args["heptamer"]["fasta"])
                @info "Loaded $(length(db)) query sequences"
                table = load_demultiplex(parsed_args["heptamer"]["tsv"])
                data = heptamer.extract_heptamers(table, db, heptamers[chain]; max_dist=parsed_args["heptamer"]["maxdist"], b=parsed_args["heptamer"]["begin"]+1, e=parsed_args["heptamer"]["end"])
                output = cli.always_gz(parsed_args["heptamer"]["output"])
                CSV.write(output, data, compress=true, delim='\t')
                @info "Heptamer data saved in compressed $output file"
                @info "Summarizing"
                data_df = CSV.File(output, delim='\t') |> DataFrame
                unique = heptamer.summarize(data_df, db, ratio=parsed_args["heptamer"]["ratio"],count=parsed_args["heptamer"]["mincount"])
                CSV.write(parsed_args["heptamer"]["summary"], unique, delim='\t')
                @info "Summary heptamer data saved in compressed $(parsed_args["heptamer"]["summary"]) file"
            end

            if get(parsed_args,"%COMMAND%","") == "diff"
                fasta = parsed_args["diff"]["fasta"]
                fasta_files = [(file=file, records=immunediscover.load_fasta.(file)) for file in fasta]
                sets = [(file=x, set=KeyedSet(reverse.(y))) for (x,y) in fasta_files]
                for i in eachindex(sets)
                    for j in i+1:length(sets)
                        # All vs vall comparisons
                        println("Comparing $(sets[i].file) vs $(sets[j].file)")
                        println("Union: $(length(union(sets[i].set, sets[j].set)))")
                        println("Intersection: $(length(intersect(sets[i].set, sets[j].set)))")
                        ab = last.(collect(setdiff(sets[i].set, sets[j].set)))
                        ba = last.(collect(setdiff(sets[j].set, sets[i].set)))
                        println("Difference ($(length(ab))): \n$(join(ab, "\n"))")
                        println("Difference ($(length(ba))): \n$(join(ba, "\n"))")
                    end
                end
            end
            if get(parsed_args,"%COMMAND%","") == "trim"
                @info "Trimming"
                table = load_demultiplex(parsed_args["trim"]["input"])
                db = [r[2] for r in load_fasta(parsed_args["trim"]["fasta"])]
                weights = parsed_args["trim"]["weights"]
                pos = parsed_args["trim"]["position"]
                mlen = parsed_args["trim"]["length"]

                prof_start, prof_stop = trim.trim_profiles(db, weights)
                ok, fragments, region = trim.find_segments(table.genomic_sequence, prof_start, prof_stop, mlen)

                subtable = table[ok,:]
                if pos
                    start = [r[1] for r in region]
                    stop = [r[2] for r in region]
                    subtable[:,:start] = start
                    subtable[:,:stop] = stop
                else
                    subtable[:,:trimmed_sequence] = fragments
                end
                output = cli.always_gz(parsed_args["trim"]["output"])
                CSV.write(output, subtable, compress=true, delim='\t')
                @info "Trimmed data saved in compressed $output file"
            end

            if get(parsed_args, "%COMMAND%", "") == "blast"
                parsed_args = apply_blast_presets!(parsed_args)
                gene = parsed_args["blast"]["gene"]
                if haskey(BLAST_PRESETS, gene)
                    @info "Initial $gene preset parameters (can be overridden by command line arguments):"
                    for (param, value) in BLAST_PRESETS[gene]
                        @info "  --$param = $value"
                    end
                end
                @info "Running with the following parameters:"
                for (param, value) in parsed_args["blast"]
                    @info "  --$param = $value"
                end
                @info "Discovery with BLAST assignments"

                # Load paths and parameters
                fasta_path = parsed_args["blast"]["fasta"]
                affixes_path = fasta_path * ".affixes"
                DB = load_fasta(fasta_path, validate=false)

                # Parameters
                verbose = parsed_args["blast"]["verbose"]
                overwrite = parsed_args["blast"]["overwrite"]
                minquality = parsed_args["blast"]["minquality"]
                forward_extension = parsed_args["blast"]["forward"]
                reverse_extension = parsed_args["blast"]["reverse"]
                keep_failed = parsed_args["blast"]["keep-failed"]
                min_corecov = parsed_args["blast"]["min-corecov"]

                # Arbitrary choice what is short, but it's just a warning
                if (forward_extension < 7) && (forward_extension > 0)
                    @warn "Forward extension $forward_extension is short and may lead to false positives due to problem with ambigous trimming alignment - you've been warned!"
                end
                if (reverse_extension < 7) && (reverse_extension > 0)
                    @warn "Reverse extension $reverse_extension is short and may lead to false positives due to problem with ambigous trimming alignment - you've been warned!"
                end

                # Process pseudo genes if provided
                db_p = Vector{Tuple{String, String}}()
                pseudo = parsed_args["blast"]["pseudo"]
                if !isempty(pseudo)
                    for (name, seq) in blast.read_fasta(pseudo)
                        push!(db_p, ("P"*name, seq))
                    end
                end

                # Add regular genes
                for (name, seq) in blast.read_fasta(fasta_path)
                    push!(db_p, (name, seq))
                end

                # Save combined database
                combined_fasta_path = blast.replace_extension(fasta_path, "fasta", tag="-combined")
                blast.save_to_fasta(db_p, combined_fasta_path)

                # Handle sequence extension based on extension parameters
                if forward_extension == 0 && reverse_extension == 0
                    # If no extension is needed, use the combined fasta directly
                    ext_fasta_path = combined_fasta_path
                    @info "No sequence extension requested, using original sequences"

                    # Create empty affixes file for consistency
                    empty_affixes = [(name="", prefix="", suffix="")]
                    CSV.write(affixes_path, DataFrame(empty_affixes), delim='\t')
                else
                    # Extend gene sequences if either extension is non-zero
                    ext_fasta_path = blast.replace_extension(combined_fasta_path, "fasta", tag="-extended")
                    if !isfile(ext_fasta_path) || overwrite
                        @info "Extending gene sequences by $forward_extension forward and $reverse_extension reverse nucleotides"
                        demux = blast.load_csv(parsed_args["blast"]["input"])

                        extended = accumulate_affixes(DB, demux,
                            forward_extension=forward_extension,
                            reverse_extension=reverse_extension)

                        affixes = save_extended(extended, ext_fasta_path)
                        CSV.write(affixes_path, DataFrame(affixes, [:name, :prefix, :suffix]), delim='\t')
                        @info "Saved affixes in $affixes_path"
                    else
                        @info "Using existing extended sequences from $ext_fasta_path"
                    end
                end

                # Run BLAST discovery
                blast_clusters = blast_discover(
                    parsed_args["blast"]["input"],
                    ext_fasta_path,
                    max_dist=parsed_args["blast"]["maxdist"],
                    min_edge=parsed_args["blast"]["edge"],
                    min_scov=parsed_args["blast"]["subjectcov"],
                    args=parsed_args["blast"]["args"],
                    verbose=verbose,
                    overwrite=overwrite
                )

                # Process affixes only if extensions were used
                if forward_extension != 0 || reverse_extension != 0
                    if isfile(affixes_path)
                        @info "Loading affixes from $affixes_path"
                        affixes = Tuple.(collect.(eachrow(CSV.read(affixes_path, DataFrame, delim='\t'))))

                        if !any(occursin.("prefix", names(blast_clusters)))
                            @info "Using affixes to trim extended sequences"
                            affixes_df = DataFrame(affixes, [:sseqid, :prefix, :suffix])
                            leftjoin!(blast_clusters, affixes_df, on=:sseqid)

                            # Replace missing values with empty strings
                            blast_clusters.prefix = coalesce.(blast_clusters.prefix, "")
                            blast_clusters.suffix = coalesce.(blast_clusters.suffix, "")
                        end
                    end
                else
                    # Add empty prefix/suffix columns when no extension was used
                    if !any(occursin.("prefix", names(blast_clusters)))
                        blast_clusters.prefix .= ""
                        blast_clusters.suffix .= ""
                    end
                end

                # Re-align and process gene coress
                @info "Re-aligning gene cores and processing alleles"
                DBdict = Dict(DB)
                transform!(blast_clusters, :sseqid => ByRow(x -> DBdict[coalesce(String(x), "")]) => :db_seq)

                stats = AlignmentStats()
                transform!(blast_clusters,
                    [:qseq, :prefix, :suffix, :db_seq, :sseqid] =>
                    ByRow((qseq, prefix, suffix, db_seq, sseqid) ->
                        blast.trim_and_align_sequence(
                            coalesce(String(qseq), ""),
                            coalesce(String(prefix), ""),
                            coalesce(String(suffix), ""),
                            coalesce(String(db_seq), ""),
                            stats,
                            min_quality=minquality,
                            sseqid=coalesce(String(sseqid), "")
                        )
                    ) => [:aln_qseq, :aln_mismatch]
                )

                # Drop failed/empty alignments by default (prevents "nothing" strings)
                if !keep_failed
                    before = nrow(blast_clusters)
                    filter!(x -> !ismissing(x.aln_qseq) && (x.aln_qseq != "") && (x.aln_qseq != "nothing") && x.aln_mismatch >= 0, blast_clusters)
                    @info "Dropped $(before - nrow(blast_clusters)) rows with failed/empty alignments"
                end

                # Core coverage filter to suppress short local matches
                transform!(blast_clusters, [:aln_qseq, :db_seq] => ByRow((aq, ds) -> begin
                    aq_s = coalesce(String(aq), "")
                    ds_s = coalesce(String(ds), "")
                    if isempty(aq_s) || isempty(ds_s)
                        0.0
                    else
                        length(aq_s) / length(ds_s)
                    end
                end) => :corecov)
                filter!(x -> x.corecov >= min_corecov, blast_clusters)
                # Info
                for (gene, failures) in sort(collect(stats.failed_genes))
                    if !isempty(failures)
                        failure_counts = counter(failures)
                        for (reason, count) in sort(collect(failure_counts), by=x->x[2], rev=true)
                            trimmed_count = length(stats.successful_trims[gene])
                            @warn "Cluster for $gene: $count out of $trimmed_count sequences: $reason"
                        end
                    end
                end
                # Create lookup dict from sequence to allele ID
                seq_to_id = Dict(seq => id for (id, seq) in DBdict)

                # Helper funciton to remove BLAST gaps before re-aligning
                nogaps(s) = replace(s,'-'=>"")

                # Add column checking if aligned sequence exists in DB
                transform!(blast_clusters, :aln_qseq => ByRow(seq -> begin
                    query_seq = coalesce(nogaps(String(seq)), "")
                    any(db_seq -> occursin(query_seq, db_seq), keys(seq_to_id))
                end) => :isin_db)

                # Filter out truncated substrings of known alleles called incorrectly as novel
                if parsed_args["blast"]["isin"]
                    @info "Filtering out truncated substrings of known alleles"
                    truncation = filter(x -> (x.isin_db && x.aln_mismatch > 0), blast_clusters)
                    @info "Found $(nrow(truncation)) truncated substrings"
                    filter!(x -> !(x.isin_db && x.aln_mismatch > 0), blast_clusters)
                end

                # Compute trimmed counts and frequencies
                transform!(groupby(blast_clusters, [:well, :case, :aln_qseq]), :full_count => sum => :count)
                transform!(groupby(blast_clusters, [:well, :case, :gene]), :count => (x -> x ./ maximum(x)) => :frequency)

                # Optional filters
                min_fullcount = parsed_args["blast"]["minfullcount"]
                min_fullfrequency = parsed_args["blast"]["minfullfreq"]
                min_length = parsed_args["blast"]["length"]
                filter!(x->x.full_count >= min_fullcount, blast_clusters)
                @info "Clusters after filtering by full count >= $min_fullcount: $(nrow(blast_clusters)) rows"
                if verbose
                    @info "Saving clusters after filtering by full_count to $(blast_tsv)-clusters-filtered-fullcount.tsv"
                    CSV.write(blast_tsv*"-clusters-filtered-fullcount.tsv", blast_clusters)
                end
                filter!(x->x.full_frequency >= min_fullfrequency, blast_clusters)
                @info "Clusters after filtering by full frequency >= $min_fullfrequency: $(nrow(blast_clusters)) rows"
                if verbose
                    @info "Saving clusters after filtering by full frequency to $(blast_tsv)-clusters-filtered-fullfreq.tsv"
                    CSV.write(blast_tsv*"-clusters-filtered-fullfreq.tsv", blast_clusters)
                end
                filter!(x->length(x.qseq) >= min_length, blast_clusters)
                @info "Clusters after filtering by minimum length >= $min_length: $(nrow(blast_clusters)) rows"
                if verbose
                    @info "Saving clusters after filtering by minimum length to $(blast_tsv)-clusters-filtered-minlen.tsv"
                    CSV.write(blast_tsv*"-clusters-filtered-minlen.tsv", blast_clusters)
                end

                # Filter results
                if !keep_failed
                    filter!(x -> x.aln_mismatch >= 0, blast_clusters) # Drop failed trimmings
                end
                filter!(x -> x.aln_mismatch <= parsed_args["blast"]["maxdist"], blast_clusters)

                # Generate allele names
                blast_clusters[:, :allele_name] = map(
                    x -> (x.aln_mismatch == 0) ? x.sseqid : unique_name(x.sseqid, x.aln_qseq)*" Novel",
                    eachrow(blast_clusters)
                )

                # Save results
                output = cli.always_gz(parsed_args["blast"]["output"])
                @info "Saving discovered gene sequences in compressed $output file"
                CSV.write(output, blast_clusters, compress=true, delim='\t')
            end

            if get(parsed_args,"%COMMAND%","") == "exact"
                @info "Exact search"
                extension = parsed_args["exact"]["extension"]
                limit = parsed_args["exact"]["limit"]
                refgenes = parsed_args["exact"]["refgene"]
                if length(refgenes) > 0
                    @info "Using reference genes $refgenes"
                end
                table = load_demultiplex(parsed_args["exact"]["tsv"])
                if limit > 0
                    @info "Limiting number of demultiplexed reads to $limit"
                    table = table[1:limit,:]
                end
                db = load_fasta(parsed_args["exact"]["fasta"], validate=false)
                mincount = parsed_args["exact"]["mincount"]
                minratio = parsed_args["exact"]["minratio"]
                if mincount < TRUST_MINCOUNT
                    @warn "Decreasing mincount below $TRUST_MINCOUNT may lead to false positives due to sequencing errors - you've been warned!"
                end

                top = parsed_args["exact"]["top"]
                affix = parsed_args["exact"]["affix"]
                local rss
                if extension !== nothing
                    @info "Using extension mode with length $extension instead of RSS elements"
                    rss = String[]  # Empty RSS when using extension
                else
                    rss = split(parsed_args["exact"]["rss"], ',')
                    validate_types(rss)
                    @info "Extract RSS: $(join(rss,','))"
                end
                if top != 1
                    @info "Uncollapsed mode enabled; at most $top full records will be returned."
                end
                gene = parsed_args["exact"]["gene"]
                expect = parsed_args["exact"]["expect"]
                deletion = parsed_args["exact"]["deletion"]
                min_allele_case_medratio = parsed_args["exact"]["min-allele-mratio"]
                min_gene_case_medratio = parsed_args["exact"]["min-gene-mratio"]

                # Load expect if provided
                expect_df = DataFrame(name=[], ratio=[])
                if expect !== nothing
                    expect_df = CSV.read(expect, DataFrame, delim='\t')
                    @assert all([name in names(expect_df) for name in ["name","ratio"]]) "File must contain following columns: name, ratio"
                end
                if nrow(expect_df) > 0
                    @info "Using expect file with $(nrow(expect_df)) entries"
                end
                expect_dict = Dict(zip(expect_df.name, expect_df.ratio))

                # Load deletions if provided
                deletion_df = DataFrame(name=[], ratio=[])
                if deletion !== nothing
                    deletion_df = CSV.read(deletion, DataFrame, delim='\t')
                    @assert all([name in names(deletion_df) for name in ["name","ratio"]]) "File must contain following columns: name, ratio"
                end
                if nrow(deletion_df) > 0
                    @info "Using deletion_df file with $(nrow(expect_df)) entries"
                end
                deletion_df = Dict(zip(deletion_df.name, deletion_df.ratio))

                raw = parsed_args["exact"]["raw"]
                locus = parsed_args["exact"]["locus"]
                counts_df = exact.exact_search(table, db, gene, mincount=mincount, minratio=minratio, expect_dict=expect_dict, affix=affix, rss=rss, extension=extension, N=top, raw=raw)
                sort!(counts_df, [:case, :db_name])
                if !parsed_args["exact"]["noplot"]
                    if nrow(counts_df) > 0
                        println(plotgenes(counts_df))
                    else
                        @warn "No exact matches to plot"
                    end
                end
                # Add gene count frequency (aggregation with alleles starting with locus only)
                @info "Excluding genes not starting with $locus for correct allele and gene count frequency calculation"

                # Calculate gene counts per well/case/gene
                transform!(groupby(counts_df, [:well, :case, :gene])) do group
                    filtered_group = filter(row -> startswith(row.db_name, locus), group)
                    gene_count = isempty(filtered_group) ? 0 : sum(filtered_group.count)
                    return DataFrame(gene_count = fill(gene_count, nrow(group)))
                end

                # Calculate total counts per well/case
                transform!(groupby(counts_df, [:well, :case])) do group
                    filtered_group = filter(row -> startswith(row.db_name, locus), group)
                    case_count = isempty(filtered_group) ? 0 : sum(filtered_group.count)
                    return DataFrame(case_count = fill(case_count, nrow(group)))
                end

                # Calculate median count per gene across ALL cases
                transform!(groupby(counts_df, [:gene])) do group
                    filtered_group = filter(row -> startswith(row.db_name, locus), group)
                    cross_case_median_count = isempty(filtered_group) ? 0 : median(filtered_group.count)
                    return DataFrame(cross_case_median_count = fill(cross_case_median_count, nrow(group)))
                end

                # Calculate median gene_count for each gene across all cases
                transform!(groupby(counts_df, [:gene])) do group
                    filtered_group = filter(row -> startswith(row.db_name, locus), group)
                    cross_case_median_gene_count = isempty(filtered_group) ? 0 : median(filtered_group.gene_count)
                    return DataFrame(cross_case_median_gene_count = fill(cross_case_median_gene_count, nrow(group)))
                end

                # Calculate median allele count for each allele (db_name) across all cases
                transform!(groupby(counts_df, [:db_name])) do group
                    filtered_group = filter(row -> startswith(row.db_name, locus), group)
                    cross_case_median_allele_count = isempty(filtered_group) ? 0 : median(filtered_group.count)
                    return DataFrame(cross_case_median_allele_count = fill(cross_case_median_allele_count, nrow(group)))
                end

                # Frequencies within a case
                counts_df[:,:allele_case_freq] = counts_df.count ./ counts_df.case_count
                counts_df[:,:gene_case_freq] = counts_df.gene_count ./ counts_df.case_count

                # Ratios compared to cross-case medians (for detecting outliers across cases)
                counts_df[:,:allele_to_cross_case_median_ratio] = counts_df.count ./ counts_df.cross_case_median_allele_count
                counts_df[:,:gene_to_cross_case_median_ratio] = counts_df.gene_count ./ counts_df.cross_case_median_gene_count

                # Calculate allele frequency within gene (per case)
                transform!(groupby(counts_df, [:well, :case, :gene]), :count => (x->x./sum(x)) => :allele_freq)

                # Apply original filters
                filter!(x -> (x.allele_freq >= get_ratio_threshold(expect_dict, x, type="allele_freq")), counts_df)
                filter!(x -> (x.gene_case_freq >= get_ratio_threshold(deletion_df, x, type="gene_case_freq")), counts_df)

                # Use the existing min-allele-mratio and min-gene-mratio parameters for cross-case comparisons
                min_allele_cross_case_ratio = parsed_args["exact"]["min-allele-mratio"]
                min_gene_cross_case_ratio = parsed_args["exact"]["min-gene-mratio"]

                # Filter out alleles/genes that are much lower than expected across all cases
                filter!(x -> (x.allele_to_cross_case_median_ratio >= min_allele_cross_case_ratio), counts_df)
                filter!(x -> (x.gene_to_cross_case_median_ratio >= min_gene_cross_case_ratio), counts_df)

                output = cli.always_gz(parsed_args["exact"]["output"])
                # Compute refgene ratios
                if length(refgenes) > 0
                    for refgene in refgenes
                        counts_df = grouped_ratios(counts_df, refgene, count_col=:count)
                        # Transform df to sum gene count per gene
                        transform!(groupby(counts_df, [:well, :case, :gene]), :count => sum => :ref_gene_count)
                        counts_df = grouped_ratios(counts_df, refgene, count_col=:ref_gene_count)
                    end
                end
                CSV.write(output, counts_df, compress=true, delim='\t')
                @info "Exact search data saved in compressed $output file"
            end
            if get(parsed_args,"%COMMAND%","") == "exclude"
                @info "Exclude"
                db = load_fasta(parsed_args["exclude"]["fasta"])
                colname = parsed_args["exclude"]["colname"]
                @info "Using $colname column"
                colseq = parsed_args["exclude"]["colseq"]
                @info "Using $colseq column"
                data = CSV.File(parsed_args["exclude"]["input"],delim='\t') |> DataFrame
                discard = Set{Tuple{String,String,String}}()
                for row in eachrow(data)
                    for (name, seq) in db
                        if occursin(row[colseq], seq)
                            @info "Allele $(row[colname]) is a substring of $(name)"
                            push!(discard, (name, row[colname], row[colseq]))
                        end
                        if occursin(seq, row[colseq])
                            @info "Allele $(name) is a substring of $(row[colname])"
                            push!(discard, (name, row[colname], row[colseq]))
                        end
                    end
                end
                for (name, allele_name, seq) in discard
                    @info "Discarding fasta $allele_name with name $name"
                    filter!(x->x[colseq] != seq, data)
                end
                output = cli.always_gz(parsed_args["exclude"]["output"])
                CSV.write(output, data, compress=true, delim='\t')
            end
            if get(parsed_args,"%COMMAND%","") == "nwpattern"
                @info "Running nwpattern"
                tsv  = parsed_args["nwpattern"]["input"]
                fasta = parsed_args["nwpattern"]["fasta"]
                candidates_tsv = parsed_args["nwpattern"]["output"]
                filtered_candidates_tsv = first(split(parsed_args["nwpattern"]["output"], '.')) * "-final.tsv.gz"
                k = parsed_args["nwpattern"]["kmer"]
                n = parsed_args["nwpattern"]["sample"]
                max_distance = parsed_args["nwpattern"]["maxdist"]
                flank_length =  parsed_args["nwpattern"]["window"]
                min_length = parsed_args["nwpattern"]["length"]
                min_count = parsed_args["nwpattern"]["mincount"]
                min_frequency = parsed_args["nwpattern"]["minfreq"]
                max_combinations = parsed_args["nwpattern"]["max-combinations"]
                max_attempts = parsed_args["nwpattern"]["max-attempts"]
                discover_alleles(tsv, fasta, candidates_tsv, filtered_candidates_tsv,
                k=k, n=n, max_distance=max_distance,
                flank_length=flank_length,
                min_length=min_length,
                min_count=min_count,
                min_frequency=min_frequency,
                max_combinations=max_combinations,
                max_attempts=max_attempts)
            end
            if get(parsed_args,"%COMMAND%","") == "pattern"
                @info "Pattern search"
                noprofile = parsed_args["pattern"]["noprofile"] == true
                table = load_demultiplex(parsed_args["pattern"]["input"])
                db = load_fasta(parsed_args["pattern"]["fasta"])
                db_df = DataFrame(db, [:id, :seq])
                db_df[!,:gene] = map(x->replace(first(split(x.id,'*')), r"D$"=>""), eachrow(db_df)) # Remove D suffixes of genes
                if length(Set(db_df.gene)) < 20
                    @warn "Less than 20 unique genes in database may lead to poor results due to lack of sufficient gene coverage"
                end
                blacklist = DataFrame([[],[],[]], [:id, :seq,:gene])
                if parsed_args["pattern"]["blacklist"] !== nothing
                    blacklist = DataFrame(load_fasta(parsed_args["pattern"]["blacklist"]),[:id, :seq])
                    blacklist[!,:gene] = map(x->replace(first(split(x.id,'*')), r"D$"=>""), eachrow(blacklist)) # Remove D suffixes of genes
                end
                final = patterns.search_patterns(table, blacklist, db_df, fragment_size=parsed_args["pattern"]["kmer"], max_fragment_size=parsed_args["pattern"]["maxkmer"], max_fragments=parsed_args["pattern"]["sample"], weights=parsed_args["pattern"]["weights"], mlen=parsed_args["pattern"]["length"], min_freq=parsed_args["pattern"]["minfreq"], min_count=parsed_args["pattern"]["mincount"], max_dist=parsed_args["pattern"]["maxdist"], noprofile=noprofile, usehalf=parsed_args["pattern"]["usehalf"])
                # Output path
                output = cli.always_gz(parsed_args["pattern"]["output"])
                # Save top candidates
                output_basename, output_extension = split(output,'.',limit=2)
                top_output = "$(output_basename)-top.$output_extension"
                # Mandatory length filter
                filter!(x->x.length >= parsed_args["pattern"]["length"], final)
                top = parsed_args["pattern"]["top"]
                if top > 0
                    @info "Saving top $top candidates per allele_name and case in $top_output"
                    top_per_gene = combine(groupby(final, [:gene, :allele_name, :case]), df -> sort(df, :count, rev=true)[1:min(top, nrow(df)), :])
                    dropmissing!(top_per_gene)
                    CSV.write(top_output, top_per_gene, compress=true, delim='\t')
                end
                # Filter by criteria
                filter!(x->x.count >= parsed_args["pattern"]["mincount"], final)
                filter!(x->x.dist <= parsed_args["pattern"]["maxdist"], final)
                filter!(x->x.ratio >= parsed_args["pattern"]["minfreq"], final)
                # Filter pseudo genes
                patterns.filter_and_report!(final, blacklist)

                CSV.write(output, final, compress=true, delim='\t')
                @info "Pattern search data saved in compressed $output file"
            end
            if get(parsed_args,"%COMMAND%","") == "linkage"
                @info "Linkage analysis"
                tsv  = parsed_args["linkage"]["input"]
                edges_out  = cli.always_gz(parsed_args["linkage"]["edges"])
                case_col = parsed_args["linkage"]["case-col"]
                allele_col = parsed_args["linkage"]["allele-col"]
                min_donors = parsed_args["linkage"]["min-donors"]
                min_support = parsed_args["linkage"]["min-support"]
                jaccard_threshold = parsed_args["linkage"]["min-jaccard"]
                sim_mode = parsed_args["linkage"]["similarity"] == "r2" ? :r2 : :r
                sim_thresh = parsed_args["linkage"]["threshold"]
                min_cluster_size = parsed_args["linkage"]["min-cluster-size"]
                clusters_path = get(parsed_args["linkage"], "clusters", nothing)
                df = CSV.File(tsv, delim='\t') |> DataFrame
                @info "Loaded $(nrow(df)) rows from input file"
                
                edges, clusters, clusters_detailed = linkage.compute_edges_and_clusters(
                    df; case_col=case_col, allele_col=allele_col, min_donors=min_donors,
                    min_support=min_support, jaccard_threshold=jaccard_threshold,
                    similarity=sim_mode, similarity_threshold=sim_thresh, min_cluster_size=min_cluster_size
                )
                
                CSV.write(edges_out, edges, compress=true, delim='\t')
                @info "Edges saved in compressed $edges_out file ($(nrow(edges)) pairs)"
                
                if clusters_path !== nothing
                    # Create detailed clusters output with donor information and statistics
                    all_cluster_rows = Vector{NamedTuple}()
                    for (gid, cluster_df) in enumerate(clusters_detailed)
                        for row in eachrow(cluster_df)
                            push!(all_cluster_rows, (
                                group_id = gid,
                                allele = row.allele,
                                donors = row.donors,
                                n_donors = row.n_donors,
                                mean_n11 = row.mean_n11,
                                mean_n10 = row.mean_n10,
                                mean_n01 = row.mean_n01,
                                mean_n00 = row.mean_n00,
                                mean_r = row.mean_r,
                                max_r = row.max_r,
                                min_r = row.min_r
                            ))
                        end
                    end
                    
                    # Also create simple clusters file for backward compatibility
                    simple_rows = Vector{NamedTuple}()
                    for (gid, comp) in enumerate(clusters)
                        for allele in comp
                            push!(simple_rows, (; group_id=gid, allele))
                        end
                    end
                    
                    # Save detailed clusters
                    CSV.write(clusters_path, DataFrame(all_cluster_rows), delim='\t')
                    @info "Detailed clusters saved in $(clusters_path) ($(length(clusters)) clusters, $(length(all_cluster_rows)) alleles)"
                    
                    # Save simple clusters as well (with .simple suffix)
                    simple_path = replace(clusters_path, r"\.tsv$" => ".simple.tsv")
                    CSV.write(simple_path, DataFrame(simple_rows), delim='\t')
                    @info "Simple clusters saved in $(simple_path) for backward compatibility"
                end
            end
            if get(parsed_args,"%COMMAND%","") == "novel"
                df = CSV.File(parsed_args["novel"]["tsv"],delim='\t') |> DataFrame
                # Assert that the input file has a column named allele_name and seq
                @assert "allele_name" ∈ names(df) "Input file must have allele_name column"
                @assert "seq" ∈ names(df) "Input file must have seq column"
                for row in eachrow(df[contains.(df.allele_name, "Novel"),:])
                    allele_name = rstrip(replace(row.allele_name, "Novel" => ""))
                    println(">$allele_name\n$(row.seq)")
                end
            end

            if get(parsed_args,"%COMMAND%","") == "regex"
                @info "Regex method to extract alleles"
                db = load_fasta(parsed_args["regex"]["fasta"])
                @info "Loaded $(length(db)) query sequences"
                table = load_demultiplex(parsed_args["regex"]["tsv"])
                table = table[.!occursin.("N",table.genomic_sequence),:]
                prefix_sorted, suffix_sorted = search_frequent_flanks(table, db, min_frequency=parsed_args["regex"]["flank-frequency"], min_length=parsed_args["regex"]["insert-minlen"], min_count=parsed_args["regex"]["flank-mincount"], min_well=parsed_args["regex"]["flank-minwell"], nprefix=parsed_args["regex"]["nprefix"], nsuffix=parsed_args["regex"]["nsuffix"])
                @info "Discarding regex flanks that are part of the known sequence"
                # Collect heptamers that are not part of D gene itself (this is a problem because some Ds are substring of other Ds)
                plist = filter(x->!any(occursin.(first(x),map(y->y[2], db))), prefix_sorted)
                slist = filter(x->!any(occursin.(first(x),map(y->y[2], db))), suffix_sorted)
                @assert length(plist) > 0 "No prefix flanks found - interrupting"
                @assert length(slist) > 0 "No suffix flanks found - interrupting"
                @info "Found prefixes $plist"
                @info "Found suffixes $slist"
                @info "Regex searching reads"
                data = search_regex(table, plist, slist, minlen=parsed_args["regex"]["insert-minlen"], maxlen=parsed_args["regex"]["insert-maxlen"])
                table[:,:hit] = map(x->x.hit,data)
                table[:,:hit_prefix] = map(x->x.prefix,data)
                table[:,:hit_suffix] = map(x->x.suffix,data)
                # Discard uninformative sequences that didn't match
                agumented_table = filter(x->(x.hit_suffix !== nothing) && (x.hit_prefix !== nothing), table)
                @info "Out of $(nrow(table)) sequences, $(nrow(agumented_table)) matched regex ($(round((nrow(agumented_table)/nrow(table))*100,digits=2))%)"
                # Annotate and collapse
                if nrow(agumented_table) == 0
                    @warn "No sequences matched regex - interrupting"
                    return
                end
                result = annotate(agumented_table, db, nprefix=parsed_args["regex"]["nprefix"], nsuffix=parsed_args["regex"]["nsuffix"])
                collapsed = combine(groupby(result, [:well, :case, :best_name, :distance, :prefix, :best_aln, :suffix]), nrow => :count)
                sort!(collapsed, [:well, :case,:best_name, :count], rev=[false, false, false, true])
                transform!(collapsed, :best_name => ByRow(x -> first(split(x, '*'))) => :gene)
                transform!(groupby(collapsed, [:well, :case, :gene]), :count => (x->x./maximum(x))=> :frequency)
                result_df = filter(r->(r.frequency >= parsed_args["regex"]["frequency"]) & (r.count > parsed_args["regex"]["mincount"]) & (r.distance <= parsed_args["regex"]["maxdist"]), collapsed)
                update_names!(result_df)
                output = cli.always_gz(parsed_args["regex"]["output"])
                CSV.write(output, result_df, compress=true, delim='\t')
                @info "Extracted Ds saved in compressed $output file"
            end
            if get(parsed_args,"%COMMAND%","") == "hash"
                @info "Hasing alleles"
                fastain = parsed_args["hash"]["fastain"]
                db = load_fasta(fastain)
                for (name, seq) in db
                    newname = unique_name(name, seq)
                    println(">$newname\n$seq")
                end
            end
            if get(parsed_args,"%COMMAND%","") == "collect"
                @info "Collecting TSV files"
                pattern = parsed_args["collect"]["pattern"]
                output = parsed_args["collect"]["output"]
                files = glob(pattern)
                @info "Found $(length(files)) files matching pattern $pattern"
                collected = []
                if length(files) > 0
                    @info "Collecting files"
                    first_file_columns = nothing
                    for file in files
                        df = CSV.File(file, delim='\t') |> DataFrame
                        df[!,:file] .= file
                        if first_file_columns === nothing
                            first_file_columns = names(df)
                        else
                            @assert first_file_columns == names(df) "Column names in $file do not match first file"
                        end
                        push!(collected, df)
                    end
                    @info "Saving collected data in $output"
                    collectd_df = vcat(collected...)
                    CSV.write(output, collectd_df, delim='\t', compress=true)
                else
                    @warn "No files found matching pattern $pattern"
                end
            end
            if get(parsed_args,"%COMMAND%","") == "bwa"
                @info "BWA search to filter candidates if they match correct chromosome"
                df = CSV.File(parsed_args["bwa"]["tsv"], delim='\t') |> DataFrame
                chromosome_name = parsed_args["bwa"]["chromosome"]
                outtsv = parsed_args["bwa"]["output"]
                genome = parsed_args["bwa"]["genome"]
                colname = parsed_args["bwa"]["colname"]
                colseq = parsed_args["bwa"]["colseq"]
                tag = parsed_args["bwa"]["tag"]
                @info "Using columns: $colname, $colseq"
                @info "Using genome: $genome"
                @info "Filter by chromosome: $chromosome_name and $tag tag"
                sequences = map(eachrow(df)) do row
                    concatenated_sequence = concatenate_columns(row, colseq)
                    (row[colname], concatenated_sequence)
                end
                indices, position = bwa_sequences(genome, sequences, chromosome_name, tag=tag)
                df[!,:position] = position
                @info "$(length(sequences))"
                @info "$(nrow(df)) sequences matched chromosome $chromosome_name"
                @info "$(length(indices)) sequences matched chromosome $chromosome_name"
                @info "$(indices[1:10])"
                @info "$(df[1:10,:])"

                CSV.write(outtsv, df[indices,:], compress=true, delim='\t')
                @info "Filtered result saved in $outtsv"
            end
            if get(parsed_args,"%COMMAND%","") == "hamming"
                @info "Hamming distance window search"
                limit = parsed_args["hamming"]["limit"]
                assignments = parsed_args["hamming"]["assignments"]
                if limit > 0
                    @info "Limiting number of reads to $limit"
                    table = load_demultiplex(parsed_args["hamming"]["tsv"], limit=limit)
                else
                    table = load_demultiplex(parsed_args["hamming"]["tsv"])
                end
                db = load_fasta(parsed_args["hamming"]["fasta"])
                umi = parsed_args["hamming"]["umi"]
                column = parsed_args["hamming"]["column"]
                @info "Using $column column"
                checkbounds = true
                if column != "genomic_sequence"
                    checkbounds = false
                end
                assignments_df = hamming_search(table, db, max_dist=parsed_args["hamming"]["maxdist"], column=column, check_bounds=checkbounds, umi=umi)
                if assignments !== nothing
                    @info "Saving intermediate assignments"
                    CSV.write(assignments, assignments_df, delim='\t', compress=true)
                end
                summarized = hamming.summarize(assignments_df, db, min_count=parsed_args["hamming"]["mincount"], cluster_ratio=parsed_args["hamming"]["ratio"])
                output = cli.always_gz(parsed_args["hamming"]["output"])
                CSV.write(output, summarized, compress=true, delim='\t')
                @info "Hamming search data saved in compressed $output file"
                if !parsed_args["hamming"]["noplot"] & nrow(summarized) > 0
                    println(plotgenes(summarized))
                end
            end
        end
    end

    """
        julia_main(args=[])::Cint

    Entry point C call for main function for Julia
    """
    function julia_main()::Cint
        try
            real_main(ARGS)
        catch
            Base.invokelatest(Base.display_error, Base.catch_stack())
            return 1
        end
        return 0
    end
end
