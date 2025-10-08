module immunediscover
    include("cli.jl")
    include("demultiplex.jl")
    include("simulate.jl")
    include("data.jl")
    include("profile.jl")
    # include("trim.jl")  # deprecated command; not compiled
    include("exact.jl")
    include("heptamer.jl")
    # include("hamming.jl")  # deprecated command; not compiled
    # include("patterns.jl")  # deprecated command; not compiled
    # include("regex.jl")  # deprecated command; not compiled
    include("hsmm.jl")
    include("bwa.jl")
    # include("nwpattern.jl")  # deprecated command; not compiled
    include("blast.jl")
    include("keyedsets.jl")
    include("linkage.jl")
    include("haplotype.jl")
    include("fasta.jl")
    include("merge.jl")
    include("table.jl")

    using .cli
    using .demultiplex
    using .simulate
    using .data
    using .profile
    # using .trim  # deprecated
    using .exact
    using .heptamer
    # using .hamming  # deprecated
    # using .nwpattern  # deprecated
    # using .patterns  # deprecated
    # using .regex  # deprecated
    using .HSMM
    using .bwa
    using .blast
    using .keyedsets
    using .linkage
    using .haplotype
    using .fasta
    using .merge
    using .table

    using CSV
    using DataFrames
    using UnicodePlots
    using Glob
    using Statistics
    using DataStructures
    using FASTX

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
            # demultiplex (top-level)
            if get(parsed_args,"%COMMAND%","") == "demultiplex"
                demultiplex.handle_demultiplex(parsed_args, cli.always_gz)
            end

            # search group
            if get(parsed_args,"%COMMAND%","") == "search"
                subcmd = get(parsed_args["search"], "%COMMAND%", "")
                if subcmd == "heptamer"
                    heptamer.handle_heptamer(parsed_args, immunediscover, cli.always_gz)
                elseif subcmd == "exact"
                    exact.handle_exact(parsed_args, immunediscover, cli.always_gz)
                elseif subcmd == "hsmm"
                    HSMM.handle_hsmm(parsed_args)
                elseif subcmd == "blast"
                    parsed_args = apply_blast_presets!(parsed_args)
                    gene = parsed_args["search"]["blast"]["gene"]
                    if haskey(BLAST_PRESETS, gene)
                        @info "Initial $gene preset parameters (can be overridden by command line arguments):"
                        for (param, value) in BLAST_PRESETS[gene]
                            @info "  --$param = $value"
                        end
                    end
                    @info "Running with the following parameters:"
                    for (param, value) in parsed_args["search"]["blast"]
                        @info "  --$param = $value"
                    end
                    @info "Discovery with BLAST assignments"
                    # Load paths and parameters
                    fasta_path = parsed_args["search"]["blast"]["fasta"]
                    affixes_path = fasta_path * ".affixes"
                    DB = load_fasta(fasta_path, validate=false)
                    # Parameters
                    verbose = parsed_args["search"]["blast"]["verbose"]
                    overwrite = parsed_args["search"]["blast"]["overwrite"]
                    minquality = parsed_args["search"]["blast"]["minquality"]
                    forward_extension = parsed_args["search"]["blast"]["forward"]
                    reverse_extension = parsed_args["search"]["blast"]["reverse"]
                    keep_failed = parsed_args["search"]["blast"]["keep-failed"]
                    min_corecov = parsed_args["search"]["blast"]["min-corecov"]
                    if (forward_extension < 7) && (forward_extension > 0)
                        @warn "Forward extension $forward_extension is short and may lead to false positives due to problem with ambigous trimming alignment - you've been warned!"
                    end
                    if (reverse_extension < 7) && (reverse_extension > 0)
                        @warn "Reverse extension $reverse_extension is short and may lead to false positives due to problem with ambigous trimming alignment - you've been warned!"
                    end
                    # Process pseudo genes if provided
                    db_p = Vector{Tuple{String, String}}()
                    pseudo = parsed_args["search"]["blast"]["pseudo"]
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
                        ext_fasta_path = combined_fasta_path
                        @info "No sequence extension requested, using original sequences"
                        empty_affixes = [(name="", prefix="", suffix="")]
                        CSV.write(affixes_path, DataFrame(empty_affixes), delim='\t')
                    else
                        ext_fasta_path = blast.replace_extension(combined_fasta_path, "fasta", tag="-extended")
                        if !isfile(ext_fasta_path) || overwrite
                            @info "Extending gene sequences by $forward_extension forward and $reverse_extension reverse nucleotides"
                            demux = blast.load_csv(parsed_args["search"]["blast"]["input"])
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
                        parsed_args["search"]["blast"]["input"],
                        ext_fasta_path,
                        max_dist=parsed_args["search"]["blast"]["maxdist"],
                        min_edge=parsed_args["search"]["blast"]["edge"],
                        min_scov=parsed_args["search"]["blast"]["subjectcov"],
                        args=parsed_args["search"]["blast"]["args"],
                        verbose=verbose,
                        overwrite=overwrite
                    )
                    if forward_extension != 0 || reverse_extension != 0
                        if isfile(affixes_path)
                            @info "Loading affixes from $affixes_path"
                            affixes = Tuple.(collect.(eachrow(CSV.read(affixes_path, DataFrame, delim='\t'))))
                            if !any(occursin.("prefix", names(blast_clusters)))
                                @info "Using affixes to trim extended sequences"
                                affixes_df = DataFrame(affixes, [:sseqid, :prefix, :suffix])
                                leftjoin!(blast_clusters, affixes_df, on=:sseqid)
                                blast_clusters.prefix = coalesce.(blast_clusters.prefix, "")
                                blast_clusters.suffix = coalesce.(blast_clusters.suffix, "")
                            end
                        end
                    else
                        if !any(occursin.("prefix", names(blast_clusters)))
                            blast_clusters.prefix .= ""
                            blast_clusters.suffix .= ""
                        end
                    end
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
                    if !keep_failed
                        before = nrow(blast_clusters)
                        filter!(x -> !ismissing(x.aln_qseq) && (x.aln_qseq != "") && (x.aln_qseq != "nothing") && x.aln_mismatch >= 0, blast_clusters)
                        @info "Dropped $(before - nrow(blast_clusters)) rows with failed/empty alignments (remaining: $(nrow(blast_clusters)))"
                    end
                    transform!(blast_clusters, [:aln_qseq, :db_seq] => ByRow((aq, ds) -> begin
                        aq_s = coalesce(String(aq), "")
                        ds_s = coalesce(String(ds), "")
                        if isempty(aq_s) || isempty(ds_s)
                            0.0
                        else
                            length(aq_s) / length(ds_s)
                        end
                    end) => :corecov)
                    before_corecov = nrow(blast_clusters)
                    filter!(x -> x.corecov >= min_corecov, blast_clusters)
                    @info "Dropped $(before_corecov - nrow(blast_clusters)) rows with core coverage < $min_corecov (remaining: $(nrow(blast_clusters)))"
                    # Lookup and flags
                    seq_to_id = Dict(seq => id for (id, seq) in DBdict)
                    nogaps(s) = replace(s,'-'=>"")
                    transform!(blast_clusters, :aln_qseq => ByRow(seq -> begin
                        query_seq = coalesce(nogaps(String(seq)), "")
                        any(db_seq -> occursin(query_seq, db_seq), keys(seq_to_id))
                    end) => :isin_db)
                    if parsed_args["search"]["blast"]["isin"]
                        @info "Filtering out truncated substrings of known alleles"
                        truncation = filter(x -> (x.isin_db && x.aln_mismatch > 0), blast_clusters)
                        @info "Found $(nrow(truncation)) truncated substrings"
                        filter!(x -> !(x.isin_db && x.aln_mismatch > 0), blast_clusters)
                    end
                    transform!(groupby(blast_clusters, [:well, :case, :aln_qseq]), :full_count => sum => :count)
                    transform!(groupby(blast_clusters, [:well, :case, :gene]), :count => (x -> x ./ maximum(x)) => :frequency)
                    min_fullcount = parsed_args["search"]["blast"]["minfullcount"]
                    min_fullfrequency = parsed_args["search"]["blast"]["minfullfreq"]
                    min_length = parsed_args["search"]["blast"]["length"]
                    filter!(x->x.full_count >= min_fullcount, blast_clusters)
                    filter!(x->x.full_frequency >= min_fullfrequency, blast_clusters)
                    filter!(x->length(x.qseq) >= min_length, blast_clusters)
                    if !keep_failed
                        filter!(x -> x.aln_mismatch >= 0, blast_clusters)
                    end
                    filter!(x -> x.aln_mismatch <= parsed_args["search"]["blast"]["maxdist"], blast_clusters)
                    blast_clusters[:, :allele_name] = map(
                        x -> (x.aln_mismatch == 0) ? x.sseqid : unique_name(x.sseqid, x.aln_qseq)*" Novel",
                        eachrow(blast_clusters)
                    )
                    output = cli.always_gz(parsed_args["search"]["blast"]["output"])
                    @info "Saving discovered gene sequences in compressed $output file"
                    CSV.write(output, blast_clusters, compress=true, delim='\t')
                end
            end

            # analyze group
            if get(parsed_args, "%COMMAND%", "") == "analyze"
                subcmd = get(parsed_args["analyze"], "%COMMAND%", "")
                if subcmd == "linkage"
                    linkage.handle_linkage(parsed_args, cli.always_gz)
                elseif subcmd == "haplotype"
                    haplotype.handle_haplotype(parsed_args)
                elseif subcmd == "bwa"
                    bwa.handle_bwa(parsed_args, immunediscover, cli.always_gz)
                end
            end

            # table group
            if get(parsed_args,"%COMMAND%","") == "table"
                table.handle_table(parsed_args, immunediscover, cli.always_gz)
            end

            # fasta group
            if get(parsed_args, "%COMMAND%", "") == "fasta"
                subcmd = get(parsed_args["fasta"], "%COMMAND%", "")
                if subcmd == "merge"
                    merge.handle_merge(parsed_args)
                elseif subcmd == "diff"
                    fasta.handle_fasta_diff(parsed_args, immunediscover)
                elseif subcmd == "hash"
                    fasta.handle_fasta_hash(parsed_args, immunediscover)
                end
            end
        end  # close if parsed_args !== nothing
    end  # close function real_main

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
