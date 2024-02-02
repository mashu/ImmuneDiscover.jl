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

    using .cli
    using .demultiplex
    using .simulate
    using .data
    using .profile
    using .trim
    using .exact
    using .heptamer
    using .hamming
    using .patterns
    using .regex
    using .bwa

    using CSV
    using DataFrames
    using UnicodePlots
    using Glob

    export load_fasta

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

    """
        real_main(args=[])
    
    Main function for Julia
    """
    function real_main(args=[])
        parsed_args = parse_commandline(args)
        if parsed_args != nothing
            if get(parsed_args,"%COMMAND%","") == "demultiplex"
                @info "Demultiplexing"
                table, stats = immunediscover.demultiplex.demux(parsed_args["demultiplex"]["fastq"],
                                                         parsed_args["demultiplex"]["indices"],
                                                         min_length=parsed_args["demultiplex"]["length"],
                                                         save_fastq_files=parsed_args["demultiplex"]["split"])
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

            if get(parsed_args,"%COMMAND%","") == "exact"
                @info "Exact search"
                table = load_demultiplex(parsed_args["exact"]["tsv"])
                db = load_fasta(parsed_args["exact"]["fasta"], validate=false)
                mincount = parsed_args["exact"]["mincount"]
                minfreq = parsed_args["exact"]["minfreq"]
                full_mincount = parsed_args["exact"]["full-mincount"]
                full_minfreq = parsed_args["exact"]["full-minfreq"]
                top = parsed_args["exact"]["top"]
                affix = parsed_args["exact"]["affix"]
                types = split(parsed_args["exact"]["types"], ',')
                validate_types(types)
                @info "Extract RSS types: $(join(types,','))"
                if top != 1
                    @info "Uncollapsed mode enabled; at most $top full records will be returned."
                end
                gene = parsed_args["exact"]["gene"]
                counts_df = exact.exact_search(table, db, gene, mincount=mincount, minfreq=minfreq, full_mincount=full_mincount, full_minfreq=full_minfreq, affix=affix, rss=types, N=top)
                sort!(counts_df, [:case, :db_name])
                if !parsed_args["exact"]["noplot"]
                    println(plotgenes(counts_df))
                end
                output = cli.always_gz(parsed_args["exact"]["output"])
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
                        if occursin(row[colseq], seq) || occursin(seq, row[colseq])
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
                final = patterns.search_patterns(table, blacklist, db_df, fragment_size=parsed_args["pattern"]["kmer"], max_fragment_size=parsed_args["pattern"]["maxkmer"], max_fragments=parsed_args["pattern"]["sample"], weights=parsed_args["pattern"]["weights"], mlen=parsed_args["pattern"]["length"], min_freq=parsed_args["pattern"]["minfreq"], min_count=parsed_args["pattern"]["mincount"], max_dist=parsed_args["pattern"]["maxdist"], noprofile=noprofile)
                # Output path
                output = cli.always_gz(parsed_args["pattern"]["output"])
                # Save top candidates
                output_basename, output_extension = split(output,'.',limit=2)
                top_output = "$(output_basename)-top.$output_extension"
                filter!(x->x.length >= parsed_args["pattern"]["length"], final)
                filter!(x->x.count >= parsed_args["pattern"]["mincount"], final)
                top_per_gene = combine(groupby(final, [:gene, :case]), df -> sort(df, :count, rev=true)[1:min(parsed_args["pattern"]["top"], nrow(df)), :])
                dropmissing!(top_per_gene)
                CSV.write(top_output, top_per_gene, compress=true, delim='\t')
                # Filter by criteria
                filter!(x->x.dist <= parsed_args["pattern"]["maxdist"], final)
                filter!(x->x.ratio >= parsed_args["pattern"]["minfreq"], final)
                # Filter pseudo genes
                patterns.filter_and_report!(final, blacklist)
 
                CSV.write(output, final, compress=true, delim='\t')
                @info "Pattern search data saved in compressed $output file"
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
                    CSV.write(output, collectd_df, delim='\t')
                else
                    @warn "No files found matching pattern $pattern"
                end
            end
            if get(parsed_args,"%COMMAND%","") == "bwa"
                @info "BWA search to filter candidates if they match correct chromosome"
                df = CSV.File(parsed_args["bwa"]["tsv"],delim='\t') |> DataFrame
                chromosome_name = parsed_args["bwa"]["chromosome"]
                tsv = CSV.File(parsed_args["bwa"]["tsv"],delim='\t') |> DataFrame
                outtsv = parsed_args["bwa"]["output"]
                genome = parsed_args["bwa"]["genome"]
                colname = parsed_args["bwa"]["colname"]
                colseq = parsed_args["bwa"]["colseq"]
                @info "Using columns: $colname, $colseq"
                @info "Using genome: $genome"
                @info "Filter by chromosome: $chromosome_name"
                sequences = map(eachrow(df)) do row
                    concatenated_sequence = concatenate_columns(row, colseq)
                    (row[colname], concatenated_sequence)
                end
                indices, position = bwa_sequences(genome, sequences, chromosome_name)
                tsv[!,:position] = position
                CSV.write(tsv[indices,:], outtsv, compress=true, delim='\t')
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
                    CSV.write(assignments, assignments_df, delim='\t')
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
