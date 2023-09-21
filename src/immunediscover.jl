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

    using .cli
    using .demultiplex
    using .simulate
    using .data
    using .profile
    using .trim
    using .exact
    using .heptamer
    using .hamming

    using CSV
    using DataFrames
    using JSON
    export load_fasta

    function real_main(args=[])
        parsed_args = parse_commandline(args)
        if parsed_args != nothing
            if get(parsed_args,"%COMMAND%","") == "demultiplex"
                @info "Demultiplexing"
                table = immunediscover.demultiplex.demux(parsed_args["demultiplex"]["fastq"],
                                                         parsed_args["demultiplex"]["indices"],
                                                         min_length=parsed_args["demultiplex"]["length"])
                output = cli.always_gz(parsed_args["demultiplex"]["output"])
                CSV.write(output, table, compress=true, delim='\t')
                @info "Demultiplexed data saved in compressed $output file"
            end

            if get(parsed_args,"%COMMAND%","") == "heptamer"
                @info "Extracting heptamers"
                if !isfile(parsed_args["heptamer"]["json"])
                    @info "No heptamer file provided, creating one"
                    heptamers = Dict("IGKV" => ["CACAGTG","CACACTG","CACTGTG","CACGGTG","CACAATG","CACATTG"],
                                     "IGLV" => ["CACAGTG","CACGGTG","CATGGTG","CACGCTG","CACAGCG","CACAGTA","CATAGTG","CACAATG"],
                                     "IGHV" => ["CACAGTG","CACAATG","CACAGAG","CACGGTG","CACAGCG"])
                    open(parsed_args["heptamer"]["json"],"w") do f
                        JSON.print(f, JSON.json(heptamers))
                    end
                    @info "Created heptamer JSON $(parsed_args["heptamer"]["json"]) file"
                end
                if isfile(parsed_args["heptamer"]["json"])
                    heptamers = JSON.parsefile(parsed_args["heptamer"]["json"])
                    @assert all([k ∈ ["IGKV","IGLV","IGHV"] for k in keys(heptamers)]) "JSON file with heptamers must contain IGKV, IGLV and IGHV keys"
                end
                @info "Loaded heptamers for $(join([k for k in keys(heptamers)],','))"
                chain = parsed_args["heptamer"]["chain"]
                @info "Parsing for $(chain)"
                @info "Using heptamers $(join(heptamers[chain],','))"
                db = load_fasta(parsed_args["heptamer"]["fasta"])
                @info "Loaded $(length(db)) query sequences"
                table = CSV.File(parsed_args["heptamer"]["tsv"],delim='\t') |> DataFrame
                data = extract_heptamers(table, db, heptamers[chain]; max_dist=parsed_args["heptamer"]["maxdist"], b=parsed_args["heptamer"]["begin"]+1, e=parsed_args["heptamer"]["end"])
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
                table = CSV.File(parsed_args["trim"]["input"], delim='\t') |> DataFrame
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
                table = CSV.File(parsed_args["exact"]["tsv"], delim='\t') |> DataFrame
                db = data.load_fasta(parsed_args["exact"]["fasta"])
                counts_df = exact.exact_search(table, db)
                sort!(counts_df, [:case, :db_name])
                output = cli.always_gz(parsed_args["exact"]["output"])
                CSV.write(output, counts_df, compress=true, delim='\t')
                @info "Exact search data saved in compressed $output file"
            end

            if get(parsed_args,"%COMMAND%","") == "hamming"
                @info "Hamming distance window search"
                limit = parsed_args["hamming"]["limit"]
                assignments = parsed_args["hamming"]["assignments"]
                if limit > 0
                    @info "Limiting number of reads to $limit"
                    table = CSV.File(parsed_args["hamming"]["tsv"], delim='\t', types=Dict(:case => String), limit=limit) |> DataFrame
                else
                    table = CSV.File(parsed_args["hamming"]["tsv"], delim='\t', types=Dict(:case => String)) |> DataFrame
                end
                db = load_fasta(parsed_args["hamming"]["fasta"])
                umi = parsed_args["hamming"]["umi"]
                column = parsed_args["hamming"]["column"] == "trimmed_sequence" ? :trimmed_sequence : :genomic_sequence
                @info "Using $column column"
                checkbounds = true
                if column != :genomic_sequence
                    checkbounds = false
                end
                found = hamming_search(table, db, max_dist=parsed_args["hamming"]["maxdist"], column=column, check_bounds=checkbounds, umi=umi)
                if assignments !== nothing
                    @info "Saving intermediate assignments"
                    assignments_df = DataFrame(vcat(found...))
                    rename!(assignments_df,[:closest_name, :distance,:start,:stop,:prefix,:middle,:suffix,:db_sequence,:case,:barcode])
                    CSV.write(assignments, assignments_df, delim='\t')
                end
                summarized = hamming.summarize(found, db, min_count=parsed_args["hamming"]["mincount"], cluster_ratio=parsed_args["hamming"]["ratio"])
                output = cli.always_gz(parsed_args["hamming"]["output"])
                CSV.write(output, summarized, compress=true, delim='\t')
                @info "Hamming search data saved in compressed $output file"
                if parsed_args["hamming"]["plot"] !== nothing
                    @info "Saving plot"
                    hamming.plot(summarized, parsed_args["hamming"]["plot"])
                    @info "Plot saved in $(parsed_args["hamming"]["plot"]) file"
                end
            end
        end
    end

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
