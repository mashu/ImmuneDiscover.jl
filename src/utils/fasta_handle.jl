function handle_fasta_diff(parsed_args, immunediscover_module)
    fasta_paths = parsed_args["fasta"]["diff"]["fasta"]
    fasta_files = [(file=file, records=immunediscover_module.load_fasta.(file)) for file in fasta_paths]
    sets = [(file=x, set=immunediscover_module.KeyedSet(reverse.(y))) for (x,y) in fasta_files]
    for i in eachindex(sets)
        for j in i+1:length(sets)
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

function handle_fasta_hash(parsed_args, immunediscover_module)
    @info "Hashing alleles"
    fastain = parsed_args["fasta"]["hash"]["fastain"]
    db = immunediscover_module.load_fasta(fastain)
    for (name, seq) in db
        newname = immunediscover_module.unique_name(name, seq)
        println(">$(newname)\n$(seq)")
    end
end
