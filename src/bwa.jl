module bwa
    using BurrowsWheelerAligner
    using BurrowsWheelerAligner.LibBWA: mem_aln_t
    using BurrowsWheelerAligner: Aligner
    using FASTX
    import FASTX.description
    using ProgressMeter
    using DataStructures

    function description(aln::BurrowsWheelerAligner.LibBWA.mem_aln_t, aligner::BurrowsWheelerAligner.Aligner)
        anns = BurrowsWheelerAligner.LibBWA.unsafe_load(aligner.index.bns).anns
        anno = BurrowsWheelerAligner.LibBWA.unsafe_load(anns, aln.rid+1).anno
        BurrowsWheelerAligner.LibBWA.unsafe_string(anno)
    end

    function create_aligner(genome_path)::Vector{Tuple{String, BurrowsWheelerAligner.Aligner}}
        return [(path, BurrowsWheelerAligner.Aligner(path)) for path in genome_path]
    end

    function bwa_sequences(genome_path, sequences, chromosome_name)
        result = zeros(Bool, length(sequences))
        aligners = create_aligner(genome_path)
        discard = Accumulator{Tuple{String,String,String}, Int}()
        position = fill("",length(sequences))
        @showprogress for (nallele, (name, sequence)) in enumerate(sequences)
            record = FASTA.Record(name, sequence)
            for (genome_file, aligner) in aligners
                alns = BurrowsWheelerAligner.align(aligner, record)
                # Pick alignments with the highest score
                if length(alns) == 0
                    @info "$name does **NOT** align to the $genome_file uniquely or at all (skipping)"
                    push!(result, false)
                    continue
                end
                max_score = maximum(map(x->(x.score), alns))
                alns = filter(x->x.score == max_score, alns)
                descriptions = map(x->description(x, aligner), alns)
                match = filter(x->occursin(chromosome_name, x), descriptions)
                nomatch = filter(x->!occursin(chromosome_name, x), descriptions)
                if length(match) > 0
                    result[nallele] = true
                    position[nallele] = "$(BurrowsWheelerAligner.refname(first(alns), aligner)):$(first(alns).pos)"
                else
                    bad_chromosome = first(nomatch)
                    push!(discard, (genome_file, name, bad_chromosome))
                    result[nallele] = false
                end
            end
        end
        for ((genome_file, name, chr),n) in discard
            @info "Discarded $name matching $chr in $genome_file (total $n)"
        end
        return result, position
    end

    export bwa_sequences
end