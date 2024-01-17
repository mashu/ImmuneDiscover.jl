module bwa
    using BurrowsWheelerAligner
    using BurrowsWheelerAligner.LibBWA: mem_aln_t
    using BurrowsWheelerAligner: Aligner
    using FASTX
    import FASTX.description
    using ProgressMeter
    
    function description(aln::BurrowsWheelerAligner.LibBWA.mem_aln_t, aligner::BurrowsWheelerAligner.Aligner)
        anns = BurrowsWheelerAligner.LibBWA.unsafe_load(aligner.index.bns).anns
        anno = BurrowsWheelerAligner.LibBWA.unsafe_load(anns, aln.rid+1).anno
        BurrowsWheelerAligner.LibBWA.unsafe_string(anno)
    end

    function bwa_sequences(genome_path, sequences, chromosome_name)
        result = Vector{Bool}()
        aligner = BurrowsWheelerAligner.Aligner(genome_path)
        discard = Set{Tuple{String,String}}()
        @showprogress for (name, sequence) in sequences
            record = FASTA.Record(name, sequence)
            alns = BurrowsWheelerAligner.align(aligner, record)
            # Pick alignments with the highest score
            if length(alns) == 0
                @error "$name does **NOT** align to the genome uniquely or at all"
                push!(result, false)
                continue
            end
            max_score = maximum(map(x->(x.score), alns))
            alns = filter(x->x.score == max_score, alns)
            descriptions = map(x->description(x, aligner), alns)
            match = filter(x->occursin(chromosome_name, x), descriptions)
            nomatch = filter(x->!occursin(chromosome_name, x), descriptions)
            if length(match) > 0
                push!(result, true)
            else
                bad_chromosome = first(nomatch)
                push!(discard, (name, bad_chromosome))
                push!(result, false)
            end
        end
        for (name, chr) in discard
            @info "Discarded $name matching $chr"
        end
        return result
    end

    export bwa_sequences
end