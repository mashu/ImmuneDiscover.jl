module trim
    include("profile.jl")
    using .profile
    using Folds
    using DataFrames
    using ProgressMeter

    """
    Find the most probable position in a read given a profile
    """
    function find_maxprob_pos(read, prof)
        prof_length = size(prof, 2)
        prob_vec = []
        ret = 1
        for i in 1:(length(read)-prof_length)
            prob = profile.motif_prob(read[i:i+prof_length-1], prof)
            push!(prob_vec, prob)
        end
        prob_max = maximum(prob_vec)
        for i in eachindex(prob_vec)
            if prob_vec[i] == prob_max
                ret = i
                break
            end
        end
        return ret
    end

    """
    Create start and stop trim profiles
    """
    function trim_profiles(motifs, len)
        start_motifs = [r[1:len] for r in motifs]
        stop_motifs = [r[end-len+1:end] for r in motifs]
        prof_start = profile.motif_profile(start_motifs)
        prof_stop = profile.motif_profile(stop_motifs)
        return prof_start, prof_stop
    end

    """
    Search read and return start and stop position that maximize likelihood `P(data|profile)`
    """
    function trim_read(read, prof_start, prof_stop)
        len = size(prof_start, 2)
        start = trim.find_maxprob_pos(read, prof_start)
        residual = read[start+len:end]
        if length(residual) > len
            stop = trim.find_maxprob_pos(residual, prof_stop)-1+start+len+len
            if stop <= start
                stop = -1
            end
        else
            stop = -1
        end
        return start, stop
    end

    """
    Process table in parallel and perform trimming
    """
    function find_segments(table, prof_start, prof_stop)
        len = size(prof_start, 2)
        p = Progress(nrow(table))
        Folds.map(eachrow(table)) do row
            sequence = row.genomic_sequence
            ProgressMeter.next!(p)
            start, stop = trim.trim_read(sequence, prof_start, prof_stop)
            if stop > 0
                prefix = sequence[start-6:start]
                center = sequence[start:stop-1]
                suffix = sequence[stop:stop-1+7]
                (prefix, center, suffix)
            else
                ("", "", "")
            end
        end
    end
end
