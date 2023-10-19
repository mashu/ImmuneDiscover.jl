module trim
    include("profile.jl")
    using .profile
    using Folds
    using DataFrames
    using ProgressMeter

    """
        find_maxprob_pos(read, prof)

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
        trim_profiles(motifs, len)

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
        trim_read(read, prof_start, prof_stop; min_length=0)

    Search read and return start and stop position that maximize likelihood `P(data|profile)`
    """
    function trim_read(read, prof_start, prof_stop; min_length=0)
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
        round_down_to_second_digit(x::Float64)
    
    Round down to second digit
    """
    function round_down_to_second_digit(x::Float64)
        return floor(x * 1e2) / 1e2
    end

    """
        find_segments(genomic_sequence, prof_start, prof_stop, minlen)

    Process table in parallel and perform trimming
    """
    function find_segments(genomic_sequence, prof_start, prof_stop, minlen)
        p = Progress(length(genomic_sequence))
        len = size(prof_start, 2)
        region = Folds.map(genomic_sequence) do seq
            ProgressMeter.next!(p)
            trim.trim_read(seq, prof_start, prof_stop, min_length=minlen)
        end
        n = length(genomic_sequence)
        ok = zeros(Bool, n)
        fragments = Vector{String}()
        for i in 1:n
            (start, stop) = region[i]
            if (start+len < stop) & (stop-start >= minlen)
                ok[i] = true
                fragment = genomic_sequence[i][start:stop-1]
                push!(fragments, fragment)
            end
        end
        valid = round_down_to_second_digit(sum(ok)/length(region)*100)
        @info "Trimmed $(sum(ok)) out of $(length(region)) ($valid %) sequences"
        return ok, fragments, region
    end
end
