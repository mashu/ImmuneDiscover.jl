module Spans
    # Exact substring locations and clipped DNA windows. Shared by exact search, heptamer
    # search, and BLAST affix collection so each caller does not reimplement findfirst-once.

    export each_exact_span, flank_slice

    """
        each_exact_span(needle, haystack) -> Vector{UnitRange{Int}}

    Non-overlapping exact occurrences of `needle` in `haystack`, left to right.
    Empty needles are ignored (they would otherwise match at every position).
    """
    function each_exact_span(needle::AbstractString, haystack::AbstractString)
        isempty(needle) && return UnitRange{Int}[]
        spans = UnitRange{Int}[]
        start = firstindex(haystack)
        stop = lastindex(haystack)
        while start <= stop
            m = findnext(needle, haystack, start)
            m === nothing && return spans
            push!(spans, m)
            nxt = nextind(haystack, last(m))
            nxt > stop && return spans
            start = nxt
        end
        return spans
    end

    """
        flank_slice(seq, start, stop) -> String

    Inclusive window clipped to `seq`. Returns "" when the window does not overlap the sequence,
    instead of clamping both ends to 1 (which used to produce overlapping garbage RSS windows).
    """
    function flank_slice(seq::AbstractString, start::Int, stop::Int)
        n = lastindex(seq)
        lo = max(firstindex(seq), start)
        hi = min(n, stop)
        lo > hi && return ""
        return String(seq[lo:hi])
    end
end
