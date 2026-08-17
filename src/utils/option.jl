module Option
    # Present/absent values for optional I/O. Dispatch on `Absent` vs `Present` replaces
    # `Union{T,Nothing}` in domain code. ArgParse may still yield `nothing`; convert at
    # the handler with `optional`. Named to avoid clashing with `Base.Some`.

    export Absent, Present, absent, optional, unwrap, or_default

    struct Absent end
    const absent = Absent()

    struct Present{T}
        value::T
    end

    optional(::Absent) = absent
    optional(s::Present) = s
    optional(::Nothing) = absent
    optional(x::AbstractString) = isempty(strip(x)) ? absent : Present(String(x))
    optional(x) = Present(x)

    unwrap(s::Present) = s.value

    or_default(::Absent, default) = default
    or_default(s::Present, _) = s.value
end
