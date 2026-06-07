module KeyedSets
    export KeyedPair, KeyedSet

    struct KeyedPair
        key::String
        value::String
    end

    struct KeyedSet
        data::Dict{String, String}
    end

    KeyedSet() = KeyedSet(Dict{String, String}())

    function KeyedSet(pairs::Vector{Tuple{String, String}})
        ks = KeyedSet()
        for pair in pairs
            push!(ks, pair)
        end
        return ks
    end

    function _register_sequence_name!(data::Dict{String,String}, sequence::String, name::String; strict::Bool=true)
        if haskey(data, sequence)
            existing = data[sequence]
            existing != name && throw(ArgumentError(
                "Identical sequence appears under different names: \"$existing\" and \"$name\""))
            strict && throw(ArgumentError("Duplicate sequence entry for \"$name\""))
            return nothing
        end
        data[sequence] = name
        return nothing
    end

    function Base.push!(ks::KeyedSet, pair::KeyedPair)
        _register_sequence_name!(ks.data, pair.key, pair.value)
        return ks
    end

    Base.push!(ks::KeyedSet, pair::Tuple{String, String}) = push!(ks, KeyedPair(pair...))

    Base.in(key::String, ks::KeyedSet) = haskey(ks.data, key)
    Base.length(ks::KeyedSet) = length(ks.data)
    Base.iterate(ks::KeyedSet, state...) = iterate(keys(ks.data), state...)
    Base.getindex(ks::KeyedSet, key::String) = ks.data[key]
    Base.:(==)(ks1::KeyedSet, ks2::KeyedSet) = ks1.data == ks2.data
    Base.show(io::IO, ks::KeyedSet) = print(io, "KeyedSet(size=$(length(ks)))")

    function Base.union(ks1::KeyedSet, ks2::KeyedSet)
        result = KeyedSet(copy(ks1.data))
        for (sequence, name) in ks2.data
            _register_sequence_name!(result.data, sequence, name; strict=false)
        end
        return result
    end

    function Base.intersect(ks1::KeyedSet, ks2::KeyedSet)
        result = KeyedSet()
        for (sequence, name1) in ks1.data
            haskey(ks2.data, sequence) || continue
            name2 = ks2.data[sequence]
            name1 != name2 && throw(ArgumentError(
                "Identical sequence has different names: \"$name1\" vs \"$name2\""))
            push!(result, (sequence, name1))
        end
        return result
    end

    function Base.setdiff(ks1::KeyedSet, ks2::KeyedSet)
        result = KeyedSet()
        for (sequence, name) in ks1.data
            haskey(ks2.data, sequence) || push!(result, (sequence, name))
        end
        return result
    end

    Base.collect(ks::KeyedSet) = [(k, v) for (k, v) in ks.data]
end
