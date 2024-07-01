module keyedsets
    export KeyedPair, KeyedSet

    struct KeyedPair
        key::String
        value::String
    end

    struct KeyedSet
        data::Dict{String, String}
    end

    # Constructor
    KeyedSet() = KeyedSet(Dict{String, String}())

    # New constructor for Vector of Tuples
    function KeyedSet(pairs::Vector{Tuple{String, String}})
        ks = KeyedSet()
        for pair in pairs
            push!(ks, pair)
        end
        return ks
    end

    # Add a key-value pair
    function Base.push!(ks::KeyedSet, pair::KeyedPair)
        if !haskey(ks.data, pair.key)
            ks.data[pair.key] = pair.value
        elseif ks.data[pair.key] != pair.value
            @warn "Key $(pair.key) with value $(pair.value) already exists in KeyedSet with value $(ks.data[pair.key])"
        else
            @info "Duplicate key $(pair.key) with value $(pair.value) already exists in KeyedSet"
        end
        return ks
    end

    # Add a tuple of strings
    Base.push!(ks::KeyedSet, pair::Tuple{String, String}) = push!(ks, KeyedPair(pair...))

    # Implement set operations
    Base.in(key::String, ks::KeyedSet) = haskey(ks.data, key)
    Base.length(ks::KeyedSet) = length(ks.data)
    Base.iterate(ks::KeyedSet, state...) = iterate(keys(ks.data), state...)

    # Set operations
    Base.union(ks1::KeyedSet, ks2::KeyedSet) = KeyedSet(merge(ks1.data, ks2.data))
    Base.intersect(ks1::KeyedSet, ks2::KeyedSet) = KeyedSet(Dict(k => ks1.data[k] for k in keys(ks1.data) if k in keys(ks2.data)))
    Base.setdiff(ks1::KeyedSet, ks2::KeyedSet) = KeyedSet(Dict(k => v for (k, v) in ks1.data if !(k in keys(ks2.data))))

    # Equality
    Base.:(==)(ks1::KeyedSet, ks2::KeyedSet) = keys(ks1.data) == keys(ks2.data)

    # Get index
    #Base.getindex(ks::KeyedSet, key::String) = get(ks.data, key, error("No key $key in KeyedSet"))
    function Base.getindex(ks::KeyedSet, key::String)
        return get(ks.data, key) do
            error("No key $key in KeyedSet")
        end
    end
    # Collect
    Base.collect(ks::KeyedSet) = [(k, ks.data[k]) for k in keys(ks.data)]

    # Show method
    Base.show(io::IO, ks::KeyedSet) = print(io, "KeyedSet(size=", length(ks.data), ")")
end