module KeyedSets
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
    function Base.union(ks1::KeyedSet, ks2::KeyedSet)
        result = KeyedSet(copy(ks1.data))
        for (k, v2) in ks2.data
            if haskey(result.data, k)
                v1 = result.data[k]
                if v1 != v2
                    @info "Key $k exists in both sets with different names. Using name '$(v1)' from the first set instead of '$(v2)' from the second set."
                end
            else
                result.data[k] = v2
            end
        end
        return result
    end

    function Base.intersect(ks1::KeyedSet, ks2::KeyedSet)
        result = KeyedSet()
        for (k, v1) in ks1.data
            if haskey(ks2.data, k)
                v2 = ks2.data[k]
                if v1 != v2
                    @info "Sequence with key $k is identical but has different names: $(v1) in set 1, $(v2) in set 2"
                end
                push!(result, (k, v1))
            end
        end
        return result
    end

    function Base.setdiff(ks1::KeyedSet, ks2::KeyedSet)
        result = KeyedSet()
        for (k1, v1) in ks1.data
            if !haskey(ks2.data, k1)
                push!(result, (k1, v1))
            end
        end
        return result
    end

    Base.collect(ks::KeyedSet) = [(k, v) for (k, v) in ks.data]
end
