"""
$(TYPEDEF)

Store the values reported by one cable datasheet without imposing a fixed
catalogue schema.

Field names and values are retained in the supplied `NamedTuple`. Property
access delegates to that record, so `info.resistance` returns the value stored
under `:resistance`.

$(TYPEDFIELDS)
"""
struct DatasheetInfo{D <: NamedTuple}
    "Named datasheet values in their stated engineering units."
    data::D

    function DatasheetInfo(data::D) where {D <: NamedTuple}
        :data in keys(data) && throw(ArgumentError(
            "datasheet field :data is reserved"
        ))
        return new{D}(data)
    end
end

"""
$(TYPEDSIGNATURES)

Construct cable datasheet information from named values.

# Keywords

- `values...`: Datasheet fields. Names and value types are retained exactly.

# Returns

- A [`DatasheetInfo`](@ref) record.
"""
DatasheetInfo(; values...) = DatasheetInfo((; values...))

Base.propertynames(info::DatasheetInfo, private::Bool = false) =
    private ? (:data, keys(getfield(info, :data))...) : keys(getfield(info, :data))

function Base.getproperty(info::DatasheetInfo, name::Symbol)
    name === :data && return getfield(info, :data)
    return getproperty(getfield(info, :data), name)
end

Base.keys(info::DatasheetInfo) = keys(getfield(info, :data))
Base.pairs(info::DatasheetInfo) = pairs(getfield(info, :data))
Base.length(info::DatasheetInfo) = length(getfield(info, :data))
Base.getindex(info::DatasheetInfo, name::Symbol) = getindex(getfield(info, :data), name)
