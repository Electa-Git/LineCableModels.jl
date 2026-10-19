"""
$(TYPEDEF)

Read a saved result with Julia's `Serialization`. The fields before `sources` are the state
that `Serialization` keeps while it reads. An extension for uncertain values gives each
uncertainty source of the file a new identity and records it in `sources`.

$(TYPEDFIELDS)
"""
mutable struct ResultReader{I <: IO} <: Serialization.AbstractSerializer
    "Stream of the saved result."
    io::I
    "Object counter of `Serialization`."
    counter::Int
    "Objects read so far, by position."
    table::IdDict{Any, Any}
    "Positions of the objects still being read."
    pending_refs::Vector{Int}
    "Types and methods that the file defines."
    known_object_data::Dict{UInt64, Any}
    "Format version that the file header gives."
    version::Int
    "New uncertainty source for each source of the file."
    sources::Dict{Any, Any}
    function ResultReader(io::I) where {I <: IO}
        return new{I}(io, 0, IdDict{Any, Any}(), Int[], Dict{UInt64, Any}(), 0,
            Dict{Any, Any}())
    end
end

"""
$(TYPEDSIGNATURES)

Save a computed result with Julia's `Serialization`. `import_data(path, :result)` reads it
back. The file is for the same code. Another version of LineCableModels or of Julia can fail
to read it, as can another version of a package whose values the result contains.

# Arguments

- `result`: a computed result or a space of results.
- `path`: destination file.

# Returns

- Absolute path of the written file.
"""
function save(result::Commons.AbstractProblemResult, path::AbstractString)
    destination = abspath(path)
    Serialization.serialize(destination, result)
    return destination
end

"""
$(TYPEDSIGNATURES)

Read a result that [`save`](@ref) wrote, with the same code. The uncertainty sources of the
loaded values are new ones, one for each source of the file: the values keep their
correlations with one another and are independent of the sources of this process.

# Arguments

- `path`: file that `save` wrote.
- `:result`: selects a saved result.

# Returns

- The saved result.
"""
function import_data(path::AbstractString, ::Val{:result})
    return open(io -> Serialization.deserialize(ResultReader(io)),
        path)::Commons.AbstractProblemResult
end

# The selector after the path states what the file contains.
import_data(path::AbstractString, selector::Symbol) = import_data(path, Val(selector))
