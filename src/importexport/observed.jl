# Each archive owns one uncertainty-source context spanning all points and the
# reference. No process-global source registry or behavior-generation tag exists.
_observed_encoding() = (indices=Dict{Any,Int}(),sources=Any[])
"""
$(TYPEDSIGNATURES)

Encode detached values using one archive-wide uncertainty-source context.
Extensions register independent sources in `context.sources` and reuse their
indices in `context.indices`, preserving dependencies across all observations.
"""
encode_observation(value,context) = serialize_value(value,Val(:scientific))
function encode_observation(value::NamedTuple,context)
    Dict("__type__"=>"NamedTuple","names"=>string.(collect(keys(value))),
        "values"=>[encode_observation(item,context) for item in values(value)])
end
encode_observation(value::Tuple,context) = Dict("__type__"=>"Tuple",
    "values"=>[encode_observation(item,context) for item in value])
encode_observation(value::AbstractArray,context) = Dict("__type__"=>"Array","size"=>collect(size(value)),
    "values"=>[encode_observation(item,context) for item in vec(value)])
encode_observation(value::AbstractDict,context) = Dict("__type__"=>"Dictionary",
    "entries"=>[encode_observation((key,item),context) for (key,item) in value])
encode_observation(value::Complex,context) = Dict("__type__"=>"Complex",
    "re"=>encode_observation(real(value),context),"im"=>encode_observation(imag(value),context))
encode_observation(value::Grammar.ObservedResult,context) = Dict("__type__"=>"ObservedResult",
    "fields"=>encode_observation((value.gridpoint,value.quantities,value.errors,value.timings),context))

function serialize_value(value::Union{Grammar.ObservedResult,AbstractVector{<:Grammar.ObservedResult}})
    return serialize_value(value,Val(:observed))
end
function serialize_value(value,::Val{:observed})
    context=_observed_encoding()
    payload=encode_observation(value,context)
    return Dict("__type__"=>"ObservedArchive","payload"=>payload,
        "sources"=>serialize_value(context.sources,Val(:scientific)))
end
function serialize_value(artifact::ReportBuilder.ReportArtifact)
    return serialize_value((observed=artifact.observed,reference=artifact.reference),Val(:observed))
end

_decode_observed(value,context) = deserialize_value(value)
function _decode_observed(value::AbstractDict,context)
    marker=get(value,"__type__",nothing)
    if marker=="ObservedResult"
        return Grammar.ObservedResult(_decode_observed(value["fields"],context)...)
    elseif marker=="ObservedMeasurement"
        return decode_observation_measurement(value,context)
    elseif marker=="NamedTuple"
        return NamedTuple{Tuple(Symbol.(value["names"]))}(Tuple(_decode_observed(item,context) for item in value["values"]))
    elseif marker=="Tuple"
        return Tuple(_decode_observed(item,context) for item in value["values"])
    elseif marker=="Array"
        elements=map(item -> _decode_observed(item,context),value["values"])
        if !isempty(elements)
            T=typeof(first(elements))
            if all(item -> item isa T,elements)
                elements=collect(T,elements)
            elseif all(item -> item isa Number || ismissing(item),elements)
                scalar_type=foldl((left,right) -> Union{left,right},typeof.(elements))
                elements=collect(scalar_type,elements)
            end
        end
        return reshape(elements,Tuple(Int.(value["size"])))
    elseif marker=="Dictionary"
        return Dict(_decode_observed(entry,context) for entry in value["entries"])
    elseif marker=="Complex"
        return complex(_decode_observed(value["re"],context),_decode_observed(value["im"],context))
    end
    return deserialize_value(value)
end

"""Restore one uncertain scalar from its archived sensitivities and shared sources."""
function decode_observation_measurement end
function deserialize_extension(::Val{:ObservedArchive},record)
    sources=deserialize_value(record["sources"])
    context=isempty(sources) ? () : observation_sources(sources,Val(:measurements))
    return _decode_observed(record["payload"],context)
end
"""Restore the archive's independent uncertainty sources once, before its values."""
function observation_sources(records,::Val{:measurements})
    throw(ArgumentError("restoring uncertain observations requires using Measurements"))
end

serialize_value(value::UUIDs.UUID) = Dict("__type__"=>"UUID","value"=>string(value))
deserialize_extension(::Val{:UUID},record) = UUIDs.UUID(record["value"])
serialize_value(::Colon) = Dict("__type__"=>"Colon")
deserialize_extension(::Val{:Colon},record) = Colon()
serialize_value(value::Units.Quantity{Q}) where {Q} = Dict("__type__"=>"Quantity","name"=>string(Q))
deserialize_extension(::Val{:Quantity},record) = Units.Quantity{Symbol(record["name"])}()
serialize_value(value::Units.Unit) = Dict("__type__"=>"Unit","name"=>string(value.name),"prefix"=>string(value.prefix))
deserialize_extension(::Val{:Unit},record) = Units.Unit(Symbol(record["name"]),Symbol(record["prefix"]))
serialize_value(value::Units.UnitExpr) = Dict("__type__"=>"UnitExpr",
    "numerator"=>serialize_value(value.numerator,Val(:scientific)),"denominator"=>serialize_value(value.denominator,Val(:scientific)))
deserialize_extension(::Val{:UnitExpr},record) = Units.UnitExpr(deserialize_value(record["numerator"]),deserialize_value(record["denominator"]))
function serialize_value(value::Type)
    owner=parentmodule(value)
    path=Base.fullname(owner)
    first(path) in (:LineCableModels,:Base,:Core) || throw(ArgumentError("no portable type identity for $value"))
    getfield(owner,nameof(value))===value || throw(ArgumentError("parametric runtime types are not portable formulation identities"))
    return Dict("__type__"=>"ScientificType","path"=>string.(path),"name"=>string(nameof(value)))
end
function deserialize_extension(::Val{:ScientificType},record)
    path=Symbol.(record["path"])
    root=first(path)
    owner=root===:LineCableModels ? LineCableModels : root===:Base ? Base : root===:Core ? Core :
        throw(ArgumentError("unknown scientific type owner"))
    for name in path[2:end]
        owner=getfield(owner,name)
        owner isa Module || throw(ArgumentError("scientific type owner must be a module"))
    end
    value=getfield(owner,Symbol(record["name"]))
    value isa Type || throw(ArgumentError("saved scientific identity does not name a type"))
    return value
end

"""
$(TYPEDSIGNATURES)

Save current observed data as JSON or a native Julia archive. One source context
preserves shared uncertainties across quantities, points, and a report reference;
BigFloat values retain their precision. Report tables and figures are rebuilt
from these observations after loading, rather than persisted as competing data.
"""
function save(value::Union{Grammar.ObservedResult,AbstractVector{<:Grammar.ObservedResult},ReportBuilder.ReportArtifact},path::AbstractString)
    encoded=serialize_value(value)
    extension=lowercase(splitext(path)[2])
    extension in (".json",".jls") || throw(ArgumentError("observed archives require .json or .jls"))
    destination=abspath(path)
    temporary=tempname(dirname(destination))
    try
        open(temporary,"w") do io
            extension==".json" ? JSON3.write(io,encoded) : Serialization.serialize(io,encoded)
        end
        mv(temporary,destination;force=true)
    finally
        isfile(temporary) && rm(temporary)
    end
    return destination
end

function import_data(::Val{:observed},path::AbstractString)
    extension=lowercase(splitext(path)[2])
    encoded=extension==".json" ? JSON3.read(read(path,String)) : extension==".jls" ?
        open(Serialization.deserialize,path) : throw(ArgumentError("observed archives require .json or .jls"))
    get(encoded,"__type__",nothing)=="ObservedArchive" || throw(ArgumentError("file is not a current observed archive"))
    return deserialize_value(encoded)
end
