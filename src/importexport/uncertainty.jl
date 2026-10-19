# Scientific result records are distinct from executable computation checkpoints.
serialize_value(value, ::Val{:scientific}) = serialize_value(value)
# JSON number readers need not preserve UInt64 values above typemax(Int64).
# Scientific seeds must retain all 64 bits, including within formulation records.
serialize_value(value::UInt64, ::Val{:scientific}) =
    Dict("__type__"=>"UInt64", "value"=>string(value))
deserialize_extension(::Val{:UInt64},record) = parse(UInt64,record["value"])
function serialize_value(value::NamedTuple, ::Val{:scientific})
    return Dict("__type__"=>"NamedTuple", "names"=>string.(collect(keys(value))),
        "values"=>[serialize_value(item, Val(:scientific)) for item in values(value)])
end
function serialize_value(value::Tuple, ::Val{:scientific})
    return Dict("__type__"=>"Tuple", "values"=>[serialize_value(item, Val(:scientific))
                                                for item in value])
end
function serialize_value(value::AbstractArray, ::Val{:scientific})
    return Dict("__type__"=>"Array", "size"=>collect(size(value)),
        "values"=>[serialize_value(item, Val(:scientific)) for item in vec(value)])
end
serialize_value(::Missing) = Dict("__type__"=>"Missing")
function serialize_value(value::AbstractArray{T, N}) where {T, N}
    return Dict("__type__"=>"Array", "size"=>collect(size(value)),
        "values"=>map(serialize_value, vec(value)))
end

# The selectors that saved requests and quantities name.
const _PORTABLE_SELECTORS = (
        Engine.Z, Engine.Y, Engine.R, Engine.X, Engine.L, Engine.G, Engine.B, Engine.C,
        ModalAnalysis.gamma, ModalAnalysis.alpha, ModalAnalysis.beta, ModalAnalysis.velocity,
        ModalAnalysis.Zc, ModalAnalysis.Yc, ModalAnalysis.H,
        ModalAnalysis.Tv, ModalAnalysis.Ti,
        frequencies, UQ.statistics, UQ.samples, UQ.histograms, Statistics.mean, Statistics.std,
        Statistics.median, minimum, maximum, abs, angle, real, imag,
        LinearAlgebra.diag, Engine.absolute_error, Engine.relative_error)

function serialize_value(selector::Function)
    selector in _PORTABLE_SELECTORS || throw(ArgumentError(
        "no portable scientific selector codec for $selector"))
    return Dict("__type__"=>"Observable", "name"=>string(nameof(selector)))
end
function serialize_value(selector::Base.Fix2{typeof(Statistics.quantile)})
    Dict("__type__"=>"Quantile", "probability"=>serialize_value(selector.x))
end
function serialize_value(selector::Base.Fix2{F}) where {F<:Union{typeof(ModalAnalysis.Zc),
        typeof(ModalAnalysis.Yc),typeof(ModalAnalysis.H)}}
    selector in Commons.observables(ModalAnalysis.PropagationParameters) ||
        throw(ArgumentError("unsupported modal representation selector"))
    return Dict("__type__"=>"ModalRepresentation",
        "selector"=>string(nameof(selector.f)),
        "field"=>selector.f===ModalAnalysis.H ? string(selector.x.field) : nothing)
end

function deserialize_extension(::Val{:Observable}, record)
    selected=filter(selector -> string(nameof(selector)) == record["name"], _PORTABLE_SELECTORS)
    length(selected)==1 || throw(ArgumentError("unknown saved scientific selector"))
    return only(selected)
end
function deserialize_extension(::Val{:Quantile}, record)
    Base.Fix2(Statistics.quantile, deserialize_value(record["probability"]))
end
function deserialize_extension(::Val{:ModalRepresentation},record)
    name=record["selector"]
    field=get(record,"field",nothing)
    if name=="H" && field in ("voltage","current")
        return Base.Fix2(ModalAnalysis.H,(domain=Engine.PhaseDomain,field=Symbol(field)))
    elseif name in ("Zc","Yc") && field===nothing
        selector=name=="Zc" ? ModalAnalysis.Zc : ModalAnalysis.Yc
        return Base.Fix2(selector,(domain=Engine.PhaseDomain,))
    end
    throw(ArgumentError("unknown saved modal representation"))
end
