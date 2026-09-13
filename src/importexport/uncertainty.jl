# Scientific result records are distinct from executable computation checkpoints.
serialize_value(value, ::Val{:scientific}) = serialize_value(value)
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

function serialize_value(selector::Function)
    selector in (
        Engine.Z, Engine.Y, Engine.R, Engine.X, Engine.L, Engine.G, Engine.B, Engine.C,
        frequencies, UQ.statistics, UQ.samples, UQ.histograms, Statistics.mean, Statistics.std,
        Statistics.median, minimum, maximum, abs, angle) || throw(ArgumentError(
        "no portable scientific selector codec for $selector"))
    return Dict("__type__"=>"Observable", "name"=>string(nameof(selector)))
end
function serialize_value(selector::Base.Fix2{typeof(Statistics.quantile)})
    Dict("__type__"=>"Quantile", "probability"=>serialize_value(selector.x))
end

function deserialize_extension(::Val{:Observable}, record)
    selectors=(
        Engine.Z, Engine.Y, Engine.R, Engine.X, Engine.L, Engine.G, Engine.B, Engine.C,
        frequencies, UQ.statistics, UQ.samples, UQ.histograms, Statistics.mean, Statistics.std,
        Statistics.median, minimum, maximum, abs, angle)
    selected=filter(selector -> string(nameof(selector)) == record["name"], selectors)
    length(selected)==1 || throw(ArgumentError("unknown saved scientific selector"))
    return only(selected)
end
function deserialize_extension(::Val{:Quantile}, record)
    Base.Fix2(Statistics.quantile, deserialize_value(record["probability"]))
end

function serialize_value(value::UQ.SampleSummary)
    return Dict("__type__"=>"SampleSummary", "values"=>serialize_value(Tuple(NamedTuple(value))))
end
function deserialize_extension(::Val{:SampleSummary}, record)
    UQ.SampleSummary(deserialize_value(record["values"])...)
end

function serialize_value(value::UQ.HistogramDensity)
    record=NamedTuple(value)
    return Dict("__type__"=>"HistogramDensity", "edges"=>serialize_value(record.edges),
        "density"=>serialize_value(record.density))
end
function deserialize_extension(::Val{:HistogramDensity}, record)
    UQ.HistogramDensity(deserialize_value(record["edges"]), deserialize_value(record["density"]))
end

function serialize_value(value::LineParameters)
    retained=LineCableModels.details(value)
    return Dict("__type__"=>"LineParameters", "Z"=>serialize_value(observe(value, Z)),
        "Y"=>serialize_value(observe(value, Y)), "frequencies"=>serialize_value(frequencies(value)),
        "basis"=>string(LineCableModels.basis(value)), "domain"=>string(nameof(Engine.domain(value))),
        "coordinates"=>serialize_value(get(retained, :coordinates, nothing)),
        "comparison_unsupported"=>serialize_value(get(retained, :comparison_unsupported, (;))))
end
function deserialize_extension(::Val{:LineParameters}, record)
    record["domain"] == "PhaseDomain" ||
        throw(ArgumentError("unsupported saved result domain"))
    coordinates=deserialize_value(get(record, "coordinates", nothing))
    unsupported=deserialize_value(get(record, "comparison_unsupported", Dict()))
    detail=(comparison_unsupported = (; (Symbol(k)=>v for (k, v) in pairs(unsupported))...),)
    coordinates === nothing || (detail=merge(detail, (; coordinates)))
    return LineParameters(
        Engine.PhaseDomain, deserialize_value(record["Z"]), deserialize_value(record["Y"]),
        deserialize_value(record["frequencies"]); basis = Symbol(record["basis"]), details = detail)
end

"""
$(TYPEDSIGNATURES)

Encode retained MC products or first-order results without solving a model.
The versioned record retains scientific data, not executable formulations.
Measurement-bearing LEP results use the extension's shared-source codec.
"""
function serialize_value(value::Union{UQ.MonteCarloResult, UQ.LinearErrorResult})
    record=NamedTuple(value)
    formulation=record.formulation isa NamedTuple ? record.formulation :
                NamedTuple(record.formulation)
    return Dict(
        "__type__"=>value isa UQ.MonteCarloResult ? "MonteCarloResult" :
                    "LinearErrorResult",
        "version"=>1, "formulation"=>serialize_value(formulation, Val(:scientific)),
        "points"=>map(serialize_value, record.values), "details"=>serialize_value(
            record.details, Val(:scientific)),
        "statistics"=>serialize_value(get(record, :statistics, nothing), Val(:scientific)),
        "samples"=>serialize_value(get(record, :samples, nothing), Val(:scientific)),
        "histograms"=>serialize_value(get(record, :histograms, nothing), Val(:scientific)),
        "root_seed"=>get(record, :root_seed, nothing), "point_seeds"=>get(record, :point_seeds, nothing),
        "trial_counts"=>get(record, :trial_counts, nothing))
end

function deserialize_extension(kind::Union{Val{:MonteCarloResult}, Val{:LinearErrorResult}}, record)
    record["version"] == 1 ||
        throw(ArgumentError("unsupported scientific UQ record version"))
    decoded=deserialize_value(record["formulation"])
    formulation=(; (Symbol(k)=>v for (k, v) in pairs(decoded))...)
    options=(; (Symbol(k)=>v for (k, v) in pairs(formulation.options))...)
    formulation=merge(formulation, (; options))
    points=map(deserialize_value, record["points"])
    details=deserialize_value(record["details"])
    details=(; (Symbol(k)=>v for (k, v) in pairs(details))...)
    kind isa Val{:LinearErrorResult} &&
        return UQ.LinearErrorResult(formulation, points, details)
    products=map(("statistics", "samples", "histograms")) do key
        values=deserialize_value(record[key])
        values === nothing && return nothing
        [(; (Symbol(k)=>v for (k, v) in pairs(point))...) for point in values]
    end
    return UQ.MonteCarloResult(
        formulation, points, products..., UInt64(record["root_seed"]),
        UInt64.(record["point_seeds"]), Int.(record["trial_counts"]), details)
end
