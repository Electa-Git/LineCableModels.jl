# Scientific result records are distinct from executable computation checkpoints.
serialize_value(value, ::Val{:scientific}) = serialize_value(value)
# JSON number readers need not preserve UInt64 values above typemax(Int64).
# Scientific seeds must retain all 64 bits, including within formulation records.
serialize_value(value::UInt64, ::Val{:scientific}) =
    Dict("__type__"=>"UInt64", "value"=>string(value))
deserialize_extension(::Val{:UInt64},record) = parse(UInt64,record["value"])
function serialize_value(value::FormulationOptions)
    return Dict("__type__"=>"FormulationOptions", "value"=>serialize_value(value.data, Val(:scientific)))
end
function serialize_value(value::ComputationOptions)
    return Dict("__type__"=>"ComputationOptions", "value"=>serialize_value(value.data, Val(:scientific)))
end
function serialize_value(value::ComputationDetails)
    return Dict("__type__"=>"ComputationDetails", "value"=>serialize_value(value.data, Val(:scientific)))
end
deserialize_extension(::Val{:FormulationOptions}, record) =
    FormulationOptions(deserialize_value(record["value"]))
deserialize_extension(::Val{:ComputationOptions}, record) =
    ComputationOptions(deserialize_value(record["value"]))
deserialize_extension(::Val{:ComputationDetails}, record) =
    ComputationDetails(deserialize_value(record["value"]))
function serialize_value(::Val{Value}, ::Val{:scientific}) where {Value}
    return Dict("__type__"=>"Val", "value"=>serialize_value(Value, Val(:scientific)))
end
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
        Statistics.median, minimum, maximum, abs, angle, LinearAlgebra.diag) || throw(ArgumentError(
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
        Statistics.median, minimum, maximum, abs, angle, LinearAlgebra.diag)
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
    retained=LineCableModels.details(value).data
    return Dict("__type__"=>"LineParameters", "Z"=>serialize_value(observe(value, Z)),
        "Y"=>serialize_value(observe(value, Y)), "frequencies"=>serialize_value(frequencies(value)),
        "basis"=>string(LineCableModels.basis(value)), "domain"=>string(nameof(Engine.domain(value))),
        "coordinates"=>serialize_value(get(retained, :coordinates, nothing)),
        "formulations"=>serialize_value(get(retained, :formulations, nothing), Val(:scientific)),
        "shunt_model"=>serialize_value(get(retained, :shunt_model, nothing),Val(:scientific)),
        "comparison_unsupported"=>serialize_value(get(retained, :comparison_unsupported, (;))),
        "gridpoint_description"=>serialize_value((; (key=>retained[key] for key in
            (:inputs,:gridpoint,:selections,:formulation_fields,:uncertainty,:uncertainty_descriptions,:modal) if haskey(retained,key))...),Val(:scientific)))
end
function deserialize_extension(::Val{:LineParameters}, record)
    record["domain"] in ("PhaseDomain","ModalDomain") ||
        throw(ArgumentError("unsupported saved result domain"))
    coordinates=deserialize_value(get(record, "coordinates", nothing))
    unsupported=deserialize_value(get(record, "comparison_unsupported", Dict()))
    detail=(comparison_unsupported = (; (Symbol(k)=>v for (k, v) in pairs(unsupported))...),)
    coordinates === nothing || (detail=merge(detail, (; coordinates)))
    formulations=deserialize_value(get(record, "formulations", nothing))
    formulations === nothing || (detail=merge(detail,
        NamedTuple{(:formulations,),Tuple{NamedTuple}}((formulations,))))
    shunt_model=deserialize_value(get(record,"shunt_model",nothing))
    shunt_model === nothing || (detail=merge(detail,
        NamedTuple{(:shunt_model,),Tuple{NamedTuple}}((shunt_model,))))
    retained_description=deserialize_value(get(record,"gridpoint_description",serialize_value((;),Val(:scientific))))
    detail=merge(detail,retained_description)
    # Reading a primary result binds its passive selections to their owners.
    # Capture current compact descriptions here, before any observation exists;
    # retained ObservedResult loading and plotting never reopen this path.
    if formulations !== nothing
        selected=deserialize_value(Val(:formulation),formulations)
        ismissing(selected) || (detail=merge(detail,Engine.completed_formulation(selected,formulations)))
    end
    return LineParameters(
        getfield(Engine,Symbol(record["domain"])), deserialize_value(record["Z"]), deserialize_value(record["Y"]),
        deserialize_value(record["frequencies"]); basis = Symbol(record["basis"]), details = Engine.completion_details(detail))
end

function serialize_value(value::Engine.CableConstants)
    return Dict("__type__"=>"CableConstants", "cores"=>serialize_value(value.cores),
        "R"=>serialize_value(value.R),"L"=>serialize_value(value.L),
        "C"=>serialize_value(value.C),"G"=>serialize_value(value.G),
        "frequency"=>serialize_value(value.frequency),
        "details"=>serialize_value(value.details.data,Val(:scientific)))
end
function deserialize_extension(::Val{:CableConstants},record)
    retained=deserialize_value(record["details"])
    return Engine.CableConstants(Symbol.(deserialize_value(record["cores"])),
        (deserialize_value(record[key]) for key in ("R","L","C","G","frequency"))...,
        Engine.completion_details(retained))
end

"""
$(TYPEDSIGNATURES)

Encode retained MC products or first-order results without solving a model.
The versioned record retains scientific data, not executable formulations.
Measurement-bearing MC and LEP results use the extension's shared-source codec.
"""
function serialize_value(value::Union{UQ.MonteCarloResult, UQ.LinearErrorResult})
    return serialize_value(value, map(serialize_value, value.values), nothing)
end

"""
$(TYPEDSIGNATURES)

Encode a UQ result envelope with already encoded core `points`. Optional
`sources` are shared Measurement source records supplied by the Measurements
extension; their point records retain signed sensitivities. This boundary keeps
empirical products, provenance and details under the scientific result codec.
"""
function serialize_value(value::Union{UQ.MonteCarloResult,UQ.LinearErrorResult},
        points::AbstractVector, sources)
    record=NamedTuple(value)
    formulation=record.formulation isa NamedTuple ? record.formulation :
                NamedTuple(record.formulation)
    retained = record.details.data
    portable_details = if isempty(retained)
        retained
    elseif value isa UQ.LinearErrorResult
        (points=map(point -> point.data, retained.points),)
    else
        merge(retained, (trials=map(trials -> map(trial -> trial.data, trials), retained.trials),))
    end
    return Dict(
        "__type__"=>value isa UQ.MonteCarloResult ? "MonteCarloResult" :
                    "LinearErrorResult",
        "version"=>2, "formulation"=>serialize_value(formulation, Val(:scientific)),
        "points"=>points, "sources"=>serialize_value(sources,Val(:scientific)), "details"=>serialize_value(
            portable_details, Val(:scientific)),
        "statistics"=>serialize_value(get(record, :statistics, nothing), Val(:scientific)),
        "samples"=>serialize_value(get(record, :samples, nothing), Val(:scientific)),
        "histograms"=>serialize_value(get(record, :histograms, nothing), Val(:scientific)),
        "root_seed"=>serialize_value(get(record, :root_seed, nothing),Val(:scientific)),
        "point_seeds"=>serialize_value(get(record, :point_seeds, nothing),Val(:scientific)),
        "trial_counts"=>get(record, :trial_counts, nothing))
end

function deserialize_extension(kind::Union{Val{:MonteCarloResult}, Val{:LinearErrorResult}}, record)
    record["version"] in (1,2) ||
        throw(ArgumentError("unsupported scientific UQ record version"))
    decoded=deserialize_value(record["formulation"])
    formulation=(; (Symbol(k)=>v for (k, v) in pairs(decoded))...)
    options=(; (Symbol(k)=>v for (k, v) in pairs(formulation.options))...)
    formulation=merge(formulation, (; options))
    points=if get(record,"sources",nothing) === nothing
        map(deserialize_value, record["points"])
    else
        Base.get_extension(LineCableModels,:LineCableModelsMeasurementsExt) === nothing &&
            throw(ArgumentError("restoring uncertainty-bearing results requires `using Measurements`"))
        deserialize_extension(Val(:MeasurementPoints),record)
    end
    details=deserialize_value(record["details"])
    retained=(; (Symbol(k)=>v for (k, v) in pairs(details))...)
    if !isempty(retained)
        if kind isa Val{:LinearErrorResult}
            records=map(retained.points,points) do detail,point
                point isa Union{Engine.CableConstants,LineParameters} ?
                    Engine.completion_details(detail) : ComputationDetails(detail)
            end
            retained=(points=records,)
        else
            records=map(retained.trials,points) do trials,point
                map(trials) do detail
                    point isa Union{Engine.CableConstants,LineParameters} ?
                        Engine.completion_details(detail) : ComputationDetails(detail)
                end
            end
            retained=merge(retained,(trials=records,))
        end
    end
    details=ComputationDetails(retained)
    kind isa Val{:LinearErrorResult} &&
        return UQ.LinearErrorResult(formulation, points, details)
    products=map(("statistics", "samples", "histograms")) do key
        values=deserialize_value(record[key])
        values === nothing && return nothing
        [(; (Symbol(k)=>v for (k, v) in pairs(point))...) for point in values]
    end
    root_seed=deserialize_value(record["root_seed"])
    point_seeds=deserialize_value(record["point_seeds"])
    root_seed isa Integer && !(root_seed isa Bool) &&
        all(seed -> seed isa Integer && !(seed isa Bool),point_seeds) ||
        throw(ArgumentError("scientific Monte Carlo seeds require exact integers; floating-point JSON seeds cannot preserve provenance"))
    if record["version"] == 1
        # Supported portable MC v1 records retained full empirical summaries
        # but only mean-valued cores. Restore their documented marginal result
        # through the UQ-owned materialization, never through native checkpoints.
        points=map(LineCableModels.materialize,points,first(products))
    end
    return UQ.MonteCarloResult(
        formulation, points, products..., UInt64(root_seed),
        UInt64.(point_seeds), Int.(record["trial_counts"]), details)
end
