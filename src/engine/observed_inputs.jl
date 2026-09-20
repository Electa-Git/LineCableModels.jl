# This projection is restricted to physical declarations owned by the package.
# It retains their named values, not constructed problems or geometry caches.
_input_record(value::Union{Number,Symbol,AbstractString,Nothing,Missing,Type}) = value
_input_record(value::NamedTuple) = map(_input_record, value)
_input_record(value::Tuple) = map(_input_record, value)
_input_record(value::AbstractArray) = map(_input_record, value)
_input_record(value::Pair) = _input_record(first(value)) => _input_record(last(value))
_input_record(value::AbstractDict) = Dict(_input_record(k) => _input_record(v) for (k,v) in value)

function _input_record(value::Union{DataModel.AbstractCablePart,DataModel.AbstractShape,
        DataModel.Pose2,DataModel.Shell,DataModel.Ring,DataModel.Polar,DataModel.Fill,
        DataModel.Lattice,DataModel.Helix,DataModel.LayRatio,DataModel.Pitch,
        DataModel.LayAngle,DataModel.FillFactor,DataModel.AssemblyMember,
        Materials.AbstractMaterial,Earth.AbstractEarthModel,Earth.AbstractEarthLayer,
        Earth.AbstractEarthMaterial})
    fields = fieldnames(typeof(value))
    return merge((kind=nameof(typeof(value)),),
        NamedTuple{fields}(map(name -> _input_record(getfield(value,name)), fields)))
end

_input_record(design::CableDesign) = (cable_id=design.cable_id,
    origin=_input_record(design.origin),
    terminal_order=copy(design.terminal_order))

function _input_record(system::LineCableSystem)
    return (system_id=system.system_id, line_length=system.line_length,
        designs=_input_record(system.designs), positions=_input_record(system.positions),
        input_positions=_input_record(system.input_positions), clearances=copy(system.clearances),
        connections=_input_record(system.connections), environment=_input_record(system.environment),
        terminal_order=copy(system.terminal_order), connection_order=copy(system.connection_order))
end

"""
$(TYPEDSIGNATURES)

Capture physical declarations from a completed package-owned problem. Store the
returned record in result details at completion, once per physical point; later
observation reads that record without evaluating problem builders.
"""
function completed_inputs(problem::LineParametersProblem)
    return NamedTuple{(:system,:temperature,:earth_props,:frequencies),
        Tuple{NamedTuple,typeof(problem.temperature),NamedTuple,typeof(problem.frequencies)}}((
        _input_record(problem.system),problem.temperature,_input_record(problem.earth_props),copy(problem.frequencies)))
end
function completed_inputs(problem::CableConstantsProblem)
    return NamedTuple{(:design,:temperature,:frequency),
        Tuple{NamedTuple,typeof(problem.temperature),typeof(problem.frequency)}}((
        _input_record(problem.design),problem.temperature,problem.frequency))
end

# Capture human descriptions while the actual selected owner objects exist.
# Consumers read these strings and never reconstruct a live formulation.
"""
$(TYPEDSIGNATURES)

Capture the actual formulation selections, controls, and owned human descriptions
when a result completes. The returned values contain no live formulation objects.
"""
function completed_formulation(formulation)
    return NamedTuple{(:formulations,:selections,:formulation_labels),
        Tuple{NamedTuple,NamedTuple,NamedTuple}}((NamedTuple(formulation),
        (Z=formula_id(formulation,Z),Y=formula_id(formulation,Y)),
        (all=description(formulation,nothing),
            Z=description(formulation,Z),Y=description(formulation,Y))))
end

# Description contents and runtime geometry do not specialize the result type.
"""
$(TYPEDSIGNATURES)

Store captured completion records without specializing the result type on runtime
geometry or nested description contents. Persistence uses this same operation.
"""
function completion_details(record::NamedTuple{names,T}) where {names,T}
    types=map(value -> value isa NamedTuple ? NamedTuple : typeof(value),values(record))
    return ComputationDetails(NamedTuple{names,Tuple{types...}}(values(record)))
end

function Grammar.observation_gridpoint(source::Union{LineParameters,CableConstants})
    retained=details(source).data
    inputs=get(retained,:inputs,nothing)
    return Grammar.detach((id=get(retained,:gridpoint,nothing), inputs,
        formulations=get(retained,:formulations,nothing),formulation_labels=get(retained,:formulation_labels,(;)),
        coordinates=get(retained,:coordinates,source isa CableConstants ? source.cores : nothing),
        uncertainty=get(retained,:uncertainty,nothing),
        transformation=get(retained,:modal,nothing),
        missing_reason=inputs===nothing ? :physical_inputs_not_supplied : nothing))
end

# Updating an association reuses numerical storage; it does not recompute or
# copy a physical problem. External result owners can supply their own method.
"""
$(TYPEDSIGNATURES)

Associate a completed result with its original gridpoint identity and additional
captured fields. Built-in results retain their numerical arrays and replace only
completion details. External result owners may extend this completion operation.
"""
retain_gridpoint(source, id; fields=(;)) = source
function retain_gridpoint(source::LineParameters, id; fields=(;))
    retained=completion_details(merge(details(source).data,fields,(gridpoint=id,)))
    return LineParameters(source.domain,source.Z,source.Y,source.f,retained)
end
function retain_gridpoint(source::CableConstants, id; fields=(;))
    retained=completion_details(merge(details(source).data,fields,(gridpoint=id,)))
    return CableConstants(source.cores,source.R,source.L,source.C,source.G,source.frequency,retained)
end

function _resolution_length(source::AbstractCoreResult)
    inputs=get(details(source).data,:inputs,nothing)
    inputs===nothing && return nothing
    system=get(inputs,:system,nothing)
    return system===nothing ? nothing : system.line_length
end
_resolution_length(source) = nothing
