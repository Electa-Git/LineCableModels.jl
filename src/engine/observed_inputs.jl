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
    return merge((kind=nameof(typeof(value)),field_descriptions=Commons.input_fields(typeof(value))),
        NamedTuple{fields}(map(name -> _input_record(getfield(value,name)), fields)))
end

_input_record(design::CableDesign) = (cable_id=design.cable_id,
    origin=_input_record(design.origin),
    terminal_order=copy(design.terminal_order))

function _input_record(system::LineCableSystem)
    return (system_id=system.system_id, field_descriptions=Commons.input_fields(typeof(system)), line_length=system.line_length,
        designs=_input_record(system.designs), positions=_input_record(system.positions),
        input_positions=_input_record(system.input_positions), clearances=copy(system.clearances),
        connections=_input_record(system.connections), environment=_input_record(system.environment),
        terminal_order=copy(system.terminal_order), connection_order=copy(system.connection_order))
end

"""
$(TYPEDSIGNATURES)

Capture physical declarations from a completed package-owned problem. Store the
returned record in result details at completion, once per physical point. Later
observation reads that record without evaluating problem builders.
"""
function completed_inputs(problem::LineParametersProblem)
    Record=NamedTuple{(:system,:temperature,:earth_props,:frequencies,:field_descriptions),
        Tuple{NamedTuple,typeof(problem.temperature),NamedTuple,typeof(problem.frequencies),NamedTuple}}
    return Record((_input_record(problem.system),problem.temperature,_input_record(problem.earth_props),
        copy(problem.frequencies),(temperature=(name="temperature",unit="°C"),)))
end
function completed_inputs(problem::CableConstantsProblem)
    Record=NamedTuple{(:design,:temperature,:frequency,:field_descriptions),
        Tuple{NamedTuple,typeof(problem.temperature),typeof(problem.frequency),NamedTuple}}
    return Record((_input_record(problem.design),problem.temperature,problem.frequency,
        (temperature=(name="temperature",unit="°C"),frequency=(name="frequency",unit="Hz"))))
end

# Capture owner-scoped meanings and names while the selected objects exist.
# Common-field omission later compares selections, never punctuation.
"""
$(TYPEDSIGNATURES)

Capture actual formulation selections, controls, and structured description
fields when a result completes.
"""
function completed_formulation(formulation, declaration::NamedTuple=NamedTuple(formulation))
    function control_fields(owner, controls, path=())
        captured=NamedTuple[]
        for (key,value) in pairs(controls)
            route=(path...,key)
            if value isa NamedTuple && key !== :equivalent_earth
                append!(captured,control_fields(owner,value,route))
            else
                text=key === :equivalent_earth ?
                     description(FormulaDefinition,Val(key),value;compact=true) :
                    description(owner,Val(key),value;compact=true)
                push!(captured,(scope=route,value=Commons.detach(value),text))
            end
        end
        return captured
    end
    function fields(quantity)
        selections=pairs(formulation;quantity)
        map(collect(selections)) do (scope,selected)
            owner,route=scope
            name=isempty(route) ? "" : description(owner,Val(first(route));compact=true)
            length(route)>1 && (name *= "("*join(string.(Base.tail(route)),",")*")")
            controls=selected isa Pair ? Commons.detach(last(selected)) : (;)
            # The description compares applied controls, including defaults.
            # Keep the original selection identity used by scientific grouping.
            effective=get(declaration,:methods,nothing)
            for key in route
                effective=effective isa NamedTuple ? get(effective,key,nothing) : nothing
            end
            applied=isempty(route) || !(effective isa NamedTuple) ? controls : merge(controls,
                (; (key=>effective[key] for key in (:parameters,:options,:equivalent_earth)
                    if haskey(effective,key) && effective[key]!==nothing)...))
            leaf=selected isa Pair ? first(selected) : selected
            control_owner=isempty(route) ? owner : typeof(leaf)
            # Routes are the scientific slots declared by pairs(owner), shared
            # by backends using the same constitutive meaning.
            (meaning=route,
                selection=(identifier=formula_id(selected),controls),name,
                value=isempty(route) ? description(selected;compact=true) : description(owner,selected;compact=true),
                summary=isempty(route) ? description(owner,formulation;quantity,compact=true) : nothing,
                control_fields=control_fields(control_owner,applied))
        end
    end
    return NamedTuple{(:formulations,:selections,:formulation_fields),
        Tuple{NamedTuple,NamedTuple,NamedTuple}}((declaration,
        (Z=formula_id(formulation,Z),Y=formula_id(formulation,Y)),
        (all=fields(nothing),Z=fields(Z),Y=fields(Y))))
end

"""
$(TYPEDSIGNATURES)

Complete `formulation` as the line-parameter computation of `workspace` applied it. Each
earth formula records the equivalent earth that its calculation in the plan consumed. A
reduction that the plan resolved through `:default` appears as that formula's
`equivalent_earth`, as an explicit one does. The requested selections remain as declared.
"""
function completed_formulation(formulation::LineParametersFormulation,
        workspace::LineParametersWorkspace)
    declaration = NamedTuple(formulation)
    length(workspace.input.earth.layers) > 2 ||
        return completed_formulation(formulation, declaration)
    slots = (:earth_impedance, :earth_admittance)
    # Each formula without an explicit `equivalent_earth` that the plan reduced, with the
    # reduction its calculation consumed.
    reduced = Pair[]
    for calculation in workspace.plan.earth.calculations
        calculation.earth isa EarthModel && continue
        for entry in (calculation.impedance, calculation.admittance)
            entry === nothing || entry.formula.equivalent_earth !== nothing ||
                push!(reduced, entry.formula => calculation.earth)
        end
    end
    function record(selected, retained)
        selected isa NamedTuple && return map(record, selected, retained)
        index = findfirst(entry -> first(entry) === selected, reduced)
        index === nothing && return retained
        return merge(retained, (equivalent_earth = NamedTuple(last(reduced[index])),))
    end
    recorded = merge(declaration.methods,
        map(record, formulation.methods[slots], declaration.methods[slots]))
    return completed_formulation(formulation,
        typeof(declaration)((declaration.backend, declaration.requested, recorded, declaration.options)))
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

function Commons.observation_gridpoint(source::Union{LineParameters,CableConstants})
    retained=details(source).data
    inputs=get(retained,:inputs,nothing)
    return Commons.detach((id=get(retained,:gridpoint,nothing), inputs,
        source_gridpoint=get(retained,:source_gridpoint,nothing),
        formulations=get(retained,:formulations,nothing),formulation_fields=get(retained,:formulation_fields,(;)),
        coordinates=get(retained,:coordinates,source isa CableConstants ? source.cores : nothing),
        uncertainty=get(retained,:uncertainty,nothing),
        uncertainty_descriptions=get(retained,:uncertainty_descriptions,nothing),
        transformation=get(retained,:modal,nothing),
        missing_reason=inputs===nothing ? :physical_inputs_not_supplied : nothing))
end

# Updating an association reuses numerical storage. It does not recompute or
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

function line_length(source::AbstractCoreResult)
    inputs=get(details(source).data,:inputs,nothing)
    inputs===nothing && return nothing
    system=get(inputs,:system,nothing)
    return system===nothing ? nothing : system.line_length
end
line_length(source) = nothing
