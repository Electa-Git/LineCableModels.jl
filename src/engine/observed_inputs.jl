# This projection is restricted to physical declarations owned by the package.
# It retains their named values, not constructed problems or geometry caches.
_input_record(value::Union{Number,Symbol,AbstractString,Nothing,Missing,Type}) = value
_input_record(value::NamedTuple) = map(_input_record, value)
_input_record(value::Tuple) = map(_input_record, value)
_input_record(value::AbstractArray) = map(_input_record, value)
_input_record(value::Pair) = _input_record(first(value)) => _input_record(last(value))
_input_record(value::AbstractDict) = Dict(_input_record(k) => _input_record(v) for (k,v) in value)

function _input_record(value::Union{DataModel.AbstractCablePart,DataModel.AbstractShape,
        DataModel.Pose2,DataModel.Shell,DataModel.Ring,
        DataModel.Helix,DataModel.LayRatio,DataModel.Pitch,
        DataModel.LayAngle,DataModel.FillFactor,DataModel.AssemblyMember,
        Materials.AbstractMaterial,Earth.AbstractEarthModel,Earth.AbstractEarthLayer,
        Earth.AbstractEarthMaterial})
    # Properties, not fields: a material's class is its type parameter, read as `kind`.
    names = propertynames(value)
    return merge((kind=nameof(typeof(value)),field_descriptions=Commons.input_fields(typeof(value))),
        NamedTuple{names}(map(name -> _input_record(getproperty(value,name)), names)))
end

_input_record(design::CableDesign) = (cable_id=design.cable_id,
    tree=_input_record(design.tree),
    terminal_order=design.terminal_order)

function _input_record(system::LineCableSystem)
    return (system_id=system.system_id, field_descriptions=Commons.input_fields(typeof(system)), line_length=system.line_length,
        designs=_input_record(system.designs), positions=_input_record(system.positions),
        connections=_input_record(system.connections), environment=_input_record(system.environment),
        terminal_order=system.terminal_order, connection_order=system.connection_order)
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
        problem.frequencies,(temperature=(name="temperature",unit="°C"),)))
end
function completed_inputs(problem::CableConstantsProblem)
    Record=NamedTuple{(:design,:temperature,:frequency,:field_descriptions),
        Tuple{NamedTuple,typeof(problem.temperature),typeof(problem.frequency),NamedTuple}}
    return Record((_input_record(problem.design),problem.temperature,problem.frequency,
        (temperature=(name="temperature",unit="°C"),frequency=(name="frequency",unit="Hz"))))
end

# Common-field omission compares selections, never punctuation.
"""
$(TYPEDSIGNATURES)

Read `formulation`, declared as `declaration`, as data with these fields:

- `formulations`: the declaration.
- `selections`: the formula identifier of each quantity, `Z` and `Y`.
- `formulation_fields`: the description of each selection and of its controls, for all
  quantities (`all`) and for each quantity.

A computed record contains the formulation, not these fields, and this method formats them
when a reader calls it.
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
            effective=get(declaration,:formulas,nothing)
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

Read `formulation` as data, for a point whose choices are empty.
"""
completed_formulation(formulation, ::Tuple{}) = completed_formulation(formulation)

"""
$(TYPEDSIGNATURES)

Read `formulation` as the line-parameter computation applied it on the earth plan `earth`.
Each earth formula records the equivalent earth that its calculation in the plan consumed. A
reduction that the plan resolved through `:default` appears as that formula's
`equivalent_earth`, as an explicit one does. The requested selections remain as declared.
"""
function completed_formulation(formulation::LineParametersFormulation, earth::EarthPlan)
    declaration = NamedTuple(formulation)
    slots = (:earth_impedance, :earth_admittance)
    # Each formula without an explicit `equivalent_earth` that the plan reduced, with the
    # reduction its calculation consumed.
    reduced = Pair[]
    for calculation in earth.calculations
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
    recorded = merge(declaration.formulas,
        map(record, formulation.formulas[slots], declaration.formulas[slots]))
    return completed_formulation(formulation,
        typeof(declaration)((declaration.backend, declaration.requested, recorded, declaration.options)))
end

"""
$(TYPEDSIGNATURES)

Read the formulation that `record` holds, with the choices its point made, as data. A
computed or loaded record holds the entries `formulation` and `choices`. A record without
them, such as that of a result built from matrices, reads as no formulation.
"""
completed_formulation(record::ComputationDetails) = haskey(record.data, :formulation) ?
    completed_formulation(record.data.formulation, record.data.choices) :
    (formulations = nothing, selections = (;), formulation_fields = (;))

function Commons.observation_gridpoint(source::Union{LineParameters,CableConstants})
    retained=details(source).data
    inputs=get(retained,:inputs,nothing)
    completed=completed_formulation(details(source))
    return Commons.detach((id=get(retained,:gridpoint,nothing), inputs,
        source_gridpoint=get(retained,:source_gridpoint,nothing),
        formulations=completed.formulations,formulation_fields=completed.formulation_fields,
        coordinates=source.terminals,
        uncertainty=get(retained,:uncertainty,nothing),
        uncertainty_descriptions=get(retained,:uncertainty_descriptions,nothing),
        transformation=get(retained,:modal,nothing),
        missing_reason=inputs===nothing ? :physical_inputs_not_supplied : nothing))
end

# The coefficient axes of a modal result are its modes, indexed and unnamed.
function Commons.observation_gridpoint(source::LineParameters{<:Any, <:Any, <:ModalDomain})
    phase=invoke(Commons.observation_gridpoint,Tuple{Union{LineParameters,CableConstants}},source)
    return merge(phase,(coordinates=string.(1:size(source.Z,1)),))
end

# Updating an association reuses numerical storage. It does not recompute or
# copy a physical problem. External result owners can supply their own method.
"""
$(TYPEDSIGNATURES)

Associate a completed result with its original gridpoint identity and additional
captured fields. Built-in results retain their numerical arrays and replace only
completion details. External result owners may extend this completion operation.
"""
retain_gridpoint(source::Union{LineParameters,CableConstants}, id; fields=(;)) =
    retain_gridpoint(source; fields=merge(fields,(gridpoint=id,)))

"""
$(TYPEDSIGNATURES)

Merge `fields` into the completion details of `source`. The gridpoint identity of `source`
is unchanged. Built-in results retain their numerical arrays and replace only completion details.
"""
function retain_gridpoint(source::LineParameters; fields)
    retained=ComputationDetails(merge(details(source).data,fields))
    return LineParameters(source.domain,source.Z,source.Y,source.f,source.terminals,
        source.line_length,retained)
end
function retain_gridpoint(source::CableConstants; fields)
    retained=ComputationDetails(merge(details(source).data,fields))
    return CableConstants(source.terminals,source.R,source.L,source.C,source.G,source.frequency,retained)
end
