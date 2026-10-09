"""
$(TYPEDEF)

Own one earth potential coefficient formulation, its indexed equation declarations and explicit
controls. Each interaction is selected by its explicit kind and source-layer and
target-layer equation signature. The layers in the signatures of its
`earth_potential_coefficient` methods declare the media it handles: methods up to layer 2 treat air
and one homogeneous earth. `parameters` stores model data. `options` stores formulation
choices and numerical controls.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions, E} <:
       EarthAdmittanceFormulation
    "Explicit physical model parameters."
    parameters::P
    "Formulation options. Projected onto required indexed equations during initialization."
    options::O
    "Independent equivalent homogeneous-earth reduction, or nothing."
    equivalent_earth::E
end

"""
Evaluate one source-owned earth potential coefficient equation.
"""
function earth_potential_coefficient end

"""
$(TYPEDSIGNATURES)

Resolve a formulation's model parameters and numerical controls. Its concrete
selection type defines the indexed equation. Material properties are evaluated
before equation execution.
Formula-specific arguments are normalized by their selected owner.
A missing physical case is unsupported.
"""
function Formula{ID}(; parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions(), equivalent_earth = nothing) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    options = formulation_options(Formula{ID}, options)
    isempty(parameters) ||
        throw(ArgumentError("earth formula :$ID has no configurable physical parameters"))
    ID in formulas(Formula) || throw(ArgumentError("unknown earthadmittance formula :$ID"))
    reduction = if equivalent_earth === nothing
        nothing
    elseif equivalent_earth isa EquivalentHomogeneous.AbstractSequence
        equivalent_earth
    elseif equivalent_earth isa FormulaDefinition
        EquivalentHomogeneous.AbstractSequence(equivalent_earth)
    else
        throw(ArgumentError("equivalent_earth must be a formula definition or explicit reduction sequence"))
    end
    return Formula{ID, typeof(parameters), typeof(options), typeof(reduction)}(
        parameters, options, reduction)
end

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(::Val{ID}; kwargs...) where {ID} = Formula{ID}(; kwargs...)
Formula(selected::EarthAdmittanceFormulation) = selected
Formula(::Nothing) = Formula(:default)

# A recipe selects one formula per earth route. Within it, a `nothing` route supplies no
# expression. Only omitting the whole slot selects the default.
function Formula(recipe::NamedTuple)
    routes = (; pairs(Formula)...)
    all(in(keys(routes)), keys(recipe)) || throw(ArgumentError(
        "$Formula selections admit only $(join(keys(routes), ", "))"))
    names = filter(in(keys(recipe)), keys(routes))
    return NamedTuple{names}(map(name -> recipe[name] === nothing ? nothing :
        Formula(recipe[name]), names))
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    return Formula{ID}(; parameters = selection.parameters,
        options = selection.options, equivalent_earth = selection.equivalent_earth)
end

formula_id(::Formula{ID}) where {ID} = ID
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
# Identity-only dispatch also describes retained selections without constructors.
description(value::Formula; compact::Bool = false) = description(typeof(value); compact)
formulation_options(value::Formula) = value.options
# Indexed numerical controls remain deferred to their consuming equations.
formulation_options(::Type{<:Formula}, options::FormulationOptions) = options

"""
$(TYPEDSIGNATURES)

Expose the selected identity, model and numerical controls, and explicit
equivalent-earth reduction as a native record.
"""
function Base.NamedTuple(value::Formula)
    return (identifier = formula_id(value),
        parameters = value.parameters, options = value.options.data,
        equivalent_earth = value.equivalent_earth === nothing ? nothing :
                           NamedTuple(value.equivalent_earth))
end

"""Iterate the independently selectable child slots admitted by this formula family."""
function Base.pairs(::Type{<:Formula}; quantity = nothing)
    pairs((air = Formula, earth = Formula, mixed = Formula))
end
