"""
Interface for a selected frequency-dependent earth-material relation.
Concrete selections expose model `parameters` and numerical `options` records.
"""
abstract type FrequencyDependentFormulation <: AbstractFormulation end

"""
$(TYPEDEF)

Select one frequency-dependent earth-material relation by its stable formula
identifier.

Each formula implements
`earth_material(selected, material, frequency, parameters, options, workspace) -> EarthMaterial`
on its concrete selection type. `:default` is a routing
alias for the explicit `:constant` pass-through.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions} <:
       FrequencyDependentFormulation
    "Resolved physical model parameters."
    parameters::P
    "Normalized numerical sections for this equation."
    options::O
end

"""
Evaluate one formula-owned frequency-dependent earth material relation.
"""
function earth_material end

"""
$(TYPEDSIGNATURES)

Construct a selected formulation with model parameters and numerical controls.
Custom formulations extend `earth_material` on their own concrete selection type.
Unknown controls fail before numerical evaluation. A formula with model parameters
defines its own identity constructor, which holds their defaults.
"""
function Formula{ID}(; parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    return Formula{ID}((;), parameters, options)
end

"""
$(TYPEDSIGNATURES)

Complete the supplied model `parameters` with a formula's `defaults`, check every
coefficient and the formula's coefficient domains, and normalize its numerical controls.
"""
function Formula{ID}(defaults::NamedTuple, parameters::NamedTuple,
        options::Union{NamedTuple, FormulationOptions}) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    unknown = setdiff(keys(parameters), keys(defaults))
    isempty(unknown) || throw(ArgumentError(
        "unknown parameters for earth-property formula :$ID: $(collect(unknown))"))
    parameters = merge(defaults, parameters)
    for (name, value) in pairs(parameters)
        value isa Real && !(value isa Bool) && isfinite(value) ||
            throw(ArgumentError("earth-property parameter :$name must be a finite real coefficient"))
    end
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    validate(selected)
    expression = Expression(selected, earth_material)
    normalized = formulation_options(expression, options)
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

"""Check model-specific coefficient domains before evaluating a material."""
validate(selected::Formula) = selected

@inline function (formula::FrequencyDependentFormulation)(
        material::EarthMaterial{T}, frequency::T; workspace = nothing
) where {T <: Real}
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency,
        "earth-property evaluation frequency must be positive and finite"
    ))
    evaluated = earth_material(
        formula, material, frequency, formula.parameters, formula.options, workspace)
    evaluated isa EarthMaterial ||
        throw(ArgumentError("a frequency-dependent earth relation must return EarthMaterial"))
    return evaluated
end

function (formula::FrequencyDependentFormulation)(
        material::EarthMaterial{T}, frequency::Real; workspace = nothing) where {T <:
                                                                                 Real}
    U = promote_type(T, typeof(float(frequency)))
    return formula(
        convert(EarthMaterial{U}, material),
        convert(U, float(frequency)); workspace
    )
end

"""
Pass static earth properties through when no constitutive relation is selected.
"""
constitutive(::Nothing, material::EarthMaterial, ::Real; workspace = nothing) = material

"""
Evaluate one registered frequency-dependent earth constitutive relation.
"""
function constitutive(
        formula::FrequencyDependentFormulation, material::EarthMaterial, frequency::Real;
        workspace = nothing)
    formula(material, frequency; workspace)
end

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(::Val{ID}; kwargs...) where {ID} = Formula{ID}(; kwargs...)
Formula(selected::FrequencyDependentFormulation) = selected

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula{ID}(; parameters = selection.parameters,
        options = selection.options)
end

"""
Return the stable formula identifier of an earth-property formula.
"""
formula_id(::Formula{ID}) where {ID} = ID
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
# Identity-only dispatch also describes retained selections without constructors.
description(value::Formula; compact::Bool = false) = description(typeof(value); compact)
formulation_options(value::Formula) = value.options

"""
$(TYPEDSIGNATURES)

Expose the selected identity, model parameters, and numerical options as a native record.
"""
function Base.NamedTuple(value::Formula)
    return (identifier = formula_id(value),
        parameters = value.parameters, options = value.options.data)
end

"""Iterate the independently selectable child slots admitted by this formula family."""
Base.pairs(::Type{<:Formula}; quantity = nothing) = pairs((;))
