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
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions} <: FrequencyDependentFormulation
    "Resolved physical/model parameters."
    parameters::P
    "Normalized numerical sections for this equation."
    options::O
end

"""
Return the stable formula identifier of an earth-property formula.
"""
formula_id(::Formula{ID}) where {ID} = ID

"""
Return the resolved physical parameters of a frequency-dependent earth formula.
"""
assumptions(formula::FrequencyDependentFormulation) = formula.parameters

"""Return the default physical parameters for a registered earth formula."""
function assumptions end

"""Return the default physical parameters of a formula identifier."""
assumptions(::Val{ID}) where {ID} = (;)

"Return vacuum permittivity represented in the scalar type of `value` \\[F/m\\]."
@inline vacuum_permittivity(value) =
    one(value) * 88541878128 * (one(value) * 10)^(-22)

"""
Evaluate one formula-owned frequency-dependent earth material relation.
"""
function earth_material end

"""
$(TYPEDSIGNATURES)

Construct a selected formulation with model parameters and numerical controls.
Custom formulations extend `earth_material` on their own concrete selection type.
Unknown controls fail before numerical evaluation.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::FrequencyDependentFormulation) = selected

function Formula(::Val{ID}; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    defaults = assumptions(Val(ID))
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
    binding = FormulaMethod(selected, earth_material)
    normalized = formulation_options(binding, options)
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

function (formula::FrequencyDependentFormulation)(material::EarthMaterial{T}, frequency::Real; workspace = nothing) where {T <:
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
constitutive(::Nothing, material::EarthMaterial, ::Real) = material

"""
Evaluate one registered frequency-dependent earth constitutive relation.
"""
function constitutive(formula::FrequencyDependentFormulation, material::EarthMaterial, frequency::Real)
    formula(material, frequency)
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula(Val(ID); parameters = selection.parameters,
        options = selection.options)
end

"""
$(TYPEDSIGNATURES)

Expose the selected identity, model parameters, and numerical options as a native record.
"""
function Base.NamedTuple(value::Formula)
    return (identifier=formula_id(value), parameters=value.parameters, options=value.options.data)
end

# Identity-only dispatch also describes retained selections without constructors.
import ...Grammar: formulation_options
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)

"""Iterate the independently selectable child slots admitted by this formula family."""
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((;))
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
formulation_options(value::Formula) = value.options
