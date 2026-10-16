"""
Interface for a selected frequency-dependent earth-material relation.
Concrete selections expose model `parameters` and numerical `options` records.
"""
abstract type FrequencyDependentFormulation <: AbstractFormulation end

"""
$(TYPEDEF)

Select one frequency-dependent earth-material relation by its stable formula
identifier.

Each formula implements `earth_material(selected, functor, workspace) -> EarthMaterial` on
its concrete selection type. The input of `functor` holds the earth material, the frequency,
the complex frequency `s` and the options. `:default` is a routing alias for the explicit `:constant` pass-through.

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
    options = FormulationOptions(options)
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

"""
$(TYPEDSIGNATURES)

Build the Functor of a frequency-dependent earth relation for one earth material at one
frequency, after checking the frequency. The state is empty.
"""
function Functor(formula::FrequencyDependentFormulation, input::NamedTuple;
        workspace = nothing)
    (; frequency) = input
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency, "earth-property evaluation frequency must be positive and finite"))
    return Functor(formula, input, (;))
end

@inline function (formula::FrequencyDependentFormulation)(
        material::EarthMaterial{T}, frequency::T; workspace = nothing,
        s = complex(zero(T), 2 * (one(T) * π) * frequency)
) where {T <: Real}
    validate(s, formula)
    functor = Functor(formula, (; material, frequency, s, options = formula.options); workspace)
    return validate(Expression(formula, earth_material)(functor, workspace), formula)
end

# The material that a frequency-dependent earth relation returns.
validate(material::EarthMaterial, formula::FrequencyDependentFormulation) = material
validate(s::Complex, ::FrequencyDependentFormulation) = s

function validate(s::Complex, formula::Union{Formula{:scott1967}, Formula{:visacro1987},
        Formula{:visacro2012}, Formula{:datsios2019}, Formula{:cigre2019}})
    iszero(real(s)) || throw(ArgumentError(
        "$(formula_id(formula)) is a fitted law that is not causal and has no value at a complex frequency; at σ ≠ 0 use a causal law: :constant, :portela1999, :alipio2014, :longmire1975 or :messier1985"))
    return s
end

function validate(material, formula::FrequencyDependentFormulation)
    throw(ArgumentError("a frequency-dependent earth relation must return EarthMaterial; " *
                        "got $(typeof(material))"))
end

function (formula::FrequencyDependentFormulation)(
        material::EarthMaterial{T}, frequency::Real; workspace = nothing,
        s = complex(zero(T), 2 * (one(T) * π) * frequency)) where {T <:
                                                                                 Real}
    U = promote_type(T, typeof(float(frequency)))
    return formula(
        convert(EarthMaterial{U}, material),
        convert(U, float(frequency)); workspace, s = convert(Complex{U}, s)
    )
end

"""
Pass static earth properties through when no constitutive relation is selected.
"""
constitutive(::Nothing, material::EarthMaterial, ::Real; workspace = nothing,
    s = nothing) = material

"""
Evaluate one registered frequency-dependent earth constitutive relation.
"""
function constitutive(
        formula::FrequencyDependentFormulation, material::EarthMaterial{T}, frequency::Real;
        workspace = nothing,
        s = complex(zero(T), 2 * (one(T) * π) * frequency)) where {T <: Real}
    formula(material, frequency; workspace, s)
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
