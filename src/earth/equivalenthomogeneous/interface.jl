"""
Abstract equivalent homogeneous-earth rule.
"""
abstract type AbstractRule <: AbstractFormulation end

"""
Abstract ordering of material frequency dependence and EquivalentHomogeneous reduction.
"""
abstract type AbstractSequence <: AbstractFormulation end

"""
$(TYPEDEF)

Select one equivalent homogeneous-earth rule by its stable formula
identifier.

The rule receives evaluated layer properties, the
static earth geometry, the interaction pair, and frequency. It returns one
artificial homogeneous [`EarthMaterial`](@ref). Formula parameters participate
in the concrete Julia type.

The `:default` formula uses the bottommost soil layer as the equivalent
homogeneous material.

$(TYPEDFIELDS)
"""
struct Formula{ID, A <: NamedTuple, O <: FormulationOptions} <: AbstractRule
    "Explicit model parameters."
    parameters::A
    "Explicit numerical sections owned by the reduction."
    options::O
end

"""
$(TYPEDEF)

Apply material frequency dependence to every physical layer before the EquivalentHomogeneous
rule constructs an equivalent material.

$(TYPEDFIELDS)
"""
struct AfterFD{R <: AbstractRule} <: AbstractSequence
    "Equivalent homogeneous-earth rule."
    rule::R
end

"""
$(TYPEDEF)

Construct an equivalent material from static layers before applying the
selected material frequency dependence to that artificial material.

$(TYPEDFIELDS)
"""
struct BeforeFD{R <: AbstractRule} <: AbstractSequence
    "Equivalent homogeneous-earth rule."
    rule::R
end

"""
Return the rule stored by an EquivalentHomogeneous composition.
"""
rule(sequence::AbstractSequence) = sequence.rule

"""
Return the stable formula identifier of an EquivalentHomogeneous formula.
"""
formula_id(::Formula{ID}) where {ID} = ID

"""
Construct one formula-owned equivalent homogeneous-earth material.
"""
function equivalent_material end

"""
$(TYPEDSIGNATURES)

Construct an equivalent-earth rule with model parameters and numerical controls.
Custom rules subtype `AbstractRule` and extend `equivalent_material` on their
own concrete type. The selected sequence owns its position relative to the
frequency-dependent material law.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::AbstractRule) = selected

function Formula(::Val{:bottommost}; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions())
    options = options isa NamedTuple ? FormulationOptions(options) : options
    isempty(parameters) || throw(ArgumentError("bottommost earth has no configurable model parameters"))
    return Formula{:bottommost, typeof(parameters), typeof(options)}(parameters, options)
end

Formula(::Val{ID}; kwargs...) where {ID} = throw(ArgumentError("unknown equivalent-earth rule :$ID"))

AfterFD(identifier::Symbol; kwargs...) = AfterFD(Formula(identifier; kwargs...))
BeforeFD(identifier::Symbol; kwargs...) = BeforeFD(Formula(identifier; kwargs...))

function description(sequence::AfterFD;compact::Bool=false)
    "$(description(sequence.rule;compact)) after layerwise FrequencyDependent"
end
function description(sequence::BeforeFD;compact::Bool=false)
    "$(description(sequence.rule;compact)) before layerwise FrequencyDependent"
end

@inline function (formula::AbstractRule)(
        rho::AbstractVector,
        eps_r::AbstractVector,
        mu_r::AbstractVector,
        model::EarthModel,
        pair,
        frequency::Real; binding = validate(formula, pair), workspace = nothing
)
    length(rho) == length(eps_r) == length(mu_r) == length(model.layers) ||
        throw(DimensionMismatch("EquivalentHomogeneous properties must align with the complete physical model"))
    isfinite(frequency) && frequency > zero(frequency) || throw(DomainError(
        frequency,
        "EquivalentHomogeneous evaluation frequency must be positive and finite"
    ))
    validate(pair, getproperty.(model.layers, :thickness))
    material = binding.equation(rho, eps_r, mu_r, model, pair, frequency,
        formula.parameters, binding.options, workspace)
    material isa EarthMaterial ||
        throw(ArgumentError("an EquivalentHomogeneous contribution must return EarthMaterial"))
    return material
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing ||
        throw(ArgumentError("a reduction cannot contain another reduction"))
    return Formula(Val(ID); parameters = selection.parameters,
        options = selection.options)
end

function AbstractSequence(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    selection.equivalent_earth === nothing ||
        throw(ArgumentError("a reduction cannot contain another reduction"))
    rule = Formula(Val(ID); parameters = selection.parameters,
        options = selection.options)
    return Order === :before ? BeforeFD(rule) : AfterFD(rule)
end

function equivalent_material(selected::AbstractRule, ::Val{Kind}, ::Val{S}, ::Val{T},
        rho, eps_r, mu_r, model, pair, frequency, parameters,
        options, workspace
) where {Kind, S, T}
    throw(ArgumentError("equivalent_material :$(formula_id(selected)) ($Kind): formula not implemented for source in layer $S and target in layer $T"))
end

validate(formula::AbstractRule, pair) = only(validate(formula, (pair,)))

function validate(formula::AbstractRule, pairs::Union{Tuple, AbstractVector})
    equations = map(pairs) do pair
        kind = pair.row == pair.column ? :self : :mutual
        FormulaMethod(formula, equivalent_material, Val(kind), Val.(pair.layers)...)
    end
    identities = unique(equations)
    defaults = map(identities) do binding
        formulation_options(binding)
    end
    admitted = union((keys(value.data) for value in defaults)...)
    isempty(setdiff(keys(formula.options.data), admitted)) ||
        throw(ArgumentError("unused equivalent-earth numerical sections"))
    resolved = map(identities, defaults) do binding, declared
        names = Tuple(intersect(keys(formula.options.data), keys(declared.data)))
        (equation = binding,
            options = formulation_options(binding, declared, FormulationOptions(formula.options.data[names])))
    end
    return map(equation -> resolved[findfirst(==(equation), identities)], equations)
end

"""Expose the reduction rule, model parameters and numerical options as a native record."""
function Base.NamedTuple(value::Formula)
    return (identifier=formula_id(value), parameters=value.parameters,
        options=value.options.data)
end

"""Expose the order of material evaluation and the selected equivalent-earth rule."""
function Base.NamedTuple(value::AbstractSequence)
    return (order=nameof(typeof(value)), rule=NamedTuple(value.rule))
end

# Identity-only dispatch also describes retained selections without constructors.
import ...Grammar: formulation_options
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)

"""Iterate the independently selectable child slots admitted by this formula family."""
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((;))
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
formulation_options(value::Formula) = value.options

"""Describe an explicit equivalent-earth rule and its requested material-law order."""
description(slot::Val{:equivalent_earth},value::FormulaDefinition;compact::Bool=true) =
    description(slot,NamedTuple(value);compact)
function description(::Val{:equivalent_earth},value::NamedTuple;compact::Bool=true)
    record=get(value,:rule,value)
    selected=Formula{record.identifier}
    text=applicable(description,selected) ? description(selected;compact) : string(record.identifier)
    order=get(value,:order,:default)
    order=order===:BeforeFD ? :before : order===:AfterFD ? :after : order
    order===:default || (text *= " "*string(order)*" FrequencyDependent")
    controls=(; (key => record[key] for key in (:parameters, :options)
        if haskey(record,key) && !isempty(record[key]))...)
    isempty(controls) || (text *= " "*description(FormulaDefinition,controls;compact))
    return text
end
description(::Val{:equivalent_earth},value::AbstractSequence;compact::Bool=true) = description(value;compact)
