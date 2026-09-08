"Abstract equivalent homogeneous-earth rule."
abstract type AbstractRule <: AbstractFormulation end

"Abstract ordering of material frequency dependence and EquivalentHomogeneous reduction."
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
struct Formula{ID, A <: NamedTuple, H <: NamedTuple, O <: NamedTuple} <: AbstractRule
    "Explicit model parameters."
    parameters::A
    "Explicit callable overrides retained for provenance."
    hooks::H
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

"Return the rule stored by an EquivalentHomogeneous composition."
rule(sequence::AbstractSequence) = sequence.rule

"Return the stable formula identifier of an EquivalentHomogeneous formula."
formula_id(::Formula{ID}) where {ID} = ID

"Construct one formula-owned equivalent homogeneous-earth material."
function equivalent_material end

"""
$(TYPEDSIGNATURES)

Construct a registered formula with separate model parameters and callable
hooks. `hooks=(contribution=f,)` replaces the complete scalar equation using
the signature `f(rho, eps_r, mu_r, model, pair, frequency, parameters, options, workspace) → EarthMaterial`. A complete replacement declares its numerical defaults with
`computation_options(binding, replacement)`. Unknown fields fail immediately.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::Formula) = selected

function Formula(::Val{ID}; parameters::NamedTuple = (;),
        hooks::NamedTuple = (;), options::NamedTuple = (;)) where {ID}
    ID in FORMULAS || throw(ArgumentError("unknown formula :$ID"))
    isempty(parameters) ||
        throw(ArgumentError("formula :$ID has no configurable model parameters"))
    isempty(setdiff(keys(hooks), (:contribution,))) ||
        throw(ArgumentError("unknown hooks for :$ID"))
    haskey(hooks, :contribution) && hooks.contribution === nothing &&
        throw(ArgumentError("a contribution hook must be callable"))
    return Formula{ID, typeof(parameters), typeof(hooks), typeof(options)}(
        parameters, hooks, options)
end

AfterFD(identifier::Symbol; kwargs...) = AfterFD(Formula(identifier; kwargs...))
BeforeFD(identifier::Symbol; kwargs...) = BeforeFD(Formula(identifier; kwargs...))

function description(sequence::AfterFD)
    "$(description(sequence.rule)) after layerwise FrequencyDependent"
end
function description(sequence::BeforeFD)
    "$(description(sequence.rule)) before layerwise FrequencyDependent"
end

@inline function (formula::Formula)(
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
    selected = get(formula.hooks, :contribution, binding.equation)
    material = selected(rho, eps_r, mu_r, model, pair, frequency,
        formula.parameters, binding.options, workspace)
    material isa EarthMaterial ||
        throw(ArgumentError("an EquivalentHomogeneous contribution must return EarthMaterial"))
    return material
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing ||
        throw(ArgumentError("a reduction cannot contain another reduction"))
    return Formula(Val(ID); parameters = selection.parameters, hooks = selection.hooks,
        options = selection.options)
end

function AbstractSequence(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    selection.equivalent_earth === nothing ||
        throw(ArgumentError("a reduction cannot contain another reduction"))
    rule = Formula(Val(ID); parameters = selection.parameters, hooks = selection.hooks,
        options = selection.options)
    return Order === :before ? BeforeFD(rule) : AfterFD(rule)
end

function equivalent_material(::Val{ID}, ::Val{Kind}, ::Val{S}, ::Val{T},
        rho, eps_r, mu_r, model, pair, frequency, parameters,
        options, workspace
) where {ID, Kind, S, T}
    throw(ArgumentError("equivalent_material :$ID ($Kind): formula not implemented for source in layer $S and target in layer $T"))
end

const EQUATION_FALLBACK = which(equivalent_material,
    Tuple{Val, Val, Val, Val, Any, Any, Any, Any, Any, Any, Any, Any, Any})

validate(formula::Formula, pair) = only(validate(formula, (pair,)))

function validate(formula::Formula{ID}, pairs::Union{Tuple, AbstractVector}) where {ID}
    equations = map(pairs) do pair
        kind = pair.row == pair.column ? :self : :mutual
        binding = FormulaMethod(Val(ID), equivalent_material, Val(kind), Val.(pair.layers)...)
        signature = Tuple{Val{ID}, typeof.(binding.arguments)...,
            Any, Any, Any, Any, Any, Any, Any, Any, Any}
        which(equivalent_material, signature) === EQUATION_FALLBACK &&
            binding(nothing, nothing, nothing, nothing, nothing,
                nothing, nothing, nothing, nothing)
        binding
    end
    identities = unique(equations)
    defaults = map(identities) do binding
        haskey(formula.hooks, :contribution) ?
        computation_options(binding, formula.hooks.contribution) :
        computation_options(binding)
    end
    admitted = union((keys(value) for value in defaults)...)
    isempty(setdiff(keys(formula.options), admitted)) ||
        throw(ArgumentError("unused equivalent-earth numerical sections"))
    resolved = map(identities, defaults) do binding, declared
        names = Tuple(intersect(keys(formula.options), keys(declared)))
        (equation = binding,
            options = computation_options(binding, declared, formula.options[names]))
    end
    return map(equation -> resolved[findfirst(==(equation), identities)], equations)
end
