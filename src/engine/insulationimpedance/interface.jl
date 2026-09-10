"""
$(TYPEDEF)

Select one insulation-impedance formula by its stable identifier.

$(TYPEDFIELDS)
"""
struct Formula{ID, R, A <: NamedTuple, H <: NamedTuple, O <: NamedTuple} <:
       InsulationImpedanceFormulation
    "Declared series-impedance equation binding."
    binding::R
    "Explicit model parameters."
    parameters::A
    "Explicit callable overrides retained for provenance."
    hooks::H
    "Normalized numerical sections for the selected contribution."
    options::O
end

"Return the stable identifier of an insulation-impedance formula."
formula_id(::Formula{ID}) where {ID} = ID

"Evaluate one formula-owned insulation-impedance route."
function insulation_impedance end

"""
$(TYPEDSIGNATURES)

Construct a registered formula with separate model parameters and callable
hooks. `hooks=(contribution=f,)` replaces the complete scalar equation using
the signature `f(r_in, r_ex, mu_r, jω, parameters, options, workspace) → complex impedance [Ω/m]`. A complete replacement declares its numerical defaults with
`computation_options(binding, replacement)`. Unknown fields fail immediately.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)
Formula(selected::Formula) = selected

function Formula(::Val{ID}; parameters::NamedTuple = (;), hooks::NamedTuple = (;),
        options::NamedTuple = (;)) where {ID}
    ID in FORMULAS || throw(ArgumentError("unknown formula :$ID"))
    isempty(parameters) ||
        throw(ArgumentError("formula :$ID has no configurable model parameters"))
    isempty(setdiff(keys(hooks), (:contribution,))) ||
        throw(ArgumentError("unknown hooks for :$ID"))
    binding = FormulaMethod(Val(ID), insulation_impedance)
    selected = get(hooks, :contribution, binding)
    selected === nothing && throw(ArgumentError("a contribution hook must be callable"))
    defaults = haskey(hooks, :contribution) ? computation_options(binding, selected) :
               computation_options(binding)
    normalized = computation_options(binding, defaults, options)
    return Formula{
        ID, typeof(binding), typeof(parameters), typeof(hooks), typeof(normalized)}(
        binding, parameters, hooks, normalized)
end

@inline function (formula::Formula)(
        r_in::T,
        r_ex::T,
        mu_r::T,
        s::Complex{T}; workspace = nothing
) where {T <: Real}
    isfinite(r_in) && isfinite(r_ex) && zero(T) <= r_in <= r_ex ||
        throw(DomainError((r_in, r_ex), "insulation radii must satisfy 0 ≤ r_in ≤ r_ex [m]"))
    isfinite(mu_r) && mu_r > zero(T) || throw(DomainError(mu_r,
        "relative insulation permeability must be positive and finite"))
    isfinite(s) || throw(DomainError(s, "jω must be finite"))
    value = get(formula.hooks, :contribution, formula.binding)(
        r_in, r_ex, mu_r, s, formula.parameters, formula.options, workspace)
    value isa Number && isfinite(value) || throw(DomainError(value,
        "insulation_impedance must return a finite scalar"))
    return value
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula(Val(ID); parameters = selection.parameters, hooks = selection.hooks,
        options = selection.options)
end

"""
$(TYPEDSIGNATURES)

Expose the selected equation, physical parameters, callable overrides and numerical
options as a native record. Callables are retained unchanged.
"""
function Base.NamedTuple(value::Formula)
    return (identifier=formula_id(value), binding=value.binding,
        parameters=value.parameters, hooks=value.hooks, options=value.options)
end
