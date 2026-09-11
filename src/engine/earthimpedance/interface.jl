"""
$(TYPEDEF)

Own one earth impedance formulation, its indexed equation declarations and explicit
customizations. `assumptions` contains scientific restrictions; each interaction
is selected by its explicit kind/source-layer/target-layer equation signature.
`parameters` and `hooks` retain the user's explicit modifications for provenance.

$(TYPEDFIELDS)
"""
struct Formula{ID, A <: NamedTuple, P <: NamedTuple, H <: NamedTuple, O <: NamedTuple, E} <:
       EarthImpedanceFormulation
    "Model class, exact medium inventory and longitudinal restriction."
    assumptions::A
    "Explicit physical/model parameters."
    parameters::P
    "Explicit callable overrides; empty for an unmodified formulation."
    hooks::H
    "Explicit numerical sections; projected onto required indexed equations at preflight."
    options::O
    "Independent equivalent homogeneous-earth reduction, or nothing."
    equivalent_earth::E
end

"""
$(TYPEDEF)

Bind one indexed interaction to evaluated physical state. Numerical options and
callable hooks remain separate from physical quantities. Mutable numerical
resources are passed when evaluating the functor.

$(TYPEDFIELDS)
"""
struct Functor{ID, B, H, S, O}
    "Selected equation and its exact source/target geometry."
    binding::B
    "Resolved callable hooks retained unchanged during evaluation."
    hooks::H
    "Evaluated physical quantities, including one authoritative Γ [1/m]."
    state::S
    "Normalized computation options."
    options::O
end

formula_id(::Formula{ID}) where {ID} = ID
assumptions(formula::Formula) = formula.assumptions
media(formula::Formula) = formula.assumptions.media

"Declare model inventory and physical restrictions, without callable behavior."
function assumptions end
"Evaluate a medium propagation law at jω [1/s], μ [H/m], σ [S/m], ε [F/m]."
function propagation end
"Evaluate the prescribed longitudinal Γ [1/m] from jω, evaluated materials and (s,t)."
function Γ end
"Evaluate one source-owned earth impedance equation."
function earth_impedance end

Γ(functor::Functor) = functor.state.Γ

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

"""
$(TYPEDSIGNATURES)

Resolve a registered formula's parameters and callable hooks once. Permitted
hooks are `Γ(jω, materials, (s,t))`, medium laws `air/earth(jω, μ, σ, ε)`,
`permeability(μ)`, and `contribution(functor, pair, workspace)`. A contribution
override remains subject to the source's declared domain. Built-in equations
never use another formula's overrides to supply a missing case.
"""
function Formula(::Val{ID}; parameters::NamedTuple = (;), hooks::NamedTuple = (;),
        options::NamedTuple = (;), equivalent_earth = nothing) where {ID}
    ID in FORMULAS || throw(ArgumentError("unknown earth-impedance formula :$ID"))
    parameters = earth_parameters(Val(ID), parameters)
    # Indexed declarations own admitted hook names. Resolve them with the actual
    # required cases, alongside case-local numerical sections, during preflight.
    any(isnothing, values(hooks)) &&
        throw(ArgumentError("an explicit hook must be callable, not nothing"))
    declared_constraints = assumptions(Val(ID))
    constraints = merge(declared_constraints, (media = Val(declared_constraints.media),))
    reduction = if equivalent_earth === nothing
        nothing
    elseif equivalent_earth isa EquivalentHomogeneous.AbstractSequence
        equivalent_earth
    elseif equivalent_earth isa FormulaDefinition
        EquivalentHomogeneous.AbstractSequence(equivalent_earth)
    else
        throw(ArgumentError("equivalent_earth must be a formula definition or explicit reduction sequence"))
    end
    reduction !== nothing && constraints.media !== Val(:homogeneous) &&
        throw(ArgumentError("a full multilayer formula cannot consume an equivalent homogeneous reduction"))
    return Formula{ID, typeof(constraints), typeof(parameters), typeof(hooks),
        typeof(options), typeof(reduction)}(
        constraints, parameters, hooks, options, reduction)
end

"""
$(TYPEDSIGNATURES)

Prepare one validated indexed interaction at fixed angular frequency. Material
vectors use physical indices for a stratified equation and exactly `(air,soil)`
for a homogeneous equation. Explicit reduction and its physical/effective
mapping are owned by the computation workspace.

# Arguments

- `resistivity`: Aligned resistivities [Ω·m].
- `permittivity`: Absolute permittivities [F/m].
- `permeability`: Absolute permeabilities [H/m].
- `jω`: Imaginary angular frequency [1/s].
- `pair`: Indexed conductor interaction; lengths [m].

# Keywords

- `Γ`: Optional explicit longitudinal constant [1/m]; conflicts with a Γ hook.
- `thickness`: Aligned layer thicknesses [m] for a stratified model.
"""
function (formula::Formula{ID})(
        resistivity::AbstractVector{T}, permittivity::AbstractVector{T},
        permeability::AbstractVector{T}, jω::Complex{T}, pair::EarthPair;
        Γ = nothing, thickness = nothing, physical_pair = pair
) where {ID, T <: Real}
    selected = validate(formula, pair)
    return formula(resistivity, permittivity, permeability, jω, pair, selected;
        Γ, thickness, physical_pair)
end

# Evaluate a declaration already bound by validate() before the frequency loop.
function (formula::Formula{ID})(
        resistivity::AbstractVector{T}, permittivity::AbstractVector{T},
        permeability::AbstractVector{T}, jω::Complex{T}, pair::EarthPair, selected;
        Γ = nothing, thickness = nothing, physical_pair = pair
) where {ID, T <: Real}
    isfinite(jω) && !iszero(jω) || throw(DomainError(jω, "jω must be finite and nonzero"))
    options = selected.options
    validate(formula, resistivity, permittivity, permeability, thickness)
    thickness === nothing || validate(pair, thickness)
    Γ !== nothing && haskey(formula.hooks, :Γ) &&
        throw(ArgumentError(
            "an explicit problem Γ conflicts with the explicit :$ID Γ hook"))
    layers = media(formula) === Val(:homogeneous) ? (1, 2) : eachindex(permeability)
    μ = map(layers) do layer
        layer == 1 ? permeability[layer] : selected.hooks.permeability(permeability[layer])
    end
    σ = map(layer -> conductivity(resistivity[layer]), layers)
    γ = map(layers) do layer
        law = layer == 1 ? selected.hooks.air : selected.hooks.earth
        law(jω, μ[layer], σ[layer], permittivity[layer])
    end
    all(value -> value isa Number && isfinite(value), γ) ||
        throw(DomainError(γ, "medium propagation laws must return finite scalars [1/m]"))
    all(value -> value isa Real && isfinite(value) && value > 0, μ) ||
        throw(DomainError(μ, "permeability hooks must return positive finite scalars [H/m]"))
    materials = (rho = resistivity, epsilon = permittivity, mu = μ, sigma = σ,
        gamma = γ, gamma_medium_squared = γ .^ 2)
    longitudinal = Γ === nothing ? selected.hooks.Γ(jω, materials, pair.layers) : Γ
    longitudinal isa Number && isfinite(longitudinal) || throw(ArgumentError(
        "Γ must be one finite scalar [1/m], not a value/square pair"))
    formula.assumptions.longitudinal === :zero && !iszero(longitudinal) &&
        throw(ArgumentError("earth-impedance :$ID fixes Γ to zero"))
    state = merge(materials, (; jω, Γ = oftype(jω, longitudinal), thickness))
    binding = (pair = pair, physical_pair = physical_pair,
        kind = selected.kind, equation = selected.equation)
    return Functor{
        ID, typeof(binding), typeof(selected.hooks), typeof(state), typeof(options)}(
        binding, selected.hooks, state, options)
end

function (functor::Functor)(workspace = nothing)
    pair = functor.binding.pair
    value = if functor.hooks.contribution === nothing
        functor.binding.equation(functor, pair, workspace)
    else
        functor.hooks.contribution(functor, pair, workspace)
    end
    value isa Number && isfinite(value) || throw(DomainError(value,
        "earth-impedance contribution must be a finite scalar"))
    return oftype(functor.state.jω, value)
end

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    return Formula(Val(ID); parameters = selection.parameters, hooks = selection.hooks,
        options = selection.options, equivalent_earth = selection.equivalent_earth)
end

Formula(selected::Formula) = selected

function FormulaMethod(formula::Formula{ID}, pair::EarthPair) where {ID}
    return FormulaMethod(Val(ID), earth_impedance,
        Val(pair.row == pair.column ? :self : :mutual), Val.(pair.layers)...)
end

function earth_impedance(::Val{ID}, ::Val{Kind}, ::Val{S}, ::Val{T},
        functor, pair, workspace) where {ID, Kind, S, T}
    throw(ArgumentError(
        "earth_impedance :$ID ($Kind): formula not implemented for source in layer $S and target in layer $T"))
end

const EQUATION_FALLBACK = which(earth_impedance, Tuple{Val, Val, Val, Val, Any, Any, Any})

function validate(binding::FormulaMethod{ID, typeof(earth_impedance), A}) where {ID, A}
    signature = Tuple{Val{ID}, A.parameters..., Any, Any, Any}
    if which(earth_impedance, signature) === EQUATION_FALLBACK
        # The diagnostic takes no physical inputs and evaluates no numerical kernel.
        binding(nothing, nothing, nothing)
    end
    return binding
end

"""
$(TYPEDSIGNATURES)

Expose the selected equations, scientific restrictions, parameters, hooks,
numerical options and explicit equivalent-earth reduction as a native record.
Callables are retained unchanged.
"""
function Base.NamedTuple(value::Formula)
    return (identifier=formula_id(value), assumptions=value.assumptions,
        parameters=value.parameters, hooks=value.hooks, options=value.options,
        equivalent_earth=value.equivalent_earth === nothing ? nothing : NamedTuple(value.equivalent_earth))
end
