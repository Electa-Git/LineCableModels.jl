"""
$(TYPEDEF)

Own one earth potential coefficient formulation, its indexed equation declarations and explicit
controls. `assumptions` contains scientific restrictions; each interaction
is selected by its explicit kind/source-layer/target-layer equation signature.
`parameters` stores model controls; `options` stores numerical controls.

$(TYPEDFIELDS)
"""
struct Formula{ID, A <: NamedTuple, P <: NamedTuple, O <: FormulationOptions, E} <:
       EarthAdmittanceFormulation
    "Model class, exact medium inventory and longitudinal restriction."
    assumptions::A
    "Explicit physical/model parameters."
    parameters::P
    "Explicit numerical sections; projected onto required indexed equations at preflight."
    options::O
    "Independent equivalent homogeneous-earth reduction, or nothing."
    equivalent_earth::E
end

"""
$(TYPEDEF)

Bind one indexed interaction to evaluated physical state. Numerical options and
physical quantities remain distinct. Mutable numerical
resources are passed when evaluating the functor.

$(TYPEDFIELDS)
"""
struct Functor{B, S, O <: FormulationOptions}
    "Selected equation and its exact source/target geometry."
    binding::B
    "Evaluated physical quantities, including one authoritative Γ [1/m]."
    state::S
    "Normalized numerical controls of the selected equation."
    options::O
end

formula_id(::Formula{ID}) where {ID} = ID
assumptions(formula::EarthAdmittanceFormulation) = formula.assumptions
media(formula::EarthAdmittanceFormulation) = formula.assumptions.media

"""
Declare model inventory and physical restrictions, without callable behavior.
"""
function assumptions end

assumptions(::Val{ID}) where {ID} = throw(ArgumentError("unknown earthadmittance formula :$ID"))
"""
Evaluate a medium propagation law at jω [1/s], μ [H/m], σ [S/m], ε [F/m].
"""
function propagation end
"""
Read the evaluated longitudinal propagation constant Γ [1/m].
"""
function Γ end
"""
Evaluate one source-owned earth potential coefficient equation.
"""
function earth_potential_coefficient end

Γ(functor::Functor) = functor.state.Γ

Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

"""
$(TYPEDSIGNATURES)

Resolve a formulation's model parameters and numerical controls. Its concrete
selection type owns the indexed equation and its medium constitutive state.
Longitudinal propagation is prescribed by the problem, not by a replacement
channel. A missing physical case is unsupported.
"""
function Formula(::Val{ID}; parameters::NamedTuple = (;),
        options::Union{NamedTuple, FormulationOptions} = FormulationOptions(), equivalent_earth = nothing) where {ID}
    options = options isa NamedTuple ? FormulationOptions(options) : options
    parameters = earth_parameters(Val(ID), parameters)
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
    return Formula{ID, typeof(constraints), typeof(parameters),
        typeof(options), typeof(reduction)}(
        constraints, parameters, options, reduction)
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

- `Γ`: Optional prescribed longitudinal constant [1/m]; zero when omitted.
- `thickness`: Aligned layer thicknesses [m] for a stratified model.
"""
function (formula::EarthAdmittanceFormulation)(
        resistivity::AbstractVector{T}, permittivity::AbstractVector{T},
        permeability::AbstractVector{T}, jω::Complex{T}, pair::EarthPair;
        Γ = nothing, thickness = nothing, physical_pair = pair
) where {T <: Real}
    selected = validate(formula, pair)
    return formula(resistivity, permittivity, permeability, jω, pair, selected;
        Γ, thickness, physical_pair)
end

# Evaluate a declaration already bound by validate() before the frequency loop.
function (formula::EarthAdmittanceFormulation)(
        resistivity::AbstractVector{T}, permittivity::AbstractVector{T},
        permeability::AbstractVector{T}, jω::Complex{T}, pair::EarthPair, selected;
        Γ = nothing, thickness = nothing, physical_pair = pair
) where {T <: Real}
    isfinite(jω) && !iszero(jω) || throw(DomainError(jω, "jω must be finite and nonzero"))
    options = selected.options
    validate(formula, resistivity, permittivity, permeability, thickness)
    thickness === nothing || validate(pair, thickness)
    layers = media(formula) === Val(:homogeneous) ? (1, 2) : eachindex(permeability)
    σ = map(layer -> conductivity(resistivity[layer]), layers)
    evaluated = map(layers) do layer
        medium = constitutive(formula, layer == 1 ? Val(:air) : Val(:earth),
            jω, permeability[layer], σ[layer], permittivity[layer])
        medium.mu isa Real && isfinite(medium.mu) && medium.mu > 0 ||
            throw(DomainError(medium.mu, "medium permeability must be positive and finite [H/m]"))
        medium.gamma isa Number && isfinite(medium.gamma) ||
            throw(DomainError(medium.gamma, "medium propagation must be finite [1/m]"))
        (mu=convert(T, medium.mu), gamma=convert(Complex{T}, medium.gamma))
    end
    μ = map(medium -> medium.mu, evaluated)
    γ = map(medium -> medium.gamma, evaluated)
    materials = (rho = resistivity, epsilon = permittivity, mu = μ, sigma = σ,
        gamma = γ, gamma_medium_squared = γ .^ 2)
    longitudinal = Γ === nothing ? zero(jω) : Γ
    longitudinal isa Number && isfinite(longitudinal) || throw(ArgumentError(
        "Γ must be one finite scalar [1/m], not a value/square pair"))
    formula.assumptions.longitudinal === :zero && !iszero(longitudinal) &&
        throw(ArgumentError("earth-potential coefficient :$(formula_id(formula)) fixes Γ to zero"))
    state = merge(materials, (; jω, Γ = oftype(jω, longitudinal), thickness))
    binding = (pair = pair, physical_pair = physical_pair,
        kind = selected.kind, equation = selected.equation)
    return Functor(binding, state, options)
end

function (functor::Functor)(workspace = nothing)
    pair = functor.binding.pair
    value = functor.binding.equation(functor, pair, workspace)
    value isa Number && isfinite(value) || throw(DomainError(value,
        "earth-potential coefficient contribution must be a finite scalar"))
    return oftype(functor.state.jω, value)
end


function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    return Formula(Val(ID); parameters = selection.parameters,
        options = selection.options, equivalent_earth = selection.equivalent_earth)
end

Formula(selected::EarthAdmittanceFormulation) = selected

function FormulaMethod(formula::EarthAdmittanceFormulation, pair::EarthPair)
    return FormulaMethod(formula, earth_potential_coefficient,
        Val(pair.row == pair.column ? :self : :mutual), Val.(pair.layers)...)
end

function earth_potential_coefficient(selected::EarthAdmittanceFormulation, ::Val{Kind}, ::Val{S}, ::Val{T},
        functor, pair, workspace) where {Kind, S, T}
    throw(ArgumentError(
        "earth_potential_coefficient :$(formula_id(selected)) ($Kind): formula not implemented for source in layer $S and target in layer $T"))
end

const EQUATION_FALLBACK = which(earth_potential_coefficient, Tuple{EarthAdmittanceFormulation, Val, Val, Val, Any, Any, Any})

function validate(binding::FormulaMethod{
        ID, typeof(earth_potential_coefficient), A}) where {ID, A}
    signature = Tuple{ID, A.parameters..., Any, Any, Any}
    if which(earth_potential_coefficient, signature) === EQUATION_FALLBACK
        # The diagnostic takes no physical inputs and evaluates no numerical kernel.
        binding(nothing, nothing, nothing)
    end
    return binding
end

"""
$(TYPEDSIGNATURES)

Expose the selected identity, scientific restrictions, model and numerical
controls, and explicit equivalent-earth reduction as a native record.
"""
function Base.NamedTuple(value::Formula)
    return (identifier=formula_id(value), assumptions=value.assumptions,
        parameters=value.parameters, options=value.options.data,
        equivalent_earth=value.equivalent_earth === nothing ? nothing : NamedTuple(value.equivalent_earth))
end

# Identity-only dispatch also describes retained selections without constructors.
import ...Grammar: formulation_options
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)

"""Iterate the independently selectable child slots admitted by this formula family."""
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((air=Formula, earth=Formula, mixed=Formula))
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
formulation_options(value::Formula) = value.options
