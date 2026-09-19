"""
$(TYPEDEF)

Select cylindrical surface equations and their model and numerical controls.
The actual conductor geometry determines which surfaces are required.

$(TYPEDFIELDS)
"""
struct Formula{ID, P <: NamedTuple, O <: FormulationOptions} <: InternalImpedanceFormulation
    "Resolved physical/model parameters."
    parameters::P
    "Normalized numerical sections indexed by surface kind."
    options::O
end

"""
$(TYPEDEF)

Retain a selected formulation and shared state for one conductor and frequency.
Leaf evaluation uses that selection's equation directly.

$(TYPEDFIELDS)
"""
struct Functor{F, S, O <: FormulationOptions}
    "Selected internal-impedance formulation."
    selection::F
    "Shared conductor state at one frequency."
    state::S
    "Normalized numerical sections for the surfaces or the current leaf."
    options::O
end

"""Return the stable identifier of a selected internal-impedance formulation."""
formula_id(::Formula{ID}) where {ID} = ID

"""Evaluate a selected cylindrical surface coefficient in Ω/m."""
function internal_impedance end

"""Evaluate required cylindrical surface impedances in Ω/m."""
function surface_impedances end

"""
$(TYPEDSIGNATURES)

Construct an internal-impedance formulation. Numerical controls are projected
onto its surface equations; unknown numerical sections are rejected.
Custom formulations subtype `InternalImpedanceFormulation`, supply their
shared-state constructor, and extend `internal_impedance` on their own type.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

function Formula(::Val{ID}; parameters::NamedTuple=(;), options::Union{NamedTuple, FormulationOptions} = FormulationOptions()) where {ID}
    ID in formulas() || throw(ArgumentError("unknown internal-impedance formula :$ID"))
    options = options isa NamedTuple ? FormulationOptions(options) : options
    isempty(parameters) || throw(ArgumentError("internal impedance :$ID has no model parameters"))
    kinds = (:inner, :outer, :transfer)
    supplied = keys(options.data)
    selected = Formula{ID, typeof(parameters), typeof(options)}(parameters, options)
    declarations = map(kinds) do kind
        binding = FormulaMethod(selected, internal_impedance, Val(kind))
        formulation_options(binding)
    end
    admitted = union((keys(value.data) for value in declarations)...)
    isempty(setdiff(supplied, admitted)) || throw(ArgumentError(
        "unused internal-impedance numerical sections"))
    normalized = FormulationOptions(NamedTuple{kinds}(map(kinds, declarations) do kind, defaults
        binding = FormulaMethod(selected, internal_impedance, Val(kind))
        names = Tuple(intersect(supplied, keys(defaults.data)))
        formulation_options(binding, defaults, FormulationOptions(options.data[names])).data
    end))
    return Formula{ID, typeof(parameters), typeof(normalized)}(parameters, normalized)
end

Formula(selected::InternalImpedanceFormulation) = selected

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula(Val(ID); parameters=selection.parameters, options=selection.options)
end

function internal_impedance(selected::InternalImpedanceFormulation, ::Val{Kind},
        functor, workspace) where {Kind}
    throw(ArgumentError(
        "internal_impedance :$(formula_id(selected)): formula not implemented for kind :$Kind"))
end

function (functor::Functor)(::Val{Kind}, workspace=nothing) where {Kind}
    options = get(functor.options.data, Kind, (;))
    leaf = Functor(functor.selection, functor.state, FormulationOptions(options))
    value = internal_impedance(functor.selection, Val(Kind), leaf, workspace)
    value isa Number && isfinite(value) || throw(DomainError(value,
        "internal_impedance must return a finite surface coefficient [Ω/m]"))
    return value
end

(functor::Functor)(kind::Symbol, workspace=nothing) = functor(Val(kind), workspace)

"""
$(TYPEDSIGNATURES)

Evaluate cylindrical surface coefficients in Ω/m. The wall operator is
`[inner transfer; transfer outer]` in the surface-current basis
`(-enclosed axial current, total axial current including the wall)`.
Assemblers own basis transformation and matrix placement.

# Arguments

- `formula`: One formulation or an explicit surface recipe covering the primitive.
- `r_in`, `r_ex`: Inner and outer conductor radii \\[m\\].
- `rho`: Conductor resistivity \\[Ω·m\\].
- `mu_r`: Relative permeability \\[dimensionless\\].
- `jω`: Imaginary angular frequency \\[1/s\\].
- `workspace`: Optional computation workspace passed to each selected equation.

# Returns

- NamedTuple of `outer` for a solid primitive, or `inner`, `outer`, and
  `transfer` for a tubular primitive \\[Ω/m\\].
"""
function surface_impedances(formula::Union{InternalImpedanceFormulation,
        NamedTuple}, r_in, r_ex, rho, mu_r, jω;
        workspace=nothing)
    return surface_impedances(formula,
        r_in > 0 ? Val((:inner,:outer,:transfer)) : Val((:outer,)),
        r_in, r_ex, rho, mu_r, jω; workspace)
end

"""Evaluate required surfaces, preparing shared state once per conductor and frequency."""
@inline function surface_impedances(formula::InternalImpedanceFormulation, ::Val{Kinds},
        r_in, r_ex, rho, mu_r, jω; workspace=nothing) where {Kinds}
    functor = formula(r_in, r_ex, rho, mu_r, jω)
    length(Kinds) == 1 && return NamedTuple{Kinds}((functor(Val(Kinds[1]),workspace),))
    length(Kinds) == 3 || throw(ArgumentError("internal surfaces require one or three kinds"))
    return NamedTuple{Kinds}((functor(Val(Kinds[1]),workspace),
        functor(Val(Kinds[2]),workspace), functor(Val(Kinds[3]),workspace)))
end

@inline function surface_impedances(selected::NamedTuple,
        ::Val{Kinds}, r_in, r_ex, rho, mu_r, jω; workspace=nothing) where {Kinds}
    if length(Kinds) == 1
        return surface_impedances(Formulation(selected, Val(first(Kinds))), Val(Kinds),
            r_in, r_ex, rho, mu_r, jω; workspace)
    end
    length(Kinds) == 3 || throw(ArgumentError("internal surfaces require one or three kinds"))
    a, b, c = map(kind -> Formulation(selected, Val(kind)), Kinds)
    first_functor = a(r_in, r_ex, rho, mu_r, jω)
    second_functor = b === a ? first_functor : b(r_in, r_ex, rho, mu_r, jω)
    third_functor = c === a ? first_functor : c === b ? second_functor :
                    c(r_in, r_ex, rho, mu_r, jω)
    return NamedTuple{Kinds}((first_functor(Val(Kinds[1]),workspace),
        second_functor(Val(Kinds[2]),workspace), third_functor(Val(Kinds[3]),workspace)))
end

"""Expose the selected identity, model parameters and numerical controls."""
Base.NamedTuple(value::Formula) = (identifier=formula_id(value),
    parameters=value.parameters, options=value.options.data)

import ...Grammar: formulation_options
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((inner=Formula, outer=Formula, transfer=Formula))
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
formulation_options(value::Formula) = value.options
