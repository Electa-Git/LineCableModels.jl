"""
$(TYPEDEF)

Select an internal-impedance recipe by its registered identifier.

Surface equations receive `(functor, workspace)`. Physical parameters, explicit
callables and numerical options have separate fields.

$(TYPEDFIELDS)
"""
struct Formula{
    ID,
    R <: NamedTuple,
    A <: NamedTuple,
    H <: NamedTuple,
    O <: NamedTuple,
    C <: Tuple
} <: InternalImpedanceFormulation
    "Surface equation bindings; an unimplemented kind has no binding."
    binding::R
    "Physical assumptions of the selected recipe."
    parameters::A
    "Surface-equation overrides supplied by the user."
    hooks::H
    "Numerical sections indexed by surface kind."
    options::O
    "Explicitly configured numerical section names, checked against actual consumers."
    configured_options::C
end

"""
$(TYPEDEF)

Store the values shared by the leaf interactions of one formula call.

The state has no common physical layout: every formula owns its state and call
methods.

$(TYPEDFIELDS)
"""
struct Functor{ID, B, H, S, O}
    "Declared equation bindings, or the single bound surface equation."
    binding::B
    "Concrete callable overrides retained unchanged from selection."
    hooks::H
    "Formula-owned shared numerical state."
    state::S
    "Normalized numerical options for the bound kind or surface collection."
    options::O
end

"""
Return the stable formula identifier of an internal-impedance formula.
"""
formula_id(::Formula{ID}) where {ID} = ID

"""
Evaluate one formula-owned internal-impedance interaction.
"""
function internal_impedance end

"""
Return the three cylindrical surface impedances supplied by one formula.
"""
function surface_impedances end

"""
$(TYPEDSIGNATURES)

Construct an internal-impedance formula from a registered identifier.

The `hooks` record replaces individual surface equations. Each replacement
declares its numerical defaults through `computation_options(binding, callable)`.
"""
Formula(identifier::Symbol; kwargs...) = Formula(Val(identifier); kwargs...)

function Formula(::Val{ID}; parameters::NamedTuple = (;),
        hooks::NamedTuple = (;), options::NamedTuple = (;)) where {ID}
    haskey(hooks, :mutual) && throw(ArgumentError(
        "internal-impedance hooks use transfer, not mutual"))
    ID in FORMULAS || throw(ArgumentError("unknown internal-impedance formula :$ID"))
    isempty(parameters) ||
        throw(ArgumentError("internal impedance :$ID has no model parameters"))
    kinds = Tuple(key for (key,_) in pairs(Formula))
    equations = (inner = FormulaMethod(Val(ID), internal_impedance, Val(:inner)),
        outer = FormulaMethod(Val(ID), internal_impedance, Val(:outer)),
        transfer = FormulaMethod(Val(ID), internal_impedance, Val(:transfer)))
    Bindings = NamedTuple{kinds,
        Tuple{Union{Nothing, typeof(equations.inner)},
            Union{Nothing, typeof(equations.outer)}, Union{
                Nothing, typeof(equations.transfer)}}}
    defaults::Bindings = Bindings(map(kinds) do kind
        which(internal_impedance, Tuple{Val{ID}, Val{kind}, Any, Any}) === EQUATION_FALLBACK ?
        nothing : getproperty(equations, kind)
    end)
    all(isnothing, values(defaults)) &&
        throw(ArgumentError("internal impedance :$ID has no implemented surface equations"))
    all(kind -> get(defaults, kind, nothing) !== nothing, keys(hooks)) ||
        throw(ArgumentError("unknown or unimplemented internal-impedance hooks"))
    any(isnothing, values(hooks)) && throw(ArgumentError("surface hooks must be callable"))
    declarations = map(keys(defaults)) do kind
        getproperty(defaults, kind) === nothing && return (;)
        binding = getproperty(equations, kind)
        replacement = get(hooks, kind, nothing)
        replacement === nothing ? computation_options(binding) :
        computation_options(binding, replacement)
    end
    admitted = union((keys(value) for value in declarations)...)
    isempty(setdiff(keys(options), admitted)) ||
        throw(ArgumentError("unused internal-impedance numerical sections"))
    normalized = NamedTuple{keys(defaults)}(map(keys(defaults), declarations) do kind,
    declared
        names = Tuple(intersect(keys(options), keys(declared)))
        computation_options(getproperty(equations, kind), declared, options[names])
    end)
    return Formula{ID, typeof(defaults), typeof(parameters), typeof(hooks),
        typeof(normalized), typeof(keys(options))}(
        defaults, parameters, hooks, normalized, keys(options))
end

Formula(selected::Formula) = selected

function Formula(selection::FormulaDefinition{ID, Order}) where {ID, Order}
    Order === :default || throw(ArgumentError("order applies only to equivalent_earth"))
    selection.equivalent_earth === nothing || throw(ArgumentError(
        "equivalent_earth applies only to external earth formulas"))
    return Formula(Val(ID); parameters = selection.parameters, hooks = selection.hooks,
        options = selection.options)
end

function internal_impedance(::Val{ID}, ::Val{Kind}, functor, workspace) where {ID, Kind}
    throw(ArgumentError("internal_impedance :$ID: formula not implemented for kind :$Kind"))
end

const EQUATION_FALLBACK = which(internal_impedance, Tuple{Val, Val, Any, Any})

function validate(binding::FormulaMethod{ID, typeof(internal_impedance), A}) where {ID, A}
    if which(internal_impedance, Tuple{Val{ID}, A.parameters..., Any, Any}) ===
       EQUATION_FALLBACK
        binding(nothing, nothing)
    end
    return binding
end

function (functor::Functor{ID})(::Val{Kind}, workspace = nothing) where {ID, Kind}
    binding = FormulaMethod(Val(ID), internal_impedance, Val(Kind))
    get(functor.binding, Kind, nothing) === nothing && return binding(nothing, workspace)
    selected = get(functor.hooks, Kind, binding)
    options = getproperty(functor.options, Kind)
    case = Functor{
        ID, typeof(binding), typeof(functor.hooks), typeof(functor.state), typeof(options)}(
        binding, functor.hooks, functor.state, options)
    value = selected(case, workspace)
    value isa Number && isfinite(value) || throw(DomainError(value,
        "internal_impedance must return a finite surface coefficient [Ω/m]"))
    return value
end

(functor::Functor)(kind::Symbol, workspace = nothing) = functor(Val(kind), workspace)

"""
$(TYPEDSIGNATURES)

Evaluate cylindrical surface coefficients using a resolved formula, including
its explicit hooks and numerical options. The returned `(inner, outer, transfer)`
coefficients have units \\[Ω/m\\]. The wall operator is `[inner transfer; transfer outer]` in the surface-current
basis `(-enclosed axial current, total axial current including the wall)`.
Assemblers supply their current-basis transformation; this action performs no
matrix placement or enclosing-pipe calculation.

# Arguments

- `formula`: Resolved internal-impedance selection.
- `r_in`, `r_ex`: Inner and outer conductor radii \\[m\\].
- `rho`: Conductor resistivity \\[Ω·m\\].
- `mu_r`: Relative permeability \\[dimensionless\\].
- `jω`: Imaginary angular frequency \\[1/s\\].
- `workspace`: Optional numerical resources, passed unchanged to every kind.
"""
function surface_impedances(formula::Union{Formula,NamedTuple{(:inner,:outer,:transfer)}},
        r_in, r_ex, rho, mu_r, jω;
        workspace = nothing)
    validate(formula, (:inner, :outer, :transfer))
    return surface_impedances(formula, Val((:inner,:outer,:transfer)),
        r_in, r_ex, rho, mu_r, jω; workspace)
end

"""
$(TYPEDSIGNATURES)

Evaluate the prevalidated surface kinds in `Kinds`, in their supplied order,
returning a NamedTuple of coefficients \\[Ω/m\\]. Shared formula state is
prepared once per distinct complete selection for this conductor and frequency.
"""
@inline function surface_impedances(formula::Formula, ::Val{Kinds},
        r_in, r_ex, rho, mu_r, jω; workspace=nothing) where {Kinds}
    functor = formula(r_in, r_ex, rho, mu_r, jω)
    length(Kinds) == 1 && return NamedTuple{Kinds}((functor(Val(Kinds[1]),workspace),))
    length(Kinds) == 3 || throw(ArgumentError("internal surfaces require one or three kinds"))
    return NamedTuple{Kinds}((functor(Val(Kinds[1]),workspace),
        functor(Val(Kinds[2]),workspace), functor(Val(Kinds[3]),workspace)))
end

@inline function surface_impedances(selected::NamedTuple{(:inner,:outer,:transfer)},
        ::Val{Kinds}, r_in, r_ex, rho, mu_r, jω; workspace=nothing) where {Kinds}
    if length(Kinds) == 1
        return surface_impedances(selected[first(Kinds)], Val(Kinds),
            r_in, r_ex, rho, mu_r, jω; workspace)
    end
    length(Kinds) == 3 || throw(ArgumentError("internal surfaces require one or three kinds"))
    a, b, c = selected[Kinds[1]], selected[Kinds[2]], selected[Kinds[3]]
    first_functor = a(r_in, r_ex, rho, mu_r, jω)
    second_functor = b === a ? first_functor : b(r_in, r_ex, rho, mu_r, jω)
    third_functor = c === a ? first_functor : c === b ? second_functor :
                    c(r_in, r_ex, rho, mu_r, jω)
    return NamedTuple{Kinds}((first_functor(Val(Kinds[1]),workspace),
        second_functor(Val(Kinds[2]),workspace), third_functor(Val(Kinds[3]),workspace)))
end

"""Validate each required internal surface against its own selected formula and controls."""
function validate(selected::NamedTuple{(:inner,:outer,:transfer)}, kinds::Tuple)
    for (kind, leaf) in pairs(selected)
        if kind in kinds
            validate(leaf, (kind,))
        else
            isempty(leaf.hooks) && isempty(leaf.parameters) && isempty(leaf.configured_options) ||
                throw(ArgumentError("explicit internal $kind controls are unused by this assembly"))
        end
    end
    return selected
end

function validate(formula::Formula{ID}, kinds::Tuple) where {ID}
    :mutual in kinds && throw(ArgumentError("internal-impedance kinds use transfer, not mutual"))
    foreach(kinds) do kind
        validate(FormulaMethod(Val(ID), internal_impedance, Val(kind)))
    end
    isempty(setdiff(keys(formula.hooks), kinds)) || throw(ArgumentError(
        "an internal surface override is unused by the required conductor interactions"))
    admitted = union((keys(getproperty(formula.options, kind)) for kind in kinds)...)
    isempty(setdiff(formula.configured_options, admitted)) || throw(ArgumentError(
        "an explicitly configured numerical section is unused by the required internal surfaces"))
    return formula
end

"""
$(TYPEDSIGNATURES)

Expose the selected equation, physical parameters, callable overrides and numerical
options as a native record. Callables are retained unchanged.
"""
function Base.NamedTuple(value::Formula)
    return (identifier=formula_id(value), binding=value.binding,
        parameters=value.parameters, hooks=value.hooks, options=value.options,
        configured_options=value.configured_options)
end

# Identity-only dispatch also describes retained selections without constructors.
import ...Grammar: formulation_options
description(value::Formula; compact::Bool=false) = description(typeof(value); compact)

"""Iterate the independently selectable child slots admitted by this formula family."""
Base.pairs(::Type{<:Formula}; quantity=nothing) = pairs((inner=Formula, outer=Formula, transfer=Formula))
formula_id(::Type{<:Formula{ID}}) where {ID} = ID
formulation_options(value::Formula) = formulation_options(typeof(value), (parameters=value.parameters, hooks=value.hooks, options=value.options))
formulation_options(::Type{<:Formula}, retained::NamedTuple) =
    formulation_options(FormulaDefinition, retained)
