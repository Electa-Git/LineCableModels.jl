include("unified/kernels.jl")
include("unified/current.jl")

function assumptions(::Val{:unified})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
end

"""
Normalize Unified's prescribed longitudinal argument Γ [1/m]. A scalar applies
at every frequency; a vector follows the computation's frequency order exactly.
This is a prescribed argument, not a solved propagation mode.
"""
function earth_parameters(::Val{:unified}, parameters::NamedTuple)
    all(==(:Γ), keys(parameters)) || throw(ArgumentError(
        "unified admits only the physical parameter Γ [1/m]"))
    argument = get(parameters, :Γ, 0)
    values = argument isa Number ? (argument,) : argument
    values isa Union{Tuple, AbstractVector} && !isempty(values) &&
    all(value -> value isa Number && !(value isa Bool) && isfinite(value), values) ||
        throw(ArgumentError("unified Γ must be a finite scalar or nonempty finite vector [1/m]"))
    argument isa Union{Number, AbstractVector} || throw(ArgumentError(
        "unified Γ must be a scalar or frequency-aligned vector [1/m]"))
    return (; Γ = argument isa AbstractVector ? copy(argument) : argument)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Circumferentially averaged, two-half-space earth return
with the complete enclosed-current constraint. Rows are receivers and columns
are sources. Air, earth and both ordered mixed directions are supported, with
independent permeabilities and caller-prescribed Γ (zero by default).

**Expression.** Assemble the source kernels and solve the complete system:

```math
P_e L=H,\\qquad Z_e L=K+\\Gamma^2 H/s,\\qquad Y_e H=sL,
\\quad L=A_r^{-1}-F_rK,\\quad K=\\mathcal Z-\\Gamma^2\\mathcal P_\\phi/s.
```

Here s=jω, Ze has units \\[Ω/m\\], Pe has units \\[m/F\\], and Ye has units
\\[S/m\\]. The shared implementation forms physical matrices before selecting
entries. Air receivers use the interface voltage; earth receivers use deep-earth
voltage. These references are fixed by receiving layer, not selectable parameters.

**Reference.** User-supplied manuscript, *Unified circumferentially averaged
framework for overhead, buried, and mixed conductor systems*. Current closure, reference conventions and the
equal-medium limit are exercised in `test/unit/engine/unified_earth_return.jl`.
These scoped controls do not establish acceptance of arbitrary complete matrices.
"""
function description(::Type{<:Formula{:unified}}; compact::Bool = false)
    compact ? "Unified" :
    "Unified circumferential earth impedance with full current closure"
end

function earth_bindings(selected::Formula{:unified},
        physical::AbstractVector{<:EarthPair}, homogeneous, indices)
    binding = invoke(earth_bindings,
        Tuple{EarthImpedanceFormulation, AbstractVector{<:EarthPair}, Any, Any},
        selected, physical, homogeneous, collect(eachindex(physical)))
    _unified_geometry(physical)
    options = first(binding.equations).declaration.options
    all(g -> isequal(g.declaration.options, options), binding.equations) ||
        throw(ArgumentError("the unified current closure requires common integration controls"))
    equations = [(declaration = group.declaration,
                     indices = intersect(group.indices, indices))
                 for group in binding.equations if !isdisjoint(group.indices, indices)]
    return merge(binding, (; equations))
end

function earth!(Z, P, selected::Formula{:unified}, ::Nothing,
        binding, ::Nothing, materials, workspace, frequency)
    _unified_current!(workspace, materials, binding, frequency)
    earth!(Z, selected, binding, materials, workspace.input.jω[frequency],
        workspace, materials.thickness)
    return workspace
end

function earth!(destination, ::Formula{:unified},
        group::NamedTuple{(:declaration, :indices)}, interactions,
        materials, jω, workspace, thickness)
    for index in group.indices
        pair=interactions[index].pair
        value=group.declaration.equation(nothing, pair, workspace)
        isfinite(value) ||
            throw(DomainError(value, "earth-impedance contribution must be finite"))
        destination[pair.row, pair.column]=oftype(jω, value)
    end
    return destination
end

function earth_impedance(::Formula{:unified}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{1}, ::Val{1}, functor, pair, workspace)
    workspace === nothing && throw(ArgumentError(
        "unified requires the complete current system; use compute for physical entries"))
    return workspace.buffers.unified.Ze[pair.row, pair.column]
end

function earth_impedance(::Formula{:unified}, ::Union{Val{:self}, Val{:mutual}},
        ::Val{2}, ::Val{2}, functor, pair, workspace)
    workspace === nothing && throw(ArgumentError(
        "unified requires the complete current system; use compute for physical entries"))
    return workspace.buffers.unified.Ze[pair.row, pair.column]
end

function earth_impedance(::Formula{:unified}, ::Val{:mutual},
        ::Val{1}, ::Val{2}, functor, pair, workspace)
    workspace === nothing && throw(ArgumentError(
        "unified requires the complete current system; use compute for physical entries"))
    return workspace.buffers.unified.Ze[pair.row, pair.column]
end

function earth_impedance(::Formula{:unified}, ::Val{:mutual},
        ::Val{2}, ::Val{1}, functor, pair, workspace)
    workspace === nothing && throw(ArgumentError(
        "unified requires the complete current system; use compute for physical entries"))
    return workspace.buffers.unified.Ze[pair.row, pair.column]
end

function formulation_options(::FormulaMethod{<:Formula{:unified},
        typeof(earth_impedance),
        A}) where {
        A <: Tuple{
        Union{Val{:self}, Val{:mutual}}, Union{Val{1}, Val{2}}, Union{Val{1}, Val{2}}}}
    return FormulationOptions((integration = (method = :quad, options = (;)),))
end

function validate(binding::FormulaMethod{<:Formula{:unified}, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:bottommost})
    return binding
end

:unified
