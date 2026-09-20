include("unified/kernels.jl")
include("unified/current.jl")

function assumptions(::Val{:unified})
    (media = :homogeneous, layers = 2:2, permittivity = :positive)
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

Here s=jω uses the ``e^{+j\\omega t}`` phasor convention, Ze has units \\[Ω/m\\], Pe has units \\[m/F\\], and Ye has units
\\[S/m\\]. The shared implementation forms physical matrices before selecting
entries. Air receivers use the interface voltage; earth receivers use deep-earth
voltage. These references are fixed by receiving layer, not selectable parameters.

Prescribe Γ \\[1/m\\] through `formula(:unified; options=(Γ=value,))`.
In the convention ``\\widehat F_- e^{\\Gamma_- x-j\\omega t}``, the same real
field has positive-time phasor ``\\widehat F_+=\\overline{\\widehat F_-}`` and
``\\Gamma_+=\\overline{\\Gamma_-}``. Consistent convention conversion also
conjugates complex material coefficients and impedance/admittance phasors.
This evaluator uses the supplied Γ verbatim in its positive-time convention;
it performs no automatic conjugation, sign change or modal iteration. Γ=0
removes longitudinal dependence, not frequency dependence.

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
    return FormulationOptions((Γ = 0, integration = (method = :quad, options = (;))))
end

function validate(binding::FormulaMethod{<:Formula{:unified}, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:bottommost})
    return binding
end

:unified
