function assumptions(::Val{:unified})
    (media = :homogeneous, layers = 2:2,
        longitudinal = :prescribed, permittivity = :positive)
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
entries. Pair evaluation requires a prepared full-system workspace. The default
voltage reference is deep earth; `:interface`, `:scalar`, or a positive reference
depth \\[m\\] are selected through the `reference` physical parameter.

**Reference.** User-supplied manuscript, *Unified circumferentially averaged
framework for overhead, buried, and mixed conductor systems*. Current closure, reference conventions and the
equal-medium limit are exercised in `test/unit/engine/unified_earth_return.jl`.
These scoped controls do not establish acceptance of arbitrary complete matrices.
"""
function description(::Type{<:Formula{:unified}}; compact::Bool=false)
    compact ? "Unified" : "Unified circumferential earth potential with full current closure"
end

system_earth(::Formula{:unified}) = true

function earth_potential_coefficient(::Formula{:unified}, ::Union{Val{:self}, Val{:mutual}},
        ::Union{Val{1}, Val{2}}, ::Union{Val{1}, Val{2}}, functor, pair, workspace)
    return unified_entry(Val(:potential), functor, pair, workspace)
end


function formulation_options(::FormulaMethod{<:Formula{:unified},
        typeof(earth_potential_coefficient),
        A}) where {
        A <: Tuple{
        Union{Val{:self}, Val{:mutual}}, Union{Val{1}, Val{2}}, Union{Val{1}, Val{2}}}}
    return FormulationOptions((integration = (method = :quad, options = (;)),))
end

function validate(binding::FormulaMethod{<:Formula{:unified}, typeof(earth_potential_coefficient)},
        ::EquivalentHomogeneous.Formula{:bottommost})
    return binding
end

"""
$(TYPEDSIGNATURES)

Evaluate this formulation's medium state: absolute permeability \\[H/m\\]
and transverse propagation constant \\[1/m\\]. Air retains its prescribed
permeability; soil follows the selected source's magnetic approximation.
"""
function constitutive(::Formula{:unified}, ::Val{:air}, jω, μ, σ, ε)
    return (mu=μ, gamma=propagation(Val(:full), jω, μ, σ, ε))
end

function constitutive(::Formula{:unified}, ::Val{:earth}, jω, μ, σ, ε)
    permeability = μ
    return (mu=permeability, gamma=propagation(Val(:full), jω, permeability, σ, ε))
end

:unified
