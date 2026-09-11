function assumptions(::Val{:default})
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
entries. Pair callbacks require a prepared full-system workspace. The default
voltage reference is deep earth; `:interface`, `:scalar`, or a positive reference
depth \\[m\\] are selected through the `reference` physical parameter.

**Reference.** User-supplied manuscript, *Unified circumferentially averaged
framework for overhead, buried, and mixed conductor systems*. The mathematical
source hash and accepted complete matrices are recorded in
`test/fixtures/reference/unified_earth_return.toml` and the locked implementation
plan in `docs/notes/unified-earth-return-implementation-plan.md`.
"""
function description(::Formula{:default})
    "Unified circumferential earth impedance with full current closure"
end

Γ(::Val{:default}, jω, materials, layers) = zero(jω)
system_earth(::Formula{:default}) = true
Formulation(::LineCableModelsCoaxial, selected::Formula{:default}) = selected

function earth_impedance(::Val{:default}, ::Union{Val{:self}, Val{:mutual}},
        ::Union{Val{1}, Val{2}}, ::Union{Val{1}, Val{2}}, functor, pair, workspace)
    return unified_entry(Val(:impedance), functor, pair, workspace)
end

function hooks(::FormulaMethod{:default,
        typeof(earth_impedance),
        A}) where {
        A <: Tuple{
        Union{Val{:self}, Val{:mutual}}, Union{Val{1}, Val{2}}, Union{Val{1}, Val{2}}}}
    return (configurable = (:Γ, :air, :earth, :permeability, :contribution),
        defaults = (Γ = FormulaMethod(Val(:default), Γ),
            air = FormulaMethod(Val(:full), propagation),
            earth = FormulaMethod(Val(:full), propagation),
            permeability = identity, contribution = nothing))
end

function computation_options(::FormulaMethod{:default,
        typeof(earth_impedance),
        A}) where {
        A <: Tuple{
        Union{Val{:self}, Val{:mutual}}, Union{Val{1}, Val{2}}, Union{Val{1}, Val{2}}}}
    return (integration = (method = :quad, options = (;)),)
end

function validate(binding::FormulaMethod{:default, typeof(earth_impedance)},
        ::EquivalentHomogeneous.Formula{:default})
    return binding
end

:default
