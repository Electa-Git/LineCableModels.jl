function description(::Type{<:Formula{:unified}}; compact::Bool = false)
    compact ? "Unified" :
    "Unified circumferential earth impedance with complete enclosed-current normalization"
end

function description(::Type{<:Formula{:unified}}, ::Val{:Γ}, value::Number; compact::Bool = false)
    "Γ="*string(value)*" m⁻¹"
end
function description(::Type{<:Formula{:unified}}, ::Val{:Γ}, value::AbstractVector; compact::Bool = false)
    "Γ=["*join(value, ", ")*"] m⁻¹ (frequency order)"
end

"""
$(TYPEDSIGNATURES)

Calculate the circumferentially averaged axial-field coefficient per source
current \\[Ω/m\\] in two homogeneous half-spaces. The exp(jωt) convention and
caller-prescribed Γ \\[1/m\\] give

```math
\\widetilde K_{ij}=\\mathcal Z_{ij}-\\Gamma^2\\mathcal P_{\\phi,ij}/(j\\omega).
```

Source columns include exp(abs(real(κⱼrⱼ))) scaling, shared with the source-potential
and enclosed-current matrices. `pair` retains source-target geometry \\[m\\].
`functor` contains evaluated media and circumferential factors. The complete
current map converts these coefficients to physical series impedance.

The electric scalar-potential contribution is retained when Γ is nonzero.
Quadrature estimates remain diagnostic warnings.

# Returns

- Scaled axial-field coefficient \\[Ω/m\\].

# Reference

User-supplied manuscript, *Unified circumferentially averaged framework for
overhead, buried, and mixed conductor systems*, complete-field current relation.
"""
function axial_field_coefficient(
        ::Union{Formula{:unified}, Val{:unified}}, kind::Union{Val{:self}, Val{:mutual}},
        source::Union{Val{1}, Val{2}}, target::Union{Val{1}, Val{2}}, functor, pair, workspace)
    u=functor.state
    medium=target === Val(1) ? 1 : 2
    hp, hq=abs(pair.heights[2]), abs(pair.heights[1])
    row, column=pair.row, pair.column
    r=u.radius[row]
    average=u.circumference_average[row]
    sp, sq=u.source_logscale[row], u.source_logscale[column]
    πT=one(u.jω)*π
    integration=functor.options.data.integration
    context=(
        formula = :unified, frequency = imag(u.jω)/(2π), receiver = row, source = column)
    direct=earth_direct(
        kind, source, target, u, pair, r, average, u.radial_argument[row], sp, sq)
    z=u.jω/πT*average*earth_spectral_term(Val(:Z), target, source, u,
        hp, hq, pair.separation, zero(r), sp+sq,
        integration.method, integration.options, workspace.buffers; context)
    z+=u.jω*u.mu[medium]/(2πT)*direct
    phi=zero(z)
    if !iszero(u.Γ)
        phi=u.jω/πT*average*earth_spectral_term(Val(:phi), target, source, u,
            hp, hq, pair.separation, zero(r), sp+sq,
            integration.method, integration.options, workspace.buffers; context)
        phi+=u.jω/(2πT*u.sh[medium])*direct
    end
    return z-u.Γ^2/u.jω*phi
end

function Expression(selected::Formula{:unified}, pair::EarthPair)
    return Expression(selected, axial_field_coefficient,
        Val(pair.row == pair.column ? :self : :mutual), Val.(layer_index(pair))...)
end

function formulation_options(::Expression{<:Formula{:unified},
        typeof(axial_field_coefficient),
        A}) where {
        A <: Tuple{
        Union{Val{:self}, Val{:mutual}}, Union{Val{1}, Val{2}}, Union{Val{1}, Val{2}}}}
    return FormulationOptions((Γ = 0, integration = (method = :quad, options = (;))))
end

function validate(reduction::EquivalentHomogeneous.Formula{:bottommost},
        ::Expression{<:Formula{:unified}, typeof(axial_field_coefficient)})
    return reduction
end

:unified
