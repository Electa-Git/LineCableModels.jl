function assumptions(::Val{:Wise1948})
    (media = :homogeneous, layers = 2:2, longitudinal = :zero, permittivity = :positive)
end

"""
$(TYPEDSIGNATURES)

**Identification.** Wideband homogeneous-earth overhead potential
coefficient.

**Expression.**

```math
P_{e,ij}=\\frac{P_{0,ij}+M_{ij}+jN_{ij}}{2\\pi\\varepsilon_0},
```

```math
M_{ij}+jN_{ij}=2\\int_0^\\infty
\\frac{e^{-H\\lambda}\\cos(y_{ij}\\lambda)}
{(\\gamma_1^2/\\gamma_0^2)\\lambda+
\\sqrt{\\lambda^2+\\gamma_1^2-\\gamma_0^2}}d\\lambda,
\\quad P_{0,ij}=\\ln(D_{ij}/d_{ij}).
```

**Reference.** W. H. Wise, “Potential Coefficients for Ground Return
Circuits,” *Bell System Technical Journal*, 27, 365–371, 1948.
"""
function description(::Formula{:Wise1948})
    "Wise1948 homogeneous-earth overhead potential coefficient"
end

Γ(::Val{:Wise1948}, jω, materials, layers) = zero(jω)

raw"""
Evaluate Wise's wideband overhead earth potential coefficient:

```math
P_{e,ij}=\frac{P_{0,ij}+M_{ij}+jN_{ij}}{2\pi\varepsilon_0},
```

```math
M_{ij}+jN_{ij}=2\int_0^\infty
\frac{e^{-(h_i+h_j)\lambda}\cos(y_{ij}\lambda)}
{(\gamma_1^2/\gamma_0^2)\lambda+
\sqrt{\lambda^2+\gamma_1^2-\gamma_0^2}}d\lambda,
\qquad P_{0,ij}=\ln(D_{ij}/d_{ij}).
```
"""
function earth_potential_coefficient(
        ::Val{:Wise1948}, ::Union{Val{:self}, Val{:mutual}}, ::Val{1}, ::Val{1},
        functor, pair, workspace
)
    state = functor.state
    geometry = _geometry(pair)
    gamma_0_squared, gamma_1_squared = state.gamma_medium_squared
    ratio = gamma_1_squared / gamma_0_squared
    integral = integrate(functor.options.integration.method,
        SpectralIntegral(Val(:cosine),
            lambda -> begin
                radial = sqrt(lambda^2 + gamma_1_squared - gamma_0_squared)
                origin = sqrt(gamma_1_squared - gamma_0_squared)
                # Exact pole extraction; rationalization avoids subtracting close kernels.
                -lambda^2 /
                ((origin + radial) * (ratio*lambda + radial) * (ratio*lambda + origin))
            end,
            (height = geometry.H, separation = geometry.y_ij),
            float(nominal(abs(sqrt(gamma_1_squared - gamma_0_squared))));
            pole = (residue = inv(ratio),
                location = sqrt(gamma_1_squared-gamma_0_squared)/ratio),
            angle = min(pi/4, atan(float(nominal(geometry.H / (2geometry.y_ij))))),
            features = retained_earth_features(state, geometry, workspace)),
        functor.options.integration.options, workspace)
    return (log(geometry.D_ij / geometry.d_ij) + 2 * integral) /
           (2π * state.epsilon[1])
end

Formulation(::LineCableModelsCoaxial, selected::Formula{:Wise1948}) = selected

function hooks(::FormulaMethod{:Wise1948, typeof(earth_potential_coefficient),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{1}, Val{1}}}
    return (configurable = (:Γ, :air, :earth, :permeability, :contribution),
        defaults = (
            Γ = FormulaMethod(Val(:Wise1948), Γ),
            air = FormulaMethod(Val(:full), propagation),
            earth = FormulaMethod(Val(:full), propagation),
            permeability = vacuum_permeability,
            contribution = nothing))
end

function computation_options(::FormulaMethod{
        :Wise1948, typeof(earth_potential_coefficient),
        A}) where {A <: Tuple{Union{Val{:self}, Val{:mutual}}, Val{1}, Val{1}}}
    (integration = (method = :quad, options = (;)),)
end

function validate(binding::FormulaMethod{:Wise1948, typeof(earth_potential_coefficient)},
        ::EquivalentHomogeneous.Formula{:default})
    binding
end

:Wise1948
