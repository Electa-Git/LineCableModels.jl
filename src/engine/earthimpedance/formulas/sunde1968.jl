function routes(identifier::Val{:Sunde1968})
    (
        self = FormulaMethod(identifier, earth_impedance, Val(:self)),
        mutual = FormulaMethod(identifier, earth_impedance, Val(:mutual)),
        homogeneous = FormulaMethod(
            identifier, earth_impedance, Val(:homogeneous)
        ),
        two_layer = FormulaMethod(
            identifier, earth_impedance, Val(:two_layer)
        ),
        multilayer = FormulaMethod(
            identifier, earth_impedance, Val(:multilayer)
        ),
        Γ = FormulaMethod(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Sunde1968})
    (
        air = _lossless,
        earth = _full,
        permeability = vacuum_permeability
    )
end

propagation(::Val{:Sunde1968}) = Val(:zero)
media(::Formula{:Sunde1968}) = Val(:stratified)
"""
$(TYPEDSIGNATURES)

**Identification.** Overhead impedance for homogeneous or horizontally
stratified earth; the number of earth layers selects the homogeneous,
two-layer, or recursive kernel.

**Expression.** The two-layer kernel is

```math
Z_{e,ij}=\\frac{j\\omega\\mu_0}{2\\pi}\\left[
\\ln\\frac{D_{ij}}{d_{ij}}+2\\int_0^\\infty F_{ij}^{S}(\\lambda)
\\cos(y_{ij}\\lambda)d\\lambda\\right],
```

```math
F_{ij}^{S}=\\frac{a_1+a_2+(a_1-a_2)e^{-2a_1d}}
{(a_1+a_2)(\\lambda+a_1)+(a_1-a_2)(\\lambda-a_1)e^{-2a_1d}}
e^{-\\lambda H}.
```

For ``N`` layers the implemented recursion uses

```math
k_{m,N}=\\frac{1-\\Gamma_{m,N}e^{-2d_ma_m}}
{1+\\Gamma_{m,N}e^{-2d_ma_m}},\\qquad
\\Gamma_{m,N}=\\frac{\\eta_m-\\eta_{m+1}k_{m+1,N}}
{\\eta_m+\\eta_{m+1}k_{m+1,N}},\\qquad k_{N,N}=1.
```

**Reference.** E. D. Sunde, *Earth Conduction Effects in Transmission
Systems*, Dover, 1968.
"""
function description(::Formula{:Sunde1968})
    "Sunde homogeneous and horizontally stratified overhead impedance (1968)"
end

function propagation_constant(::Val{:Sunde1968}, jω, permeability, permittivity)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Sunde1968})(
        rho, epsilon, mu, jω, Γ, segments = nothing
)
    length(rho) == 2 || throw(DimensionMismatch(
        ":Sunde1968 requires layer thicknesses for a stratified earth"
    ))
    return _homogeneous_functor(
        Val(:Sunde1968), formula, rho, epsilon, mu, jω, Γ, segments
    )
end

function (formula::Formula{:Sunde1968})(
        rho, epsilon, mu, jω, Γ, segments, thickness
)
    length(rho) >= 2 || throw(DimensionMismatch(
        ":Sunde1968 requires air and at least one earth layer"
    ))
    return _stratified_functor(
        Val(:Sunde1968), formula,
        rho, epsilon, mu, jω, Γ, segments, thickness
    )
end

"""
$(TYPEDSIGNATURES)

Select Sunde's earth-return kernel from the number of physical earth layers.

One earth layer uses the homogeneous formula, two use the explicit two-layer
formula, and larger models use Sunde's recursive transmission-line form. The
three leaves are one literature recipe and remain individually replaceable.
"""
function earth_impedance(
        ::Val{:Sunde1968}, ::Val{:mutual}, functor, pair
)
    _require(pair, Val(:overhead))
    count = length(functor.state.rho) - 1
    count == 1 && return functor.routes.homogeneous(functor, pair)
    count == 2 && return functor.routes.two_layer(functor, pair)
    return functor.routes.multilayer(functor, pair)
end

raw"""
Evaluate Sunde's homogeneous-earth overhead impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[\ln\frac{D_{ij}}{d_{ij}}+
2\int_0^\infty\frac{e^{-(h_i+h_j)\lambda}\cos(y_{ij}\lambda)}
{\lambda+\sqrt{\lambda^2+\gamma_1^2}}\,d\lambda\right],
```

with ``\gamma_1^2=j\omega\mu_0(\sigma_1+j\omega\varepsilon_1)``.
"""
function earth_impedance(
        ::Val{:Sunde1968}, ::Val{:homogeneous}, functor, pair
)
    state = functor.state
    geometry = _geometry(pair)
    gamma_squared = state.gamma_medium_squared[2]
    integral = _quadrature(state) do lambda
        attenuation = sqrt(lambda^2 + gamma_squared)
        exp(-geometry.H * lambda) * cos(geometry.y_ij * lambda) /
        (lambda + attenuation)
    end
    return state.jω * state.mu[1] / (2π) *
           (log(geometry.D_ij / geometry.d_ij) + 2 * integral)
end

raw"""
Evaluate Sunde's two-layer overhead earth-return impedance:

```math
Z_{e,ij}=\frac{j\omega\mu_0}{2\pi}\left[\ln\frac{D_{ij}}{d_{ij}}+
2\int_0^\infty F_{ij}^{S}(\lambda)\cos(y_{ij}\lambda)d\lambda\right],
```

```math
F_{ij}^{S}=\frac{a_1+a_2+(a_1-a_2)e^{-2a_1d}}
{(a_1+a_2)(\lambda+a_1)+(a_1-a_2)(\lambda-a_1)e^{-2a_1d}}
e^{-\lambda(h_i+h_j)},\qquad
a_m=\sqrt{\lambda^2+\gamma_m^2}.
```
"""
function earth_impedance(
        ::Val{:Sunde1968}, ::Val{:two_layer}, functor, pair
)
    state = functor.state
    geometry = _geometry(pair)
    d = state.thickness[2]
    integral = _quadrature(state) do lambda
        a_1 = sqrt(lambda^2 + state.gamma_medium_squared[2])
        a_2 = sqrt(lambda^2 + state.gamma_medium_squared[3])
        decay = exp(-2a_1 * d)
        numerator = a_1 + a_2 + (a_1 - a_2) * decay
        denominator = (a_1 + a_2) * (lambda + a_1) +
                      (a_1 - a_2) * (lambda - a_1) * decay
        numerator / denominator * exp(-lambda * geometry.H) *
        cos(lambda * geometry.y_ij)
    end
    return state.jω * state.mu[1] / (2π) *
           (log(geometry.D_ij / geometry.d_ij) + 2 * integral)
end

raw"""
Evaluate Sunde's recursive N-layer overhead formula as printed in the corpus:

```math
F_{ij}^{S,N}=\frac{k_{1,N}(\lambda)}{\sqrt{j\omega\mu_0}}
e^{-\lambda(h_i+h_j)},
```

```math
k_{m,N}=\frac{1-\Gamma_{m,N}e^{-2d_ma_m}}
{1+\Gamma_{m,N}e^{-2d_ma_m}},\qquad
\Gamma_{m,N}=\frac{\eta_m-\eta_{m+1}k_{m+1,N}}
{\eta_m+\eta_{m+1}k_{m+1,N}},
```

with ``k_{N,N}=1``, ``a_m=\sqrt{\lambda^2+\gamma_m^2}``, and
``\eta_m=\sqrt{j\omega\mu_0/a_m}``. The corpus rates this transcription low
confidence; this support leaf intentionally does not repair its dimensions.
"""
function earth_impedance(
        ::Val{:Sunde1968}, ::Val{:multilayer}, functor, pair
)
    state = functor.state
    geometry = _geometry(pair)
    N = length(state.rho) - 1
    T = typeof(state.jω)
    a = Vector{T}(undef, N)
    eta = similar(a)
    k = similar(a)
    integral = _quadrature(state) do lambda
        @inbounds for m in 1:N
            a[m] = sqrt(lambda^2 + state.gamma_medium_squared[m + 1])
            eta[m] = sqrt(state.jω * state.mu[1] / a[m])
        end
        k[N] = one(T)
        @inbounds for m in (N - 1):-1:1
            reflection = (eta[m] - eta[m + 1] * k[m + 1]) /
                         (eta[m] + eta[m + 1] * k[m + 1])
            decay = exp(-2 * state.thickness[m + 1] * a[m])
            k[m] = (1 - reflection * decay) / (1 + reflection * decay)
        end
        k[1] / sqrt(state.jω * state.mu[1]) *
        exp(-lambda * geometry.H) * cos(lambda * geometry.y_ij)
    end
    return state.jω * state.mu[1] / (2π) *
           (log(geometry.D_ij / geometry.d_ij) + 2 * integral)
end

:Sunde1968
