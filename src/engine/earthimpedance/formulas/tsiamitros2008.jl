function routes(identifier::Val{:Tsiamitros2008})
    (
        self = formula_method(identifier, earth_impedance, Val(:self)),
        mutual = formula_method(identifier, earth_impedance, Val(:mutual)),
        Γ = formula_method(identifier, propagation_constant)
    )
end

function assumptions(::Val{:Tsiamitros2008})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Tsiamitros2008}) = Val(:zero)
media(::Formula{:Tsiamitros2008}) = Val(:stratified)
"""
$(TYPEDSIGNATURES)

**Identification.** General same-layer, cross-layer, and mixed conductor
kernel for arbitrarily stratified earth.

**Expression.** For source layer ``m`` and target layer ``l``,

```math
Z_{e,ij}=\\frac{j\\omega\\mu_m}{2\\pi}\\int_0^\\infty
\\frac{\\cos(uy_{ij})}{\\bar\\alpha_m}
\\left[2^{m-l}
\\frac{(\\mu_1\\cdots\\mu_{m-1})(\\bar\\alpha_1\\cdots\\bar\\alpha_m)}
{(\\mu_1\\cdots\\mu_{l-1})(\\bar\\alpha_1\\cdots\\bar\\alpha_l)}
\\frac{e^{-\\sum_{q=l}^{m}\\bar\\alpha_qd_q}\\bar F_1\\bar F_2}
{\\overline{DTD}_0}\\right]du.
```

The ``\\overline{DTD}``, ``\\overline{DTN}``, ``\\overline{TDD}``, and
``\\overline{TDN}`` interface factors are evaluated by the paper's upward and
downward recursions.

**Reference.** D. A. Tsiamitros, G. K. Papagiannis, and P. S. Dokopoulos,
“Earth Return Impedances of Conductor Arrangements in Multilayer Soils—Part
I: Theoretical Model,” *IEEE Transactions on Power Delivery*, 23(4),
2392–2400, 2008.
"""
function description(::Formula{:Tsiamitros2008})
    "Tsiamitros et al. true arbitrary-layer impedance kernel (2008)"
end

function propagation_constant(
        ::Val{:Tsiamitros2008}, jω, permeability, permittivity
)
    return (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Tsiamitros2008})(
        rho, epsilon, mu, jω, Γ, segments, thickness
)
    length(rho) >= 2 || throw(DimensionMismatch(
        ":Tsiamitros2008 requires air and at least one earth layer"
    ))
    return _stratified_functor(
        Val(:Tsiamitros2008), formula,
        rho, epsilon, mu, jω, Γ, segments, thickness
    )
end

raw"""
Evaluate the Tsiamitros et al. true multilayer impedance kernel:

```math
Z_{e,ij}=\frac{j\omega\mu_m}{2\pi}\int_0^\infty
\frac{\cos(u y_{ij})}{\bar\alpha_m}
\left\{2^{m-l}
\frac{(\mu_1\cdots\mu_{m-1})
(\bar\alpha_1\cdots\bar\alpha_m)
(e^{-\bar\alpha_l d_l}\cdots e^{-\bar\alpha_m d_m})}
{(\mu_1\cdots\mu_{l-1})
(\bar\alpha_1\cdots\bar\alpha_l)}
\frac{\bar F_1\bar F_2}{\overline{DTD}_0}\right\}du,
```

```math
\bar F_1=\overline{DTD}_m e^{\bar\alpha_m(d_m-h_1)}+
\overline{DTN}_m e^{-\bar\alpha_m(d_m-h_1)},\qquad
\bar F_2=\overline{TDD}_{l-1}e^{\bar\alpha_lh_2}+
\overline{TDN}_{l-1}e^{-\bar\alpha_lh_2}.
```

The downward recursion terminates with
``\overline{DTN}_N=0,\overline{DTD}_N=1``; the upward recursion starts from
``\overline{TDD}_{-1}=1,\overline{TDN}_{-1}=0``. The air-to-earth route is
the paper's Eq. 14 specialization. The implementation combines exponential
factors before evaluation, so conductors may lie in the infinite bottom
layer without forming ``0e^{-\infty}``.

```math
\overline{DTN}_k=(\mu_{k+1}\bar\alpha_k-\mu_k\bar\alpha_{k+1})
\overline{DTD}_{k+1}+(\mu_{k+1}\bar\alpha_k+\mu_k\bar\alpha_{k+1})
\overline{DTN}_{k+1}e^{-2\bar\alpha_{k+1}d_{k+1}},
```

```math
\overline{DTD}_k=(\mu_{k+1}\bar\alpha_k+\mu_k\bar\alpha_{k+1})
\overline{DTD}_{k+1}+(\mu_{k+1}\bar\alpha_k-\mu_k\bar\alpha_{k+1})
\overline{DTN}_{k+1}e^{-2\bar\alpha_{k+1}d_{k+1}}.
```

The corresponding top-down recursion is

```math
\overline{TDD}_{k-1}=(\mu_{k-1}\bar\alpha_k+\mu_k\bar\alpha_{k-1})
\overline{TDD}_{k-2}+(\mu_{k-1}\bar\alpha_k-\mu_k\bar\alpha_{k-1})
\overline{TDN}_{k-2}e^{-2\bar\alpha_{k-1}d_{k-1}},
```

```math
\overline{TDN}_{k-1}=(\mu_{k-1}\bar\alpha_k-\mu_k\bar\alpha_{k-1})
\overline{TDD}_{k-2}+(\mu_{k-1}\bar\alpha_k+\mu_k\bar\alpha_{k-1})
\overline{TDN}_{k-2}e^{-2\bar\alpha_{k-1}d_{k-1}}.
```

As specified by the author, self impedance is the same-layer mutual route
with both depths equal and ``y_{ii}`` set to the conductor outer radius. The
canonical `EarthPair` already carries that radius as its self separation.
"""
function earth_impedance(
        identifier::Val{:Tsiamitros2008}, ::Val{:mutual}, functor, pair
)
    return earth_impedance(identifier, _placement(pair), functor, pair)
end

function earth_impedance(
        identifier::Val{:Tsiamitros2008}, ::Val{:overhead}, functor, pair
)
    state = functor.state
    geometry = _geometry(pair)
    count = length(state.rho)
    T = typeof(state.jω)
    attenuation = Vector{T}(undef, count)
    numerator = similar(attenuation)
    denominator = similar(attenuation)
    difference = abs(geometry.h_i - geometry.h_j)
    integral = _quadrature(state) do lambda
        downward_coefficients!(
            identifier, numerator, denominator, attenuation, state, lambda
        )
        a_0 = attenuation[1]
        reflection = numerator[1] / denominator[1]
        exact = (
            exp(-a_0 * difference) +
            reflection * exp(-a_0 * geometry.H)
        ) / a_0
        span = geometry.H - difference
        reference = iszero(lambda) ? span :
                    exp(-lambda * difference) *
                    (-expm1(-lambda * span)) / lambda
        (exact - reference) * cos(geometry.y_ij * lambda)
    end
    return state.jω * state.mu[1] / (2π) *
           (log(geometry.D_ij / geometry.d_ij) + integral)
end

function earth_impedance(
        identifier::Val{:Tsiamitros2008}, ::Val{:underground}, functor, pair
)
    state = functor.state
    source, target = conductor_order(identifier, pair)
    m = pair.layers[source] - 1
    l = pair.layers[target] - 1
    h_1 = local_layer_depth(identifier, state, pair, source, m)
    h_2 = local_layer_depth(identifier, state, pair, target, l)
    count = length(state.rho)
    T = typeof(state.jω)
    attenuation = Vector{T}(undef, count)
    numerator = similar(attenuation)
    denominator = similar(attenuation)
    same_layer = m == l
    direct = same_layer ?
             special_besselk(
        0,
        state.gamma[m + 1] * hypot(pair.separation, h_1 - h_2)
    ) : zero(T)
    integral = _quadrature(state) do lambda
        downward_coefficients!(
            identifier, numerator, denominator, attenuation, state, lambda
        )
        value = spectral_kernel(
            identifier,
            state,
            attenuation,
            numerator,
            denominator,
            l,
            m,
            h_1,
            h_2
        )
        if same_layer
            value -= exp(-attenuation[m + 1] * (h_1 - h_2)) /
                     attenuation[m + 1]
        end
        value * cos(pair.separation * lambda)
    end
    return state.jω * state.mu[m + 1] / (2π) * (direct + integral)
end

function earth_impedance(
        identifier::Val{:Tsiamitros2008}, ::Val{:mixed}, functor, pair
)
    state = functor.state
    air_position = pair.layers[1] == 1 ? 1 : 2
    earth_position = air_position == 1 ? 2 : 1
    earth_layer = pair.layers[earth_position]
    m = earth_layer - 1
    Nlayer = length(state.rho) - 1
    h_air = abs(pair.heights[air_position])
    local_depth = local_layer_depth(identifier, state, pair, earth_position, m)
    T = typeof(state.jω)
    a = Vector{T}(undef, Nlayer + 1)
    numerator = similar(a)
    denominator = similar(a)
    integral = _quadrature(state) do lambda
        downward_coefficients!(
            identifier, numerator, denominator, a, state, lambda
        )
        coefficient = exp(-a[1] * h_air) * 2^(m - 1)
        attenuation_prior = one(T)
        @inbounds for k in 1:m
            coefficient *= state.mu[k + 1]
            k < m && (coefficient *= a[k + 1])
            k < m && (attenuation_prior *= exp(
                -a[k + 1] * state.thickness[k + 1]
            ))
        end
        field = denominator[m + 1] * exp(-a[m + 1] * local_depth)
        if !iszero(numerator[m + 1])
            d_m = state.thickness[m + 1]
            field += numerator[m + 1] *
                     exp(-a[m + 1] * (2d_m - local_depth))
        end
        coefficient * attenuation_prior * field / denominator[1] *
        cos(pair.separation * lambda)
    end
    return state.jω * state.mu[1] / π * integral
end

function downward_coefficients!(
        ::Val{:Tsiamitros2008},
        numerator,
        denominator,
        attenuation,
        state,
        lambda
)
    N = length(state.rho) - 1
    T = eltype(attenuation)
    @inbounds for medium in 0:N
        attenuation[medium + 1] = sqrt(
            lambda^2 + state.gamma_medium_squared[medium + 1]
        )
    end
    numerator[N + 1] = zero(T)
    denominator[N + 1] = one(T)
    @inbounds for k in (N - 1):-1:0
        left = k + 1
        right = k + 2
        contrast = state.mu[right] * attenuation[left] -
                   state.mu[left] * attenuation[right]
        sumterm = state.mu[right] * attenuation[left] +
                  state.mu[left] * attenuation[right]
        reflected = iszero(numerator[right]) ? zero(T) :
                    numerator[right] * exp(
            -2attenuation[right] * state.thickness[right]
        )
        numerator[left] = contrast * denominator[right] + sumterm * reflected
        denominator[left] = sumterm * denominator[right] + contrast * reflected
    end
    return nothing
end

function upward_coefficients(
        ::Val{:Tsiamitros2008}, state, attenuation, layer
)
    T = eltype(attenuation)
    denominator = one(T)
    numerator = zero(T)
    @inbounds for boundary in 1:layer
        upper = boundary
        lower = boundary + 1
        sumterm = state.mu[upper] * attenuation[lower] +
                  state.mu[lower] * attenuation[upper]
        contrast = state.mu[upper] * attenuation[lower] -
                   state.mu[lower] * attenuation[upper]
        reflected = boundary == 1 || iszero(numerator) ? zero(T) :
                    numerator * exp(
            -2attenuation[upper] * state.thickness[upper]
        )
        denominator, numerator = (
            sumterm * denominator + contrast * reflected,
            contrast * denominator + sumterm * reflected
        )
    end
    return denominator, numerator
end

function spectral_kernel(
        identifier::Val{:Tsiamitros2008},
        state,
        attenuation,
        down_numerator,
        down_denominator,
        l,
        m,
        h_1,
        h_2
)
    T = eltype(attenuation)
    top_denominator, top_numerator = upward_coefficients(
        identifier, state, attenuation, l
    )
    coefficient = T(2^(m - l))
    @inbounds for layer in l:(m - 1)
        coefficient *= state.mu[layer + 1]
    end
    if l == m
        coefficient /= attenuation[m + 1]
        direct = down_denominator[m + 1] * top_denominator *
                 exp(-attenuation[m + 1] * (h_1 - h_2))
        upper_image = down_denominator[m + 1] * top_numerator *
                      exp(-attenuation[m + 1] * (h_1 + h_2))
        field = direct + upper_image
        if !iszero(down_numerator[m + 1])
            d_m = state.thickness[m + 1]
            lower_image = down_numerator[m + 1] * top_denominator *
                          exp(-attenuation[m + 1] * (2d_m - h_1 - h_2))
            double_image = down_numerator[m + 1] * top_numerator *
                           exp(-attenuation[m + 1] * (2d_m - h_1 + h_2))
            field += lower_image + double_image
        end
        return coefficient * field / down_denominator[1]
    else
        @inbounds for layer in (l + 1):(m - 1)
            coefficient *= attenuation[layer + 1]
        end
        @inbounds for layer in (l + 1):(m - 1)
            coefficient *= exp(
                -attenuation[layer + 1] * state.thickness[layer + 1]
            )
        end
        d_m = state.thickness[m + 1]
        source = down_denominator[m + 1] * exp(-attenuation[m + 1] * h_1)
        if !iszero(down_numerator[m + 1])
            source += down_numerator[m + 1] *
                      exp(-attenuation[m + 1] * (2d_m - h_1))
        end
        d_l = state.thickness[l + 1]
        target = top_denominator *
                 exp(-attenuation[l + 1] * (d_l - h_2)) +
                 top_numerator * exp(-attenuation[l + 1] * (d_l + h_2))
        return coefficient * source * target / down_denominator[1]
    end
end

function conductor_order(::Val{:Tsiamitros2008}, pair)
    first_layer, second_layer = pair.layers
    first_depth, second_depth = abs.(pair.heights)
    if first_layer > second_layer ||
       (first_layer == second_layer && first_depth >= second_depth)
        return 1, 2
    end
    return 2, 1
end

function local_layer_depth(
        ::Val{:Tsiamitros2008}, state, pair, position, layer
)
    depth = abs(pair.heights[position])
    top = layer == 1 ? zero(depth) : sum(@view state.thickness[2:layer])
    local_depth = depth - top
    thickness = state.thickness[layer + 1]
    valid = local_depth >= zero(depth) &&
            (!isfinite(thickness) || local_depth <= thickness)
    valid || throw(ArgumentError(
        ":Tsiamitros2008 conductor depth is outside its resolved earth layer"
    ))
    return local_depth
end

:Tsiamitros2008
