function routes(::Val{:Nakagawa1973})
    (
        self = nakagawa_recursive1973,
        mutual = nakagawa_recursive1973,
        Γ = nakagawa_recursive1973_gamma
    )
end

function assumptions(::Val{:Nakagawa1973})
    (
        air = _full,
        earth = _full,
        permeability = _material
    )
end

propagation(::Val{:Nakagawa1973}) = Val(:zero)
media(::Formula{:Nakagawa1973}) = Val(:stratified)
"""
$(TYPEDSIGNATURES)

**Identification.** Recursive ``N``-layer overhead-earth impedance.

**Expression.**

```math
F_{ij}^{N}=\\frac{A_1+B_1}
{(\\lambda+\\mu_0a_1/\\mu_1)A_1+(\\lambda-\\mu_0a_1/\\mu_1)B_1}
e^{-\\lambda H},
```

```math
A_{N-1}=b_{N-1}+b_N,\\qquad
B_{N-1}=(b_{N-1}-b_N)e^{-2a_{N-1}t_{N-1}},\\qquad
b_m=\\frac{a_m}{\\mu_m}.
```

The implementation evaluates the published upward ``A_m,B_m`` recursion and
inserts ``F_{ij}^{N}`` in the standard overhead integral.

**Reference.** M. Nakagawa, A. Ametani, and K. Iwamoto, “Further Studies on
Wave Propagation in Overhead Lines with Earth Return: Impedance of Stratified
Earth,” *Proceedings of the IEE*, 120, 1521–1528, 1973.
"""
function description(::Formula{:Nakagawa1973})
    "Nakagawa et al. recursive N-layer overhead impedance (1973)"
end

function nakagawa_recursive1973_gamma(jω, permeability, permittivity)
    (Γ = zero(jω), squared = zero(jω))
end

function (formula::Formula{:Nakagawa1973})(
        rho, epsilon, mu, jω, Γ, segments, thickness
)
    length(rho) >= 3 || throw(DimensionMismatch(
        ":Nakagawa1973 requires air and at least two earth layers"
    ))
    return _stratified_functor(
        Val(:Nakagawa1973), formula,
        rho, epsilon, mu, jω, Γ, segments, thickness
    )
end

raw"""
Evaluate the Nakagawa et al. recursive N-layer overhead impedance:

```math
F_{ij}^{N}=\frac{A_1+B_1}
{(\lambda+\mu_0a_1/\mu_1)A_1+
(\lambda-\mu_0a_1/\mu_1)B_1}e^{-\lambda(h_i+h_j)}.
```

With ``b_m=a_m/\mu_m`` and
``a_m=\sqrt{\lambda^2+\gamma_m^2-\gamma_0^2}``, the bottom termination is

```math
A_{N-1}=b_{N-1}+b_N,\qquad
B_{N-1}=(b_{N-1}-b_N)e^{-2a_{N-1}t_{N-1}},
```

and the exact upward ``A_m,B_m`` recursions from the corpus are retained.
The three-layer record `NakagawaEtAl1973` is not duplicated because it is the
``N=3`` specialization of this route.
"""
function nakagawa_recursive1973(functor, pair)
    _require(pair, Val(:overhead))
    state = functor.state
    geometry = _geometry(pair)
    N = length(state.rho) - 1
    T = typeof(state.jω)
    a = Vector{T}(undef, N)
    b = similar(a)
    A = similar(a)
    B = similar(a)
    depths = cumsum(@view state.thickness[2:(end - 1)])
    gamma_0_squared = state.gamma_medium_squared[1]
    integral = _quadrature(state) do lambda
        @inbounds for m in 1:N
            layer = m + 1
            a[m] = sqrt(lambda^2 + state.gamma_medium_squared[layer] - gamma_0_squared)
            b[m] = a[m] / state.mu[layer]
        end
        @inbounds begin
            A[N - 1] = b[N - 1] + b[N]
            B[N - 1] = (b[N - 1] - b[N]) *
                       exp(-2a[N - 1] * depths[N - 1])
            for m in (N - 2):-1:1
                bridge = exp(-2a[m + 1] * depths[m])
                A[m] = (b[m] + b[m + 1]) * A[m + 1] +
                       (b[m] - b[m + 1]) * B[m + 1] * bridge
                B[m] = ((b[m] - b[m + 1]) * A[m + 1] +
                        (b[m] + b[m + 1]) * B[m + 1] * bridge) *
                       exp(-2a[m] * depths[m])
            end
        end
        numerator = A[1] + B[1]
        denominator = (lambda + state.mu[1] * b[1]) * A[1] +
                      (lambda - state.mu[1] * b[1]) * B[1]
        numerator / denominator * exp(-lambda * geometry.H) *
        cos(lambda * geometry.y_ij)
    end
    return state.jω * state.mu[1] / (2π) *
           (log(geometry.D_ij / geometry.d_ij) + 2 * integral)
end

:Nakagawa1973
