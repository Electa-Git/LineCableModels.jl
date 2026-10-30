"""
$(TYPEDSIGNATURES)

Return the declared options of the earth-return formula `Formula{:martins2026}`.

# Arguments

- `Formula{:martins2026}`: earth-return formula type.

# Returns

- `integral`: the default integration selection, `integrate => :quad`.
"""
function defaults(::Type{Formula{:martins2026}})
    return (integral = integrate => :quad,)
end

"""
$(TYPEDSIGNATURES)

Check the integration controls of the earth-return formula `Formula{:martins2026}`.

# Arguments

- `options`: formula options with an integration section.

# Returns

- `options`, unchanged.

# Errors

- Throws `ArgumentError` for a negative or nonfinite tolerance, two zero tolerances, or a
  nonpositive evaluation limit.
"""
function validate(options::FormulationOptions{Formula{:martins2026}})
    validate(options.integral, integrate, options.integral.method)
    return options
end

"""
$(TYPEDSIGNATURES)

Check that each outer cable circumference lies strictly within one medium. The
circumference-current formula requires an outer radius smaller than the distance from the
cable center to the interface.

# Arguments

- `problem`: cable geometry and earth model.
- `formula`: earth-return formula `Formula{:martins2026}` that consumes the problem.

# Returns

- `problem`, unchanged.

# Errors

- Throws `DomainError` for a cable whose outer circumference touches or crosses the
  interface.
"""
function validate(problem::LineParametersProblem, formula::Formula{:martins2026})
    system = problem.system
    for (k, design) in enumerate(system.designs)
        y, R = system.placements[k].y, support(design.outer)
        R < abs(y) || throw(DomainError((y, R),
            "the formula :martins2026 needs the cable $k of $(system.id) in one medium: its " *
            "outer radius $R m reaches the interface $(abs(y)) m from its center"))
    end
    return problem
end

# The formula gives both earth-return quantities.
validate(formula::Formula{:martins2026}, ::typeof(impedance!)) = formula
validate(formula::Formula{:martins2026}, ::typeof(potential!)) = formula

"""
$(TYPEDSIGNATURES)

Allocate cable matrices and state vectors under `martins2026`, then add the integral arrays for
the selected quadrature method.

# Arguments

- `formula`: earth-return formula `Formula{:martins2026}` with integration controls.
- `LineParametersProblem`: problem whose computation needs the arrays.
- `input`: numerical cable geometry.
- `buffers`: computation arrays already allocated.

# Returns

- `buffers` extended with the `martins2026` and `integral` blocks. [`to_buffer!`](@ref) fills the
  `martins2026` arrays.
"""
function initialize_buffers(formula::Formula{:martins2026}, ::LineParametersProblem,
        input::CoaxialInput{T}, buffers) where {T}
    c = length(input.radius)
    block = (K̃ = Matrix{Complex{T}}(undef, c, c),
        P̃ = Matrix{Complex{T}}(undef, c, c),
        T = Matrix{Complex{T}}(undef, c, c),
        Zg = Matrix{Complex{T}}(undef, c, c),
        Pg = Matrix{Complex{T}}(undef, c, c),
        factor = Matrix{Complex{T}}(undef, c, c),
        rhs = Matrix{Complex{T}}(undef, c, c),
        ξ = Vector{Complex{T}}(undef, c),
        W = Vector{Complex{T}}(undef, c),
        F = Vector{Complex{T}}(undef, c),
        i0minus = Vector{Complex{T}}(undef, c),
        scaling = Vector{T}(undef, c))
    return initialize_buffers(SpectralIntegral, formula.options.integral.method, T,
        merge(buffers, (martins2026 = block,)))
end

"""
$(TYPEDSIGNATURES)

Compute the earth-return impedances and potential coefficients of the cables at the
frequency of index `k` by the formula `Formula{:martins2026}`, whose equations are in the
section "Formulas" of [`EarthReturn`](@ref). The longitudinal wavenumber is
`buffers.gamma[k]`. The results and the shared state of the two quantities stay in
`buffers.martins2026`, which the pair variants read.

# Arguments

- `f`: earth-return formula `Formula{:martins2026}` with integration controls.
- `buffers`: computation arrays with `gamma`, `integral`, and `martins2026`.
- `input`: cable centers and outer radii.
- `media`: air and earth records with resistivity `rho` in \\[Ω·m\\] and relative `eps_r`,
  `mu_r` \\[dimensionless\\].
- `s`: complex frequency, `j 2π f`, \\[1/s\\].
- `k`: frequency index.

# Returns

- `buffers`, with `martins2026.Zg` \\[Ω/m\\], `martins2026.Pg` \\[m/F\\], and their shared state
  filled.

# Errors

- Throws `DomainError` for a nonfinite weighted earth kernel, a nonfinite integral, or a
  spectral path whose geometric weight grows.
"""
function to_buffer!(f::Formula{:martins2026}, buffers, input::CoaxialInput, media, s, k)
    block = from_buffer(f, buffers)
    section = f.options.integral
    R = eltype(buffers.integral.points)
    Γ = buffers.gamma[k]
    # The medium state of the air and of the earth.
    σ̃ = map(m -> inv(m.rho) + s * (vacuum_permittivity * m.eps_r), media)
    μ = map(m -> vacuum_permeability * m.mu_r, media)
    κ² = map((m, a) -> sqrt(s * m * a)^2 - Γ^2, μ, σ̃)
    κ = map(sqrt, κ²)
    # The candidates of the squared electric denominator are poles where the denominator is
    # nonzero and the unsquared one vanishes.
    σ̃n = map(nominal, σ̃)
    κ²n = map(nominal, κ²)
    denominator = σ̃n[2]^2 - σ̃n[1]^2
    candidates = (1, -1) .* sqrt((σ̃n[1]^2 * κ²n[2] - σ̃n[2]^2 * κ²n[1]) / denominator)
    poles = map(
        candidate -> begin
            a = σ̃n[2] * sqrt(candidate^2 + κ²n[1])
            b = σ̃n[1] * sqrt(candidate^2 + κ²n[2])
            return !iszero(denominator) && abs(a + b) <= sqrt(eps(R)) * (abs(a) + abs(b)) ?
                   candidate : zero(candidate)
        end,
        candidates)
    scales = (abs(nominal(κ[2] * σ̃[1] / σ̃[2])),
        abs(nominal(κ[1] * σ̃[2] / σ̃[1])))
    # The state of each cable.
    for i in eachindex(input.radius)
        t = input.y[i] > 0 ? 1 : 2
        ξ = κ[t] * input.radius[i]
        block.ξ[i] = ξ
        block.scaling[i] = abs(real(ξ))
        block.W[i] = besselix(0, ξ)
        block.F[i] = 2 * (one(s) * π) * σ̃[t] * input.radius[i]^2 *
                     (if abs(nominal(ξ)) < 1e-3
                         begin
                             x2 = ξ * ξ
                             one(ξ) / 2 - x2 / 16 + x2 * x2 / 96 - 11x2^3 / 6144
                         end
                     else
                         besselix(1, ξ) / (ξ * besselix(0, ξ))
                     end)
        block.i0minus[i] = i0m1(ξ)
    end
    medium = (; κ², μ, σ̃)
    πT = one(s) * π
    for j in eachindex(input.radius), i in eachindex(input.radius)

        t = input.y[i] > 0 ? 1 : 2
        u = input.y[j] > 0 ? 1 : 2
        h_i = abs(input.y[i])
        h_j = abs(input.y[j])
        r = input.radius[i]
        y = abs(input.x[i] - input.x[j])
        zero_r = zero(r)
        scaling_i = block.scaling[i]
        scaling_j = block.scaling[j]
        pair = (; t, s = u, h_i, h_j, y, r, logscale = scaling_i + scaling_j,
            scaling_j, i0minus_i = block.i0minus[i])
        # The direct and image term of two cables in one medium, averaged over the
        # circumference of the receiver and scaled as the source column.
        Λ_ij = if t == u
            D′ = hypot(input.x[i] - input.x[j], h_i + h_j)
            exp(scaling_j) * (i == j ? Λ(f, Val(:self), κ[t], r, h_i) :
             Λ(f, Val(:mutual), κ[t], r,
                hypot(input.x[i] - input.x[j], input.y[i] - input.y[j]), D′))
        else
            zero(s)
        end
        # The magnetic, electric and voltage terms share the heights and the path.
        height = h_i + h_j
        sL̃ = s / πT * block.W[i] *
              first(integrate(
                  SpectralIntegral(
                      λ -> f(Val(:J), λ, medium, pair),
                      buffers,
                      κ,
                      poles,
                      scales,
                      height,
                      y + zero_r),
                  section.method,
                  section,
                  buffers))
        t == u && (sL̃ += s * μ[t] / (2πT) * Λ_ij)
        Φ̃ = zero(sL̃)
        if !iszero(Γ)
            Φ̃ = s / πT * block.W[i] *
                 first(integrate(
                     SpectralIntegral(
                         λ -> f(Val(:M), λ, medium, pair),
                         buffers,
                         κ,
                         poles,
                         scales,
                         height,
                         y + zero_r),
                     section.method,
                     section,
                     buffers))
            t == u && (Φ̃ += s / (2πT * σ̃[t]) * Λ_ij)
        end
        K̃_ij = sL̃ - Γ^2 / s * Φ̃
        combined = t == 1 && abs(real(nominal(block.ξ[i]))) < 300
        P̃_ij = if combined
            # The potential of the air receiver to the interface, its two terms combined.
            s / πT * first(integrate(
                SpectralIntegral(
                    λ -> f(Val(:Θ), λ, medium, pair),
                    buffers,
                    κ,
                    poles,
                    scales,
                    h_j - h_j / 2,
                    y + r),
                section.method,
                section,
                buffers))
        else
            s / πT * block.W[i] *
            first(integrate(
                SpectralIntegral(
                    λ -> f(Val(:P̃), λ, medium, pair),
                    buffers,
                    κ,
                    poles,
                    scales,
                    height,
                    y + zero_r),
                section.method,
                section,
                buffers))
        end
        t == u && (P̃_ij += s / (2πT * σ̃[t]) * Λ_ij)
        if t == 1 && !combined
            # The potential of the interface below the air receiver.
            P̃_ij += s / πT * first(integrate(
                SpectralIntegral(
                    λ -> f(Val(:Θ), Val(:interface), λ, medium, pair),
                    buffers,
                    κ,
                    poles,
                    scales,
                    h_j,
                    y + r),
                section.method,
                section,
                buffers))
        end
        block.K̃[i, j] = K̃_ij
        block.P̃[i, j] = P̃_ij
        block.T[i, j] = (i == j ? inv(block.W[i]) : zero(K̃_ij)) - block.F[i] * K̃_ij
    end
    # The currents close on the circumferences.
    copyto!(block.factor, transpose(block.T))
    factor = lu!(block.factor)
    copyto!(block.rhs, transpose(block.P̃))
    ldiv!(factor, block.rhs)
    copyto!(block.Pg, transpose(block.rhs))
    @. block.Zg = block.K̃ + Γ^2 / s * block.P̃
    copyto!(block.rhs, transpose(block.Zg))
    ldiv!(factor, block.rhs)
    copyto!(block.Zg, transpose(block.rhs))
    return buffers
end

function (::Formula{:martins2026})(::Val{:J}, λ, medium, pair)
    (; κ², μ) = medium
    (; t, h_i, h_j, y, logscale) = pair
    u = pair.s
    height = h_i + h_j
    a = (sqrt(λ^2 + κ²[1]), sqrt(λ^2 + κ²[2]))
    exponent = logscale -
               h_i * (iszero(κ²[t]) ? zero(a[t]) : κ²[t] / (a[t] + λ)) -
               h_j * (iszero(κ²[u]) ? zero(a[u]) : κ²[u] / (a[u] + λ))
    amplitude = μ[1] * μ[2] / (μ[2] * a[1] + μ[1] * a[2])
    value = amplitude * (exp(exponent - (height + im * y) * λ) +
             exp(exponent - (height - im * y) * λ)) / 2
    isfinite(value) ||
        throw(DomainError((λ, value), "nonfinite weighted earth kernel"))
    return value
end

function (::Formula{:martins2026})(::Val{:M}, λ, medium, pair)
    (; κ², μ, σ̃) = medium
    (; t, h_i, h_j, y, logscale) = pair
    u = pair.s
    height = h_i + h_j
    a = (sqrt(λ^2 + κ²[1]), sqrt(λ^2 + κ²[2]))
    exponent = logscale -
               h_i *
               (iszero(κ²[t]) ? zero(a[t]) : κ²[t] / (a[t] + λ)) -
               h_j *
               (iszero(κ²[u]) ? zero(a[u]) : κ²[u] / (a[u] + λ))
    amplitude = (μ[1] * a[1] + μ[2] * a[2]) /
                ((μ[2] * a[1] + μ[1] * a[2]) *
                 (σ̃[2] * a[1] + σ̃[1] * a[2]))
    value = amplitude * (exp(exponent - (height + im * y) * λ) +
             exp(exponent - (height - im * y) * λ)) / 2
    isfinite(value) ||
        throw(DomainError((λ, value), "nonfinite weighted earth kernel"))
    return value
end

function (::Formula{:martins2026})(::Val{:Θ}, λ, medium, pair)
    (; κ², σ̃) = medium
    (; h_i, h_j, y, r, scaling_j, i0minus_i) = pair
    u = pair.s
    padding = h_j / 2
    a0 = sqrt(λ^2 + κ²[1])
    ag = sqrt(λ^2 + κ²[2])
    aq = u == 1 ? a0 : ag
    exponent = scaling_j - h_j * (iszero(κ²[u]) ? zero(aq) : κ²[u] / (aq + λ))
    amplitude = ag * (if iszero(a0)
        argument = r * λ
        j0 = exp(-padding * λ + abs(imag(argument))) *
             besseljx(0, argument)
        -h_i * j0
    elseif abs(nominal(r * λ)) < 0.5
        jminus = i0m1(im * (r * λ))
        envelope = exp(-padding * λ)
        envelope *
        (i0minus_i * exp(-h_i * a0) + expm1(-h_i * a0) - jminus) /
        a0
    else
        argument = r * λ
        envelope = exp(-padding * λ + abs(imag(argument)))
        scaled = iszero(abs(nominal(envelope))) ? zero(envelope) :
                 envelope * besseljx(0, argument)
        ((one(i0minus_i) + i0minus_i) * exp(-padding * λ - h_i * a0) -
         scaled) /
        a0
    end) / (σ̃[2] * a0 + σ̃[1] * ag)
    value = amplitude *
            (exp(exponent - ((h_j - padding) + im * y) * λ) +
             exp(exponent - ((h_j - padding) - im * y) * λ)) / 2
    isfinite(value) ||
        throw(DomainError((λ, value), "nonfinite weighted earth kernel"))
    return value
end

function (::Formula{:martins2026})(::Val{:P̃}, λ, medium, pair)
    (; κ², σ̃) = medium
    (; t, h_i, h_j, y, logscale) = pair
    u = pair.s
    height = h_i + h_j
    a = (sqrt(λ^2 + κ²[1]), sqrt(λ^2 + κ²[2]))
    exponent = logscale -
               h_i * (iszero(κ²[t]) ? zero(a[t]) : κ²[t] / (a[t] + λ)) -
               h_j * (iszero(κ²[u]) ? zero(a[u]) : κ²[u] / (a[u] + λ))
    amplitude = (t == 2 ? a[1] / a[2] : a[2] / a[1]) /
                (σ̃[2] * a[1] + σ̃[1] * a[2])
    value = amplitude * (exp(exponent - (height + im * y) * λ) +
             exp(exponent - (height - im * y) * λ)) / 2
    isfinite(value) ||
        throw(DomainError((λ, value), "nonfinite weighted earth kernel"))
    return value
end

function (::Formula{:martins2026})(::Val{:Θ}, ::Val{:interface}, λ, medium, pair)
    (; κ², σ̃) = medium
    (; h_j, y, r, scaling_j) = pair
    u = pair.s
    a = (sqrt(λ^2 + κ²[1]), sqrt(λ^2 + κ²[2]))
    exponent = scaling_j -
               h_j * (iszero(κ²[u]) ? zero(a[u]) : κ²[u] / (a[u] + λ))
    amplitude = -a[2] / (a[1] * (σ̃[2] * a[1] + σ̃[1] * a[2]))
    growth = abs(imag(r * λ))
    value = amplitude * (exp(exponent - (h_j + im * y) * λ + growth) +
             exp(exponent - (h_j - im * y) * λ + growth)) / 2
    isfinite(value) ||
        throw(DomainError((λ, value), "nonfinite weighted earth kernel"))
    return (iszero(abs(nominal(value))) ? value :
            value * besseljx(0, r * λ))
end

function (f::Formula{:martins2026})(::Val{:Z}, ::Val{:self}, i, buffers, input, media, s, k)
    from_buffer(f, buffers).Zg[i, i]
end

function (f::Formula{:martins2026})(
        ::Val{:Z}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    from_buffer(f, buffers).Zg[i, j]
end

function (f::Formula{:martins2026})(::Val{:P}, ::Val{:self}, i, buffers, input, media, s, k)
    from_buffer(f, buffers).Pg[i, i]
end

function (f::Formula{:martins2026})(
        ::Val{:P}, ::Val{:mutual}, i, j, buffers, input, media, s, k)
    from_buffer(f, buffers).Pg[i, j]
end

"""
$(TYPEDSIGNATURES)

Evaluate the modified Bessel difference without subtraction near zero:

```math
I_0(z)-1=\\sum_{n=1}^{\\infty}\\frac{(z^2/4)^n}{(n!)^2}.
```

Here ``I_0`` is the modified Bessel function of order zero. For nominal ``|z|\\leq0.5``,
sum at most 1000 terms and stop when another term leaves the sum unchanged. Otherwise,
evaluate ``e^{|\\operatorname{Re}z|}I_{0,\\mathrm{scaled}}(z)-1`` with the exponentially
scaled Bessel function.

# Arguments

- `z`: Bessel argument in \\[dimensionless\\].

# Returns

- The difference ``I_0(z)-1``, in \\[dimensionless\\].
"""
function i0m1(z)
    if abs(nominal(z)) > 0.5
        return besselix(0, z) * exp(abs(real(z))) - one(z)
    end
    term = z * z / 4
    result = term
    for n in 2:1000
        term *= (z / (2n))^2
        next = result + term
        next == result && return next
        result = next
    end
    return result
end
