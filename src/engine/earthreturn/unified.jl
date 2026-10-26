"""
$(TYPEDSIGNATURES)

Return the options of the formula `:unified` with their defaults. `integral` is the section
of the method that evaluates the spectral integrals.
"""
function defaults(::Type{Formula{:unified}})
    return (integral = integrate => :quad,)
end

"""
$(TYPEDSIGNATURES)

Check the controls of the integration section of the formula `:unified`.

# Errors

- Throws `ArgumentError` for controls that `integrate` refuses.
"""
function validate(options::FormulationOptions{Formula{:unified}})
    validate(options.integral, integrate, options.integral.method)
    return options
end

"""
$(TYPEDSIGNATURES)

Check that `formula` can compute `problem`. The formula averages each source over the outer
circumference of its cable, which must lie wholly in the air or wholly in the earth. The
outer radius of each cable stays below its distance to the interface.

# Errors

- Throws `DomainError` for a cable whose outer circumference touches or crosses the
  interface.
"""
function validate(problem::LineParametersProblem, formula::Formula{:unified})
    system = problem.system
    for (k, design) in enumerate(system.designs)
        y, R = system.placements[k].y, support(design.outer)
        R < abs(y) || throw(DomainError((y, R),
            "the unified formula needs the cable $k of $(system.id) in one medium: its " *
            "outer radius $R m reaches the interface $(abs(y)) m from its center"))
    end
    return problem
end

# The formula gives both earth-return quantities.
validate(formula::Formula{:unified}, ::typeof(impedance!)) = formula
validate(formula::Formula{:unified}, ::typeof(potential!)) = formula

"""
$(TYPEDSIGNATURES)

Return `buffers` with the arrays of the formula under the key `unified`, and the engine's
integral block of its integration method. The block `unified` has the matrices `K`, `H`, `L`,
`Ze`, `Pe`, `factor` and `rhs` of the cables, and their state `x`, `A`, `F`, `i0minus` and
`scaling`.
"""
function initialize_buffers(formula::Formula{:unified}, ::LineParametersProblem,
        input::CoaxialInput{T}, buffers) where {T}
    c = length(input.radius)
    unified = (K = Matrix{Complex{T}}(undef, c, c),
        H = Matrix{Complex{T}}(undef, c, c),
        L = Matrix{Complex{T}}(undef, c, c),
        Ze = Matrix{Complex{T}}(undef, c, c),
        Pe = Matrix{Complex{T}}(undef, c, c),
        factor = Matrix{Complex{T}}(undef, c, c),
        rhs = Matrix{Complex{T}}(undef, c, c),
        x = Vector{Complex{T}}(undef, c),
        A = Vector{Complex{T}}(undef, c),
        F = Vector{Complex{T}}(undef, c),
        i0minus = Vector{Complex{T}}(undef, c),
        scaling = Vector{T}(undef, c))
    return initialize_buffers(SpectralIntegral, formula.options.integral.method, T,
        merge(buffers, (unified = unified,)))
end

"""
$(TYPEDSIGNATURES)

Fill `buffers.unified` with the earth-return state of the cables at the frequency of index
`k`, of complex frequency `s` \\[1/s\\]. `media` gives the resistivity, the relative
permittivity and the relative permeability of the air and of the earth. The longitudinal
wavenumber ``\\Gamma`` \\[1/m\\] of the frequency is `buffers.gamma[k]`, which the
formulation fills from its option `Γ`.

The formula averages the field of each source over the outer circumference of its cable.
Each medium ``m`` has the admittivity ``\\sigma_m + s \\varepsilon_m`` and the wavenumber
``k_m = \\sqrt{s \\mu_m (\\sigma_m + s \\varepsilon_m) - \\Gamma^2}``. Each ordered pair of
cables gives an impedance term ``k`` and a potential term ``h``. These are spectral
integrals along the ray of `ray!`, with the direct and image term ``K_0(k D) - K_0(k d)``
of two cables in one medium. An air receiver takes its potential to the interface below it.
The currents close on the circumferences:

```math
P_e L = H, \\qquad Z_e L = K + \\Gamma^2 H / s, \\qquad L = A^{-1} - F K,
```

for the earth-return impedances `Ze` \\[Ω/m\\] and potential coefficients `Pe` \\[m/F\\], where
``A`` holds ``I_0(x_p)`` and ``F`` the current factors of the cables.
"""
function prepare!(formula::Formula{:unified}, buffers, input::CoaxialInput, media, s, k)
    unified = buffers.unified
    section = formula.options.integral
    R = eltype(buffers.integral.points)
    Γ = buffers.gamma[k]
    # The medium state of the air and of the earth.
    sh = map(m -> inv(m.rho) + s * (vacuum_permittivity * m.eps_r), media)
    mu = map(m -> vacuum_permeability * m.mu_r, media)
    k2 = map((m, a) -> sqrt(s * m * a)^2 - Γ^2, mu, sh)
    wavenumbers = map(sqrt, k2)
    # The candidates of the squared electric denominator are poles where the denominator is
    # nonzero and the unsquared one vanishes.
    shn = map(nominal, sh)
    k2n = map(nominal, k2)
    denominator = shn[2]^2 - shn[1]^2
    candidates = (1, -1) .* sqrt((shn[1]^2 * k2n[2] - shn[2]^2 * k2n[1]) / denominator)
    poles = map(
        candidate -> begin
            a = shn[2] * sqrt(candidate^2 + k2n[1])
            b = shn[1] * sqrt(candidate^2 + k2n[2])
            return !iszero(denominator) && abs(a + b) <= sqrt(eps(R)) * (abs(a) + abs(b)) ?
                   candidate : zero(candidate)
        end,
        candidates)
    scales = (abs(nominal(wavenumbers[2] * sh[1] / sh[2])),
        abs(nominal(wavenumbers[1] * sh[2] / sh[1])))
    # The state of each cable.
    for p in eachindex(input.radius)
        medium = input.y[p] > 0 ? 1 : 2
        x = wavenumbers[medium] * input.radius[p]
        unified.x[p] = x
        unified.scaling[p] = abs(real(x))
        unified.A[p] = besselix(0, x)
        unified.F[p] = 2 * (one(s) * π) * sh[medium] * input.radius[p]^2 *
                       (if abs(nominal(x)) < 1e-3
                           begin
                               x2 = x * x
                               one(x) / 2 - x2 / 16 + x2 * x2 / 96 - 11x2^3 / 6144
                           end
                       else
                           besselix(1, x) / (x * besselix(0, x))
                       end)
        unified.i0minus[p] = i0m1(x)
    end
    πT = one(s) * π
    for q in eachindex(input.radius), p in eachindex(input.radius)

        P = input.y[p] > 0 ? 1 : 2
        Q = input.y[q] > 0 ? 1 : 2
        hp = abs(input.y[p])
        hq = abs(input.y[q])
        r = input.radius[p]
        y = abs(input.x[p] - input.x[q])
        zero_r = zero(r)
        sp = unified.scaling[p]
        sq = unified.scaling[q]
        # The direct and image term of two cables in one medium, averaged over the
        # circumference of the receiver and scaled as the source column.
        own = if P == Q
            kP = wavenumbers[P]
            d = hypot(input.x[p] - input.x[q], hp + hq)
            D = p == q ? r : hypot(input.x[p] - input.x[q], input.y[p] - input.y[q])
            if abs(nominal(kP * d)) < 0.5
                difference = iszero(kP) ? log(d / D) :
                             besselk(0, kP * D) - besselk(0, kP * d)
                if p == q
                    correction = iszero(kP) ? zero(kP) :
                                 unified.i0minus[p] * besselk(0, kP * d)
                    exp(sq) * (difference - correction)
                else
                    unified.A[p] * exp(sp + sq) * difference
                end
            elseif p == q
                direct = besselkx(0, kP * r) * exp(sq - kP * r)
                image = unified.A[p] * besselkx(0, kP * d) * exp(sp + sq - kP * d)
                direct - image
            else
                unified.A[p] * (besselkx(0, kP * D) * exp(sp + sq - kP * D) -
                 besselkx(0, kP * d) * exp(sp + sq - kP * d))
            end
        else
            zero(s)
        end
        # The magnetic, electric and voltage terms share the heights and the ray.
        height = hp + hq
        logscale = sp + sq
        ray = ray!(buffers, wavenumbers, poles, scales, height, y + zero_r)
        z = s / πT * unified.A[p] *
            first(integrate(
                SpectralIntegral(t -> begin
                    λ = ray * t
                    a = (sqrt(λ^2 + k2[1]), sqrt(λ^2 + k2[2]))
                    exponent = logscale -
                               hp * (iszero(k2[P]) ? zero(a[P]) : k2[P] / (a[P] + λ)) -
                               hq * (iszero(k2[Q]) ? zero(a[Q]) : k2[Q] / (a[Q] + λ))
                    amplitude = mu[1] * mu[2] / (mu[2] * a[1] + mu[1] * a[2])
                    value = amplitude * (exp(exponent - (height + im * y) * λ) +
                             exp(exponent - (height - im * y) * λ)) / 2
                    isfinite(value) ||
                        throw(DomainError((λ, value), "nonfinite weighted earth kernel"))
                    return ray * value
                end),
                section.method,
                section,
                buffers;
                points = buffers.integral.points))
        P == Q && (z += s * mu[P] / (2πT) * own)
        phi = zero(z)
        if !iszero(Γ)
            phi = s / πT * unified.A[p] *
                  first(integrate(
                      SpectralIntegral(t -> begin
                          λ = ray * t
                          a = (sqrt(λ^2 + k2[1]), sqrt(λ^2 + k2[2]))
                          exponent = logscale -
                                     hp *
                                     (iszero(k2[P]) ? zero(a[P]) : k2[P] / (a[P] + λ)) -
                                     hq *
                                     (iszero(k2[Q]) ? zero(a[Q]) : k2[Q] / (a[Q] + λ))
                          amplitude = (mu[1] * a[1] + mu[2] * a[2]) /
                                      ((mu[2] * a[1] + mu[1] * a[2]) *
                                       (sh[2] * a[1] + sh[1] * a[2]))
                          value = amplitude * (exp(exponent - (height + im * y) * λ) +
                                   exp(exponent - (height - im * y) * λ)) / 2
                          isfinite(value) ||
                              throw(DomainError((λ, value), "nonfinite weighted earth kernel"))
                          return ray * value
                      end),
                      section.method,
                      section,
                      buffers;
                      points = buffers.integral.points))
            P == Q && (phi += s / (2πT * sh[P]) * own)
        end
        kpq = z - Γ^2 / s * phi
        combined = P == 1 && abs(real(nominal(unified.x[p]))) < 300
        h = if combined
            # The potential of the air receiver to the interface, its two terms combined.
            padding = hq / 2
            i0minus = unified.i0minus[p]
            air_ray = ray!(buffers, wavenumbers, poles, scales, hq - padding, y + r)
            s / πT * first(integrate(
                SpectralIntegral(t -> begin
                    λ = air_ray * t
                    a0 = sqrt(λ^2 + k2[1])
                    ag = sqrt(λ^2 + k2[2])
                    aq = Q == 1 ? a0 : ag
                    exponent = sq - hq * (iszero(k2[Q]) ? zero(aq) : k2[Q] / (aq + λ))
                    amplitude = ag * (if iszero(a0)
                        argument = r * λ
                        j0 = exp(-padding * λ + abs(imag(argument))) *
                             besseljx(0, argument)
                        -hp * j0
                    elseif abs(nominal(r * λ)) < 0.5
                        jminus = i0m1(im * (r * λ))
                        envelope = exp(-padding * λ)
                        envelope * (i0minus * exp(-hp * a0) + expm1(-hp * a0) - jminus) /
                        a0
                    else
                        argument = r * λ
                        envelope = exp(-padding * λ + abs(imag(argument)))
                        scaled = iszero(abs(nominal(envelope))) ? zero(envelope) :
                                 envelope * besseljx(0, argument)
                        ((one(i0minus) + i0minus) * exp(-padding * λ - hp * a0) - scaled) /
                        a0
                    end) / (sh[2] * a0 + sh[1] * ag)
                    value = amplitude * (exp(exponent - ((hq - padding) + im * y) * λ) +
                             exp(exponent - ((hq - padding) - im * y) * λ)) / 2
                    isfinite(value) ||
                        throw(DomainError((λ, value), "nonfinite weighted earth kernel"))
                    return air_ray * value
                end),
                section.method,
                section,
                buffers;
                points = buffers.integral.points))
        else
            s / πT * unified.A[p] *
            first(integrate(
                SpectralIntegral(t -> begin
                    λ = ray * t
                    a = (sqrt(λ^2 + k2[1]), sqrt(λ^2 + k2[2]))
                    exponent = logscale -
                               hp * (iszero(k2[P]) ? zero(a[P]) : k2[P] / (a[P] + λ)) -
                               hq * (iszero(k2[Q]) ? zero(a[Q]) : k2[Q] / (a[Q] + λ))
                    amplitude = (P == 2 ? a[1] / a[2] : a[2] / a[1]) /
                                (sh[2] * a[1] + sh[1] * a[2])
                    value = amplitude * (exp(exponent - (height + im * y) * λ) +
                             exp(exponent - (height - im * y) * λ)) / 2
                    isfinite(value) ||
                        throw(DomainError((λ, value), "nonfinite weighted earth kernel"))
                    return ray * value
                end),
                section.method,
                section,
                buffers;
                points = buffers.integral.points))
        end
        P == Q && (h += s / (2πT * sh[P]) * own)
        if P == 1 && !combined
            # The potential of the interface below the air receiver.
            reference_ray = ray!(buffers, wavenumbers, poles, scales, hq, y + r)
            h += s / πT * first(integrate(
                SpectralIntegral(t -> begin
                    λ = reference_ray * t
                    a = (sqrt(λ^2 + k2[1]), sqrt(λ^2 + k2[2]))
                    exponent = sq - hq * (iszero(k2[Q]) ? zero(a[Q]) : k2[Q] / (a[Q] + λ))
                    amplitude = -a[2] / (a[1] * (sh[2] * a[1] + sh[1] * a[2]))
                    growth = abs(imag(r * λ))
                    value = amplitude * (exp(exponent - (hq + im * y) * λ + growth) +
                             exp(exponent - (hq - im * y) * λ + growth)) / 2
                    isfinite(value) ||
                        throw(DomainError((λ, value), "nonfinite weighted earth kernel"))
                    return reference_ray * (iszero(abs(nominal(value))) ? value :
                            value * besseljx(0, r * λ))
                end),
                section.method,
                section,
                buffers;
                points = buffers.integral.points))
        end
        unified.K[p, q] = kpq
        unified.H[p, q] = h
        unified.L[p, q] = (p == q ? inv(unified.A[p]) : zero(kpq)) - unified.F[p] * kpq
    end
    # The currents close on the circumferences.
    copyto!(unified.factor, transpose(unified.L))
    factor = lu!(unified.factor)
    copyto!(unified.rhs, transpose(unified.H))
    ldiv!(factor, unified.rhs)
    copyto!(unified.Pe, transpose(unified.rhs))
    @. unified.Ze = unified.K + Γ^2 / s * unified.H
    copyto!(unified.rhs, transpose(unified.Ze))
    ldiv!(factor, unified.rhs)
    copyto!(unified.Ze, transpose(unified.rhs))
    return buffers
end

"""
$(TYPEDSIGNATURES)

Write into `Zg` the earth-return impedances of the cables \\[Ω/m\\] that `prepare!` left in
`buffers.unified`.
"""
function impedance!(formula::Formula{:unified}, Zg::AbstractMatrix, buffers)
    return copyto!(Zg, buffers.unified.Ze)
end

"""
$(TYPEDSIGNATURES)

Write into `Pg` the earth-return potential coefficients of the cables \\[m/F\\] that
`prepare!` left in `buffers.unified`.
"""
function potential!(formula::Formula{:unified}, Pg::AbstractMatrix, buffers)
    return copyto!(Pg, buffers.unified.Pe)
end

"""
$(TYPEDSIGNATURES)

Return ``I_0(z) - 1`` without the cancellation of the difference. For a nominal ``|z|`` of at
most 0.5, it sums the power series of ``I_0`` from its second term. Otherwise, it uses the
exponentially scaled Bessel function.
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
