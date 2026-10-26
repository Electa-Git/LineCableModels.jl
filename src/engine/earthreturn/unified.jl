"""
$(TYPEDSIGNATURES)

Return the declared options of the unified earth-return formula.

# Arguments

- `Formula{:unified}`: earth-return formula type.

# Returns

- `integral`: the default integration selection, `integrate => :quad`.
"""
function defaults(::Type{Formula{:unified}})
    return (integral = integrate => :quad,)
end

"""
$(TYPEDSIGNATURES)

Check the integration controls of the unified earth-return formula.

# Arguments

- `options`: formula options with an integration section.

# Returns

- `options`, unchanged.

# Errors

- Throws `ArgumentError` for a negative or nonfinite tolerance, two zero tolerances, or a
  nonpositive evaluation limit.
"""
function validate(options::FormulationOptions{Formula{:unified}})
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
- `formula`: unified earth-return formula that consumes the problem.

# Returns

- `problem`, unchanged.

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

Allocate cable matrices and state vectors under `unified`, then add the integral arrays for
the selected quadrature method.

# Arguments

- `formula`: unified earth-return formula with integration controls.
- `LineParametersProblem`: problem whose computation needs the arrays.
- `input`: numerical cable geometry.
- `buffers`: computation arrays already allocated.

# Returns

- `buffers` extended with the `unified` and `integral` blocks. [`prepare!`](@ref) fills the
  `unified` arrays.
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

Compute the circumference-current state and earth-return matrices at one frequency. For each
medium ``m``, form admittivity ``\\eta_m`` \\[S/m\\] and squared radial wavenumber ``k_m^2``
\\[1/m²\\]:

```math
\\eta_m = \\frac{1}{\\rho_m} + s\\varepsilon_0\\varepsilon_{r,m}, \\qquad
k_m^2 = s \\mu_0 \\mu_{r,m} \\eta_m - \\Gamma^2.
```

Here ``\\rho_m`` is resistivity \\[Ω·m\\], ``\\varepsilon_{r,m}`` and ``\\mu_{r,m}`` are
relative permittivity and permeability \\[dimensionless\\], and ``\\varepsilon_0`` and
``\\mu_0`` are [`vacuum_permittivity`](@ref) and [`vacuum_permeability`](@ref). Read the
longitudinal wavenumber ``\\Gamma`` from `buffers.gamma[k]`, in \\[1/m\\].

For a cable ``p`` of radius ``r_p`` in \\[m\\] and medium ``m``, set ``x_p=k_m r_p`` and
``b_p=|\\operatorname{Re}x_p|``. The scaled receiver factor is ``a_p=e^{-b_p}I_0(x_p)``,
with ``I_\\nu`` the modified Bessel function of order ``\\nu``. The factor for circumference
currents is ``F_p`` in \\[S·m\\]:

```math
F_p=2\\pi\\eta_m r_p^2
\\begin{cases}
\\frac12-\\frac{x_p^2}{16}+\\frac{x_p^4}{96}-\\frac{11x_p^6}{6144},
    & |\\operatorname{nominal}(x_p)|<10^{-3},\\\\
\\frac{I_1(x_p)}{x_p I_0(x_p)}, & \\text{otherwise}.
\\end{cases}
```

Evaluate magnetic, electric, and potential spectral integrals on the ray selected by
[`ray!`](@ref). Same-medium cable pairs also include direct and image Bessel terms. Air
receivers reference potential to the interface below them. Store the scaled impedance kernel
``K`` \\[Ω/m\\] and potential kernel ``H`` \\[m/F\\], then solve

```math
L_{pq}=\\frac{\\delta_{pq}}{a_p}-F_p K_{pq}, \\qquad
P_e L=H, \\qquad Z_e L=K+\\frac{\\Gamma^2}{s}H.
```

Here ``L`` is the dimensionless current matrix, ``\\delta_{pq}`` is the Kronecker delta,
``P_e`` holds earth-return potential coefficients \\[m/F\\], and ``Z_e`` holds earth-return
impedances \\[Ω/m\\]. Transposed LU solves apply the common right factor ``L^{-1}``.

# Arguments

- `formula`: unified earth-return formula with integration controls.
- `buffers`: computation arrays with `gamma`, `integral`, and `unified`.
- `input`: cable centers and outer radii.
- `media`: air and earth records with resistivity `rho` in \\[Ω·m\\] and relative `eps_r`,
  `mu_r` \\[dimensionless\\].
- `s`: complex frequency, `j 2π f`, \\[1/s\\].
- `k`: frequency index.

# Returns

- `buffers`, with `unified.Ze` \\[Ω/m\\], `unified.Pe` \\[m/F\\], and their shared state
  filled.

# Errors

- Throws `DomainError` for a nonfinite weighted earth kernel, a nonfinite integral, or a
  spectral ray whose geometric weight grows.
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

Copy the earth-return impedances computed by [`prepare!`](@ref) into the destination.

# Arguments

- `formula`: unified earth-return formula.
- `Zg`: destination for cable impedances in \\[Ω/m\\].
- `buffers`: computation arrays with the prepared `unified.Ze` matrix.
- `input`: numerical cable geometry and connections.
- `media`: air and earth properties at the frequency.
- `s`: complex frequency, `j 2π f`, \\[1/s\\].
- `k`: frequency index.

# Returns

- `Zg`, with earth-return impedances \\[Ω/m\\].
"""
function impedance!(
        formula::Formula{:unified}, Zg::AbstractMatrix, buffers, input, media, s, k)
    return copyto!(Zg, buffers.unified.Ze)
end

"""
$(TYPEDSIGNATURES)

Copy the earth-return potential coefficients computed by [`prepare!`](@ref) into the
destination.

# Arguments

- `formula`: unified earth-return formula.
- `Pg`: destination cable potential-coefficient matrix \\[m/F\\].
- `buffers`: computation arrays with the prepared `unified.Pe` matrix.
- `input`: numerical cable geometry and connections.
- `media`: air and earth properties at the frequency.
- `s`: complex frequency, `j 2π f`, \\[1/s\\].
- `k`: frequency index.

# Returns

- `Pg`, with earth-return potential coefficients \\[m/F\\].
"""
function potential!(
        formula::Formula{:unified}, Pg::AbstractMatrix, buffers, input, media, s, k)
    return copyto!(Pg, buffers.unified.Pe)
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
