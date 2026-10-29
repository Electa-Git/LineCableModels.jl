"""
$(TYPEDEF)

Represent the integral of a complete integrand over ``[0, \\infty)``.
The function includes the weights and Jacobians of its formula.

$(TYPEDFIELDS)
"""
struct SpectralIntegral{F <: Function}
    "Complete integrand, a function of the nonnegative real coordinate."
    integrand::F
end

"""
$(TYPEDSIGNATURES)

Evaluate an integral by adaptive Gauss–Kronrod quadrature with QuadGK. Transform its
nonnegative coordinate ``\\lambda`` into ``t``:

```math
\\int_0^\\infty g(\\lambda)\\,\\mathrm{d}\\lambda
= \\int_0^1 g\\!\\left(\\frac{t}{1-t}\\right)\\frac{\\mathrm{d}t}{(1-t)^2}.
```

Here ``g`` is the complete integrand. Map each breakpoint ``\\lambda_b`` to
``t_b = \\lambda_b / (1 + \\lambda_b)`` in ``[0, 1)`` and seed one segment between adjacent
points. The error norm is the magnitude of the nominal value.

# Arguments

- `integral`: complete integrand and its integral.
- `Val(:quad)`: adaptive quadrature method tag.
- `section`: controls `rtol`, `atol`, and `maxevals` from [`defaults`](@ref).
- `buffers`: computation arrays with the `integral` block.

# Keywords

- `points`: finite, nonnegative breakpoints in the coordinate of the integrand.

# Returns

- The integral value, in the unit of the integrand times its coordinate.
- The estimated absolute error, in the unit of the integral value.

# Notes

An estimated error above `max(atol, rtol * abs(nominal(value)))`, or a nonfinite estimate,
produces a warning. The method still returns the value.

# Errors

- Throws `DomainError` for a negative or nonfinite breakpoint, or a nonfinite integral
  value.
"""
function integrate(integral::SpectralIntegral, ::Val{:quad}, section::NamedTuple, buffers;
        points)
    block = buffers.integral
    R = eltype(block.mapped)
    mapped = empty!(block.mapped)
    push!(mapped, zero(R))
    for point in points
        isfinite(point) && point >= 0 ||
            throw(DomainError(point, "integration points must be finite nonnegative real coordinates"))
        push!(mapped, R(point) / (one(R) + R(point)))
    end
    push!(mapped, one(R))
    unique!(sort!(mapped; alg = Base.Sort.QuickSort))
    seeds = empty!(block.seeds)
    sample = block.prototype
    for i in 1:(length(mapped) - 1)
        quadgk(_ -> sample, mapped[i], mapped[i + 1];
            segbuf = block.seed, norm = z -> abs(nominal(z)))
        append!(seeds, block.seed)
    end
    value,
    estimate = quadgk(t -> begin
            scale = inv(one(R) - t)
            integral.integrand(t * scale) * scale^2
        end, zero(R),
        one(R); rtol = section.rtol, atol = section.atol,
        maxevals = section.maxevals, segbuf = block.segments,
        eval_segbuf = seeds, norm = z -> abs(nominal(z)))
    isfinite(value) || throw(DomainError(value, "QuadGK returned a nonfinite integral"))
    target = max(section.atol, section.rtol * abs(nominal(value)))
    if !isfinite(estimate) || estimate > target
        @warn("QuadGK returned an estimated error above the requested target", value,
            estimate, target, rtol=section.rtol, atol=section.atol,
            maxevals=section.maxevals)
    end
    return value, estimate
end

"""
$(TYPEDSIGNATURES)

Return the declared controls of adaptive Gauss–Kronrod quadrature.

# Arguments

- `integrate`: integral evaluation function.
- `Val(:quad)`: adaptive quadrature method tag.

# Returns

- `method`: `:quad`.
- `rtol`: relative tolerance, `1e-8`.
- `atol`: absolute tolerance in the unit of the integral, `0.0`.
- `maxevals`: evaluation limit, `10^7`.
"""
function defaults(::typeof(integrate), ::Val{:quad})
    return (method = :quad, rtol = 1e-8, atol = 0.0, maxevals = 10^7)
end

# A method without a tag of its own has no section.
function defaults(::typeof(integrate), ::Val{M}) where {M}
    throw(ArgumentError("unknown integration method :$M; the methods are (:quad,)"))
end

"""
$(TYPEDSIGNATURES)

Check the quadrature controls before an integral evaluation.

# Arguments

- `section`: integration controls.
- `integrate`: consumer of the controls.
- `Val(:quad)`: adaptive quadrature method tag.

# Returns

- `section`, unchanged.

# Errors

- Throws `ArgumentError` for a negative or nonfinite tolerance, two zero tolerances, or a
  nonpositive evaluation limit.
"""
function validate(section::NamedTuple, ::typeof(integrate), ::Val{:quad})
    isfinite(section.rtol) && section.rtol >= 0 ||
        throw(ArgumentError("rtol must be finite and nonnegative"))
    isfinite(section.atol) && section.atol >= 0 ||
        throw(ArgumentError("atol must be finite and nonnegative"))
    section.rtol > 0 || section.atol > 0 ||
        throw(ArgumentError("at least one integration tolerance must be positive"))
    section.maxevals > 0 || throw(ArgumentError("maxevals must be positive"))
    return section
end

"""
$(TYPEDSIGNATURES)

Allocate the integral arrays. Coordinates use the type of the nominal value of `T`, and
integrand values use `Complex{T}`.

# Arguments

- `SpectralIntegral`: integral type whose arrays the computation needs.
- `Val(:quad)`: adaptive quadrature method tag.
- `T`: numerical type of the formula.
- `buffers`: arrays already allocated.

# Returns

- `buffers` extended with the `integral` block. It holds breakpoint and scale vectors plus
  mapped points, QuadGK segment arrays, and a zero prototype.
"""
function initialize_buffers(::Type{SpectralIntegral}, ::Val{:quad}, ::Type{T}, buffers) where {T}
    R = typeof(nominal(one(T)))
    return merge(buffers,
        (integral = (
            points = sizehint!(R[], 128),
            scales = sizehint!(R[], 16),
            segments = alloc_segbuf(R, Complex{T}, R; size = 128),
            seeds = alloc_segbuf(R, Complex{T}, R; size = 128),
            seed = alloc_segbuf(R, Complex{T}, R; size = 1),
            prototype = zero(Complex{T}),
            mapped = sizehint!(R[], 128)),))
end

"""
$(TYPEDSIGNATURES)

Build an integration section from method defaults and supplied controls, then replace
its method identifier with a `Val` tag.

# Arguments

- `default`: pair of `integrate` and the default method identifier.
- `given`: supplied controls, with or without a `method` key.

# Returns

- A named tuple with the method keys and their declared value types.

# Errors

- Throws `ArgumentError` for an unknown method or an unknown control key.
"""
function select(default::Pair{typeof(integrate), Symbol}, given::NamedTuple)
    tag = Val(merge((method = last(default),), given).method)
    return merge(select(defaults(integrate, tag), given), (method = tag,))
end

# The declared default stands for a section that the user did not give.
function select(default::Pair{typeof(integrate), Symbol},
        ::Pair{typeof(integrate), Symbol})
    return select(default, (;))
end

"""
$(TYPEDSIGNATURES)

Choose a complex spectral ray and write quadrature breakpoints in its dimensionless
coordinate ``t``. The spectral coordinate is

```math
\\lambda = R e^{\\mathrm{j}\\theta}t, \\qquad
R = \\max\\!\\left(\\max_m |k_m|,\\frac{1}{h}\\right), \\qquad
h\\cos\\theta-y\\sin\\theta > 0.
```

Here ``h`` is `height`, ``y`` is `reach`, and ``k_m`` are medium wavenumbers. The angle is

```math
\\theta = \\min\\left(\\frac{\\pi}{6},
    \\arctan\\frac{h}{2 \\max(y, \\epsilon)},
    \\min_{\\zeta} \\frac{\\arg \\zeta}{2}\\right),
```

where ``\\epsilon`` is the machine epsilon of the coordinate type and ``\\zeta`` runs over the
branch points ``\\pm\\mathrm{j}k_m`` and the poles in the open first quadrant. These
choices make the geometric weight ``e^{-h\\lambda}\\cos(y\\lambda)`` decay along the ray.

Use nominal values for all decisions. The scales are nonzero wavenumber magnitudes, finite
positive formula scales, and the reciprocal geometric decay. Each scale contributes half,
once, and twice its value. Each singularity contributes up to five positive points around
its projection. Bridge gaps by factors of four, divide points by ``R``, include zero and
one, then sort and remove duplicates.

# Arguments

- `buffers`: computation arrays with the `integral` block.
- `wavenumbers`: tuple of medium wavenumbers \\[1/m\\].
- `poles`: tuple of poles of the integrand \\[1/m\\]. Ignore zero entries.
- `scales`: tuple of formula-specific spectral scales \\[1/m\\]. Only finite positive values
  contribute.
- `height`: sum of source and receiver distances from the interface \\[m\\].
- `reach`: horizontal separation, with any receiver-radius allowance needed by the integrand
  \\[m\\].

# Returns

- The complex ray ``R e^{\\mathrm{j}\\theta}``, in \\[1/m\\]. The method also fills
  `buffers.integral.points` with dimensionless breakpoints.

# Errors

- Throws `DomainError` when the geometric weight does not decay along the ray.
"""
function ray!(buffers, wavenumbers::Tuple, poles::Tuple, scales::Tuple, height, reach)
    block = buffers.integral
    R = eltype(block.points)
    angle = min(R(π) / 6, atan(R(nominal(height)) / (2max(R(nominal(reach)), eps(R)))))
    locations = (map(k -> nominal(-im * k), wavenumbers)...,
        map(k -> nominal(im * k), wavenumbers)..., map(nominal, poles)...)
    for location in locations
        real(location) > 0 && imag(location) > 0 || continue
        angle = min(angle, atan(imag(location), real(location)) / 2)
    end
    angle = R(angle)
    empty!(block.scales)
    for k in wavenumbers
        iszero(nominal(k)) || push!(block.scales, R(abs(nominal(k))))
    end
    for scale in scales
        value = R(nominal(scale))
        isfinite(value) && value > 0 && push!(block.scales, value)
    end
    decay = R(nominal(height)) * cos(angle) - R(nominal(reach)) * sin(angle)
    decay > 0 ||
        throw(DomainError(angle, "spectral contour has a growing geometric weight"))
    push!(block.scales, inv(decay))
    unique!(sort!(block.scales; alg = Base.Sort.QuickSort))
    points = empty!(block.points)
    push!(points, zero(R))
    for scale in block.scales
        append!(points, (scale / 2, scale, 2scale))
    end
    rotation = cis(angle)
    for location in locations
        projected = location / rotation
        centre = R(real(projected))
        centre > 0 || continue
        width = max(R(abs(imag(projected))), 64eps(R) * centre)
        for offset in (-4, -1, 0, 1, 4)
            point = centre + offset * width
            point > 0 && isfinite(point) && push!(points, point)
        end
    end
    unique!(sort!(points; alg = Base.Sort.QuickSort))
    n = length(points)
    for i in 2:(n - 1)
        point = 4points[i]
        while point < points[i + 1]
            push!(points, point)
            point *= 4
        end
    end
    magnitude = max(maximum(k -> R(abs(nominal(k))), wavenumbers), inv(R(nominal(height))))
    points ./= magnitude
    push!(points, one(R))
    unique!(sort!(points; alg = Base.Sort.QuickSort))
    return magnitude * cis(angle)
end
