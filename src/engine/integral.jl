"""
$(TYPEDEF)

The integral of a function over ``[0, \\infty)``. The function is the complete integrand,
its weights and Jacobians included, and `integrate` evaluates the integral by the method
that its tag selects.

$(TYPEDFIELDS)
"""
struct SpectralIntegral{F <: Function}
    "Complete integrand, a function of the nonnegative real coordinate."
    integrand::F
end

"""
$(TYPEDSIGNATURES)

Evaluate `integral` over ``[0, \\infty)`` by the adaptive Gauss–Kronrod method of QuadGK,
with the tolerances and the evaluation limit of `section`. The change of variable

```math
\\int_0^\\infty f(\\lambda) \\, \\mathrm{d}\\lambda
= \\int_0^1 f\\!\\left(\\frac{t}{1 - t}\\right) \\frac{\\mathrm{d}t}{(1 - t)^2}
```

maps the breakpoints `points` into ``[0, 1)``, and each interval between two breakpoints
seeds one segment, so that the refinement starts where the integrand changes character. The
error norm is the magnitude of the nominal value.

# Arguments

- `integral`: the integral to evaluate.
- `section`: the section of the method, with the keys of `defaults(integrate, Val(:quad))`.
- `buffers`: arrays of the computation, with the engine's integral block `integral`.

# Keywords

- `points`: breakpoints of the integrand in its own coordinate, finite and nonnegative.

# Returns

- The value of the integral.
- The estimated absolute error, in the unit of the value.

# Notes

- An estimated error above `max(atol, rtol * abs(value))` gives a warning, and `integrate`
  still returns the value.

# Errors

- Throws `DomainError` for a breakpoint that is negative or not finite, and for a value
  that is not finite.
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

Return the section of the adaptive Gauss–Kronrod method `:quad` with its default values. Its
keys are `method`, the relative tolerance `rtol`, the absolute tolerance `atol` in the unit
of the integral, and the evaluation limit `maxevals`.
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

Check the section of the method `:quad` for the integrals that `integrate` evaluates.

# Errors

- Throws `ArgumentError` when `rtol` or `atol` is negative or not finite, when neither is
  positive, or when `maxevals` is not positive.
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

Return `buffers` with the engine's integral block under the key `integral`. Every formula
that integrates declares this block in its `initialize_buffers`. The block has:

- `points` and `scales`, the breakpoints and the scales that `ray!` writes.
- `segments`, `seeds` and `seed`, the segments of QuadGK and the segments between
  breakpoints.
- `prototype`, a value of the type of the integrand.
- `mapped`, the breakpoints after the change of variable.

The coordinates have the type of the nominal value of `T`, and the integrand values have the
type `Complex{T}`.
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

Lift the integration section of a formula's options. The default `integrate => :quad`
names the function whose method the section selects, and the default method. The tuple
`given` has the keys that the user gave, with or without `method`. The lifted section has
the keys of the method in their declared types, and `method` becomes its tag.

# Errors

- Throws `ArgumentError` for an unknown method or an unknown key of the method.
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

Return the ray of a spectral integral over ``[0, \\infty)``, the complex number `ray`
such that the integrand takes ``\\lambda = \\mathrm{ray} \\, t`` for ``t \\ge 0``. Write
into `buffers.integral.points` the coordinates ``t`` where the quadrature subdivides the
integral.

The ray leaves the real axis at an angle of at most π/6, so that the geometric factor
``e^{-h \\lambda} \\cos(y \\lambda)`` still decays along it:
``h \\cos\\theta - y \\sin\\theta > 0``, with ``y`` the `reach`. The ray also passes below
each singularity in the first quadrant. These are the branch points
``\\pm \\mathrm{j} k`` of the `wavenumbers` of the media and the `poles` of the integrand,
where a zero pole stands for none. The modulus of the ray is the largest of the moduli of
the wavenumbers and ``1 / h``. On the ray, the integrand takes its roots
``\\sqrt{\\lambda^2 + k^2}`` by Julia's `sqrt`, whose real part is nonnegative. That is the
sheet on which a field of the engine's time convention ``e^{\\mathrm{j} \\omega t}``
decays away from its source.

The points mark where the integrand changes character. The adaptive Gauss–Kronrod rule
starts its segments at them, so that its refinement begins where the integrand demands it.
A scale ``s`` contributes the points ``s/2``, ``s`` and ``2s``. The scales are the moduli
of the wavenumbers and the `scales` of the formula, with the reciprocal of the geometric
decay along the ray. Each singularity adds five points around its projection on the ray.
Points at factors of four bridge each gap between consecutive points.

`ray!` sorts the points and keeps each once, 1 among them. Without them the rule can step
over a feature narrower than its nodes and return a converged wrong value. For the
earth-return integrals with a longitudinal wavenumber, the default tolerance then misses by
up to 4.6e-8 relative, against 4e-12 with them. This function is the one tool for that
purpose: waypoints, complex paths and safe poles are other names for it.

# Arguments

- `buffers`: arrays of the computation, with the engine's integral block `integral`.
- `wavenumbers`: the wavenumbers ``k`` of the media \\[1/m\\], a tuple.
- `poles`: a tuple of the poles of the integrand \\[1/m\\]. A zero entry is no pole.
- `scales`: a tuple of the scales that the formula adds \\[1/m\\]. Only a finite positive
  scale counts.
- `height`: the sum of the heights of the receiver and the source above the interface
  \\[m\\].
- `reach`: the horizontal separation of the receiver and the source plus the radius of the
  receiver \\[m\\].

# Errors

- Throws `DomainError` when the geometric factor does not decay along the ray.
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
