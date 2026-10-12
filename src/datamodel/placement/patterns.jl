# The deferred maximum member count, `capacity()`. A physical tree resolves it from the
# member and its context.
struct _DeferredCardinality end

"""
Return the deferred maximum wire count used by placement patterns.
"""
capacity() = _DeferredCardinality()

"""
Return the maximum admissible member count for a placement declaration.
"""
function capacity end

"""
$(TYPEDEF)

Place `n` members on one circular locus. `n = capacity()` places as many members as fit,
and `r = nothing` takes the radius from the physical context. A physical tree resolves
both.

$(TYPEDFIELDS)
"""
struct Ring{N <: Union{Int, _DeferredCardinality}, T <: Real, R <: Union{Real, Nothing}}
    "Number of members \\[dimensionless\\]."
    n::N
    "Radius of the member-center locus \\[m\\]."
    r::R
    "Starting angle \\[rad\\]."
    φ0::T
    "Counter-clockwise angular span \\[rad\\]."
    span::T
    "Fractional clearance between adjacent members \\[dimensionless\\]."
    gap_frac::T

    function Ring{N, T, R}(n::N, r::R, φ0::T, span::T, gap_frac::T) where {
            N <: Union{Int, _DeferredCardinality}, T <: Real, R <: Union{Real, Nothing}}
        return validate(new{N, T, R}(n, r, φ0, span, gap_frac))
    end
end

function validate(pattern::Ring)
    validate(pattern.n, Ring)
    validate(pattern.r, Ring)
    isfinite(pattern.φ0) ||
        throw(DomainError(pattern.φ0, "ring start angle must be finite"))
    isfinite(pattern.span) && zero(pattern.span) < pattern.span <= oftype(pattern.span, 2π) ||
        throw(DomainError(pattern.span, "ring span must lie in (0, 2π]"))
    isfinite(pattern.gap_frac) && pattern.gap_frac >= zero(pattern.gap_frac) ||
        throw(DomainError(
            pattern.gap_frac, "ring gap fraction must be nonnegative and finite"
        ))
    return pattern
end

# A ring count is positive, or `capacity()`.
function validate(count::Int, ::Type{Ring})
    count > 0 || throw(ArgumentError("ring cardinality must be positive or capacity()"))
    return count
end
validate(count::_DeferredCardinality, ::Type{Ring}) = count

# A ring radius is nonnegative and finite, or `nothing`, which the context resolves.
function validate(radius::Real, ::Type{Ring})
    isfinite(radius) && radius >= zero(radius) ||
        throw(DomainError(radius, "ring radius must be nonnegative, finite, or contextual"))
    return radius
end
validate(radius::Nothing, ::Type{Ring}) = radius

"""
Maximum circular-wire count on a ring with fractional adjacent clearance.
"""
function capacity(
        ::Type{Ring},
        lay_radius::Real,
        item_radius::Real;
        gap_frac::Real = 0
)
    lay_radius > zero(lay_radius) || throw(DomainError(
        lay_radius, "lay radius must be positive"
    ))
    item_radius > zero(item_radius) || throw(DomainError(
        item_radius, "item radius must be positive"
    ))
    gap_frac >= zero(gap_frac) || throw(DomainError(
        gap_frac, "gap fraction must be nonnegative"
    ))
    ratio = item_radius * (one(gap_frac) + gap_frac) / lay_radius
    zero(ratio) < ratio < one(ratio) || return 0
    count = pi / asin(ratio)
    return max(0, floor(Int, count + 8eps(float(count))))
end

function Ring(n::Union{Int, _DeferredCardinality}, r::Real, φ0::Real, span::Real,
        gap_frac::Real)
    radius, angle, angular_span, gap = map(float, promote(r, φ0, span, gap_frac))
    return Ring{typeof(n), typeof(radius), typeof(radius)}(
        n, radius, angle, angular_span, gap)
end

function Ring(n::Union{Int, _DeferredCardinality}, r::Nothing, φ0::Real, span::Real,
        gap_frac::Real)
    angle, angular_span, gap = map(float, promote(φ0, span, gap_frac))
    return Ring{typeof(n), typeof(angle), Nothing}(n, r, angle, angular_span, gap)
end

# Any other integer count becomes an `Int`. `capacity()` always reaches a method above.
function Ring(n::Union{Integer, _DeferredCardinality}, r::Union{Real, Nothing}, φ0::Real,
        span::Real, gap_frac::Real)
    return Ring(Int(n), r, φ0, span, gap_frac)
end

function Ring(
        n::Union{Integer, _DeferredCardinality, AbstractGrid, Gridspace};
        r = nothing,
        φ0 = 0,
        span = 2π,
        gap_frac = 0,
        combine::Symbol = :product
)
    return parameterize(Ring, Ring, (n, r, φ0, span, gap_frac); combine)
end

"""
$(TYPEDEF)

Place members on `nr` concentric loci with `nφ` angular positions per locus.

$(TYPEDFIELDS)
"""
struct Polar{T <: Real}
    "Number of radial loci \\[dimensionless\\]."
    nr::Int
    "Angular positions per nonzero locus \\[dimensionless\\]."
    nφ::Int
    "First locus radius \\[m\\]."
    r0::T
    "Radial increment \\[m\\]."
    dr::T
    "Starting angle \\[rad\\]."
    φ0::T
    "Counter-clockwise angular span \\[rad\\]."
    span::T

    function Polar{T}(
            nr::Int, nφ::Int, r0::T, dr::T, φ0::T, span::T
    ) where {T <: Real}
        return validate(new{T}(nr, nφ, r0, dr, φ0, span))
    end
end

function validate(pattern::Polar)
    pattern.nr > 0 || throw(ArgumentError("polar radial count must be positive"))
    pattern.nφ > 0 || throw(ArgumentError("polar angular count must be positive"))
    isfinite(pattern.r0) && pattern.r0 >= zero(pattern.r0) || throw(DomainError(
        pattern.r0, "polar initial radius must be nonnegative and finite"))
    isfinite(pattern.dr) && pattern.dr >= zero(pattern.dr) || throw(DomainError(
        pattern.dr, "polar radial increment must be nonnegative and finite"))
    isfinite(pattern.φ0) ||
        throw(DomainError(pattern.φ0, "polar start angle must be finite"))
    isfinite(pattern.span) && zero(pattern.span) < pattern.span <= oftype(pattern.span, 2π) ||
        throw(DomainError(pattern.span, "polar span must lie in (0, 2π]"))
    return pattern
end

function Polar(
        nr::Union{Signed, Unsigned},
        nφ::Union{Signed, Unsigned},
        r0::Real,
        dr::Real,
        φ0::Real,
        span::Real
)
    values = map(float, promote(r0, dr, φ0, span))
    return Polar{typeof(first(values))}(Int(nr), Int(nφ), values...)
end

function Polar(;
        nr,
        nφ,
        r0 = 0,
        dr,
        φ0 = 0,
        span = 2π,
        combine::Symbol = :product
)
    return parameterize(Polar, Polar, (nr, nφ, r0, dr, φ0, span); combine)
end

"""
$(TYPEDEF)

Fill a circular domain with one central member and concentric tangent courses.

$(TYPEDFIELDS)
"""
struct Fill{T <: Real}
    "Outer radius available to the filled members \\[m\\]."
    r::T
    "Angular offset applied between consecutive courses \\[rad\\]."
    φ::T
    "Starting angle of the first outer course \\[rad\\]."
    φ0::T
    "Counter-clockwise angular span of each outer course \\[rad\\]."
    span::T

    Fill{T}(r::T, φ::T, φ0::T, span::T) where {T <: Real} = validate(new{T}(r, φ, φ0, span))
end

function validate(pattern::Fill)
    isfinite(pattern.r) && pattern.r > zero(pattern.r) || throw(DomainError(
        pattern.r, "fill radius must be positive and finite"
    ))
    isfinite(pattern.φ) || throw(DomainError(pattern.φ, "fill course offset must be finite"))
    isfinite(pattern.φ0) || throw(DomainError(pattern.φ0, "fill start angle must be finite"))
    isfinite(pattern.span) && zero(pattern.span) < pattern.span <= oftype(pattern.span, 2π) ||
        throw(DomainError(pattern.span, "fill span must lie in (0, 2π]"))
    return pattern
end

function Fill(r::Real, φ::Real, φ0::Real, span::Real)
    values = map(float, promote(r, φ, φ0, span))
    return Fill{typeof(first(values))}(values...)
end

function Fill(;
        r,
        φ,
        φ0 = 0,
        span = 2π,
        combine::Symbol = :product
)
    return parameterize(Fill, Fill, (r, φ, φ0, span); combine)
end

"""
$(TYPEDEF)

Place members on a centered rectangular lattice.

$(TYPEDFIELDS)
"""
struct Lattice{T <: Real}
    "Number of columns \\[dimensionless\\]."
    nx::Int
    "Number of rows \\[dimensionless\\]."
    ny::Int
    "Column spacing \\[m\\]."
    dx::T
    "Row spacing \\[m\\]."
    dy::T

    Lattice{T}(nx::Int, ny::Int, dx::T, dy::T) where {T <: Real} =
        validate(new{T}(nx, ny, dx, dy))
end

function validate(pattern::Lattice)
    pattern.nx > 0 || throw(ArgumentError("lattice column count must be positive"))
    pattern.ny > 0 || throw(ArgumentError("lattice row count must be positive"))
    (; nx, ny, dx, dy) = pattern
    isfinite(dx) && (nx == 1 ? dx >= zero(dx) : dx > zero(dx)) ||
        throw(DomainError(dx, "lattice column spacing is invalid"))
    isfinite(dy) && (ny == 1 ? dy >= zero(dy) : dy > zero(dy)) ||
        throw(DomainError(dy, "lattice row spacing is invalid"))
    return pattern
end

function Lattice(
        nx::Union{Signed, Unsigned},
        ny::Union{Signed, Unsigned},
        dx::Real,
        dy::Real
)
    values = map(float, promote(dx, dy))
    return Lattice{typeof(first(values))}(Int(nx), Int(ny), values...)
end

function Lattice(;
        nx,
        ny,
        dx,
        dy,
        combine::Symbol = :product
)
    return parameterize(Lattice, Lattice, (nx, ny, dx, dy); combine)
end

"""
Return the local member placements prescribed by a placement pattern: one
`(pose, primitive)` record per member. `primitive` is `nothing` when each member keeps its
own primitive, or the primitive that replaces it, such as a rectangle bent into a strip on
a ring.
"""
function placements end

# A ring places its members once its count and radius are resolved.
function placements(
        pattern::Union{Ring{_DeferredCardinality}, Ring{<:Any, <:Real, Nothing}}, item,
        compact)
    throw(ArgumentError("a ring with capacity() or a contextual radius places members " *
                        "only inside a physical tree, which resolves them"))
end

function _ring_step(pattern::Ring, count::Int)
    return count == 1 ? zero(pattern.span) :
           isapprox(pattern.span, oftype(pattern.span, 2π)) ?
           pattern.span / count : pattern.span / (count - 1)
end

function _ring_placements(pattern::Ring{Int, <:Real, <:Real})
    step = _ring_step(pattern, pattern.n)
    return [(pose = Pose2(
                 pattern.r * cos(pattern.φ0 + index * step),
                 pattern.r * sin(pattern.φ0 + index * step),
                 pattern.φ0 + index * step
             ), primitive = nothing)
            for index in 0:(pattern.n - 1)]
end

placements(::Nothing, item, ::Nothing) = [(pose = Pose2(0, 0, 0), primitive = nothing)]

# Ring members of this tangential width keep the requested gap along the ring.
function validate(pattern::Ring, ::typeof(placements), tangential_width::Real)
    pattern.n == 1 && return pattern
    interval_count = isapprox(pattern.span, oftype(pattern.span, 2π)) ?
                     pattern.n : pattern.n - 1
    chord = 2pattern.r * sin(pattern.span / interval_count / 2)
    required = tangential_width * (one(pattern.gap_frac) + pattern.gap_frac)
    tolerance = sqrt(eps(float(chord))) * max(one(chord), chord)
    chord + tolerance >= required || throw(DomainError(
        pattern.n,
        "ring members overlap or violate the requested gap: " *
        "count=$(pattern.n), radius=$(nominal(pattern.r)) m, " *
        "member width=$(nominal(tangential_width)) m, " *
        "gap fraction=$(nominal(pattern.gap_frac)), " *
        "required chord=$(nominal(required)) m, " *
        "available chord=$(nominal(chord)) m, " *
        "deficit=$(nominal(required - chord)) m"
    ))
    return pattern
end

function _tangential_width(definition::AbstractPrimitive)
    primitive = resolve(EmptyBoundary(), definition)
    return support(primitive, pi / 2) + support(primitive, -pi / 2)
end

function placements(pattern::Ring{Int, <:Real, <:Real}, item::AbstractPrimitive, ::Nothing)
    validate(pattern, placements, _tangential_width(item))
    return _ring_placements(pattern)
end

function placements(pattern::Ring{Int, <:Real, <:Real}, item::Disk, ::Nothing)
    validate(pattern, placements, 2item.r)
    return _ring_placements(pattern)
end

# Each rectangle bends into one strip on the ring, rotated into its angular slot.
function _strip_placements(pattern::Ring{Int, <:Real, <:Real}, inner, outer, angular_width)
    definition = BentStrip(inner, outer, angular_width)
    step = pattern.n == 1 ? zero(pattern.span) : pattern.span / pattern.n
    return [(pose = Pose2(0, 0, pattern.φ0 + index * step), primitive = definition)
            for index in 0:(pattern.n - 1)]
end

function placements(pattern::Ring{Int, <:Real, <:Real}, item::Rectangle, ::Nothing)
    occupied = pattern.n * item.w * (one(pattern.gap_frac) + pattern.gap_frac)
    available = pattern.r * pattern.span
    tolerance = sqrt(eps(float(available))) * max(one(available), available)
    occupied <= available + tolerance || throw(DomainError(
        pattern.n,
        "rectangular members exceed the available angular span"
    ))
    inner = max(zero(pattern.r), pattern.r - item.h / 2)
    outer = inner + item.h
    angular_width = 2area(item) / (outer^2 - inner^2)
    return _strip_placements(pattern, inner, outer, angular_width)
end

function placements(pattern::Ring{Int, <:Real, <:Real}, item::Sector, ::Nothing)
    iszero(pattern.r) &&
        return _origin_ring_placements(pattern, resolve(EmptyBoundary(), item), false)
    validate(pattern, placements, _tangential_width(item))
    return _ring_placements(pattern)
end

function placements(pattern::Ring{Int, <:Real, <:Real}, item::CableGeometry, ::Nothing)
    if iszero(pattern.r)
        shape = boundary(item)
        outer = findlast(region -> boundary(region.primitive) === shape, item.regions)
        # Insulated sector boundaries may touch, while bare sectors need a side clearance.
        insulated = outer !== nothing &&
                    item.regions[outer].source.material.kind === :insulator
        return _origin_ring_placements(pattern, shape, insulated)
    end
    validate(pattern, placements, 2support(boundary(item)))
    return _ring_placements(pattern)
end

# Members that share the ring's center are sectors whose sides face their neighbors.
function _origin_ring_placements(pattern::Ring, shape::SectorShape, allow_contact::Bool)
    pattern.n == 1 && return _ring_placements(pattern)
    iszero(pattern.gap_frac) || throw(DomainError(
        pattern.gap_frac,
        "origin-centered sector spacing is set by its physical side clearance, " *
        "not Ring gap_frac"
    ))
    pitch = _ring_step(pattern, pattern.n)
    angle_tolerance = sqrt(eps(float(_geometry_scalar(pitch)))) *
                      max(one(_geometry_scalar(pitch)), abs(_geometry_scalar(pitch)))
    isapprox(
        _geometry_scalar(shape.primitive.span),
        _geometry_scalar(pitch);
        atol = angle_tolerance,
        rtol = sqrt(eps(float(_geometry_scalar(pitch))))
    ) || throw(DomainError(
        shape.primitive.span,
        "an origin-centered sector span must equal its angular pitch $pitch so " *
        "neighboring wedge sides remain parallel"
    ))
    clearance = _sector_side_clearance(shape.primitive)
    clearance_tolerance = sqrt(eps(float(_geometry_scalar(shape.primitive.r_back)))) *
                          max(
        one(_geometry_scalar(shape.primitive.r_back)),
        abs(_geometry_scalar(shape.primitive.r_back))
    )
    valid_clearance = allow_contact ?
                      _geometry_scalar(clearance) >= -clearance_tolerance :
                      _geometry_scalar(clearance) > clearance_tolerance
    valid_clearance || throw(DomainError(
        clearance,
        allow_contact ?
        "sector insulation boundaries overlap" :
        "bare origin-centered sectors require a positive physical side clearance"
    ))
    return _ring_placements(pattern)
end

function _origin_ring_placements(pattern::Ring, shape, allow_contact::Bool)
    throw(ArgumentError("only sector-shaped members can share an origin"))
end

function placements(pattern::Polar, item, ::Nothing)
    courses = map(0:(pattern.nr - 1)) do radial_index
        radius = pattern.r0 + radial_index * pattern.dr
        iszero(radius) && return [(pose = Pose2(0, 0, pattern.φ0), primitive = nothing)]
        return placements(
            Ring(pattern.nφ; r = radius, φ0 = pattern.φ0, span = pattern.span),
            item,
            nothing
        )
    end
    return reduce(vcat, courses)
end

function _ring_capacity(pattern::Ring, ratio::Real)
    zero(ratio) < ratio < one(ratio) || return 0
    half_angle = asin(ratio)
    if isapprox(pattern.span, oftype(pattern.span, 2π))
        count = pi / half_angle
        return max(0, floor(Int, count + 8eps(float(count))))
    end
    count = pattern.span / (2half_angle)
    return max(0, floor(Int, count + 8eps(float(count))) + 1)
end

function capacity(pattern::Ring, item::Disk, ::Nothing)
    radius = something(
        pattern.r,
        item.r * (one(item.r) + pattern.gap_frac)
    )
    ratio = item.r * (one(pattern.gap_frac) + pattern.gap_frac) / radius
    return _ring_capacity(pattern, ratio)
end

function capacity(pattern::Ring, item::Rectangle, ::Nothing)
    radius = something(pattern.r, item.h / 2)
    radius > zero(radius) || return 0
    ratio = item.w * (one(pattern.gap_frac) + pattern.gap_frac) / (2radius)
    return _ring_capacity(pattern, ratio)
end

function capacity(pattern::Ring, item::Sector, ::Nothing)
    if pattern.r !== nothing && !iszero(pattern.r)
        width = _tangential_width(item)
        ratio = width * (one(pattern.gap_frac) + pattern.gap_frac) /
                (2pattern.r)
        return _ring_capacity(pattern, ratio)
    end
    occupied = item.span * (one(pattern.gap_frac) + pattern.gap_frac)
    return max(0, floor(Int, pattern.span / occupied + 8eps(float(pattern.span))))
end

function capacity(pattern::Ring, item::CableGeometry, ::Nothing)
    member_radius = support(boundary(item))
    radius = something(
        pattern.r,
        member_radius * (one(member_radius) + pattern.gap_frac)
    )
    ratio = member_radius * (one(pattern.gap_frac) + pattern.gap_frac) / radius
    return _ring_capacity(pattern, ratio)
end

function placements(pattern::Lattice, item, ::Nothing)
    x0 = -(pattern.nx - 1) * pattern.dx / 2
    y0 = -(pattern.ny - 1) * pattern.dy / 2
    return [(pose = Pose2(x0 + ix * pattern.dx, y0 + iy * pattern.dy, 0), primitive = nothing)
            for iy in 0:(pattern.ny - 1) for ix in 0:(pattern.nx - 1)]
end

_fill_member_extent(item::AbstractPrimitive) = support(resolve(EmptyBoundary(), item))
_fill_member_extent(item::CableGeometry) = support(boundary(item))

function _fill_placements(pattern::Fill, item, compact)
    extent = _fill_member_extent(item)
    extent > zero(extent) || throw(DomainError(
        extent, "fill members must have positive radial extent"
    ))
    tolerance = sqrt(eps(float(pattern.r))) * max(one(pattern.r), pattern.r)
    extent <= pattern.r + tolerance || throw(DomainError(
        extent, "one member does not fit inside the fill radius"
    ))

    # Courses may resolve to another scalar type than the center, as the poses did.
    placed = NamedTuple{(:pose, :primitive)}[(pose = Pose2(0, 0, pattern.φ0), primitive = nothing)]
    course = 1
    radius = 2extent
    while radius + extent <= pattern.r + tolerance
        ring = Ring(
            capacity();
            r = radius,
            φ0 = pattern.φ0 + (course - 1) * pattern.φ,
            span = pattern.span
        )
        count = capacity(ring, item, compact)
        count > 0 || break
        append!(placed, placements(
            Ring(count; r = ring.r, φ0 = ring.φ0, span = ring.span),
            item,
            compact
        ))
        course += 1
        radius = 2course * extent
    end
    return placed
end

placements(pattern::Fill, item, ::Nothing) = _fill_placements(pattern, item, nothing)
capacity(pattern::Fill, item, ::Nothing) = length(_fill_placements(pattern, item, nothing))

function placements(poses::AbstractVector{<:Pose2}, item, ::Nothing)
    isempty(poses) && throw(ArgumentError("explicit placements cannot be empty"))
    return [(pose = pose, primitive = nothing) for pose in poses]
end
