"""
$(TYPEDEF)

Place `n` members on one circular locus.

$(TYPEDFIELDS)
"""
struct _DeferredCardinality end

"Return the deferred maximum-capacity policy used by placement patterns."
capacity() = _DeferredCardinality()

"Return the maximum admissible member count for a placement declaration."
function capacity end

struct Ring{N, T <: Real, R}
    "Number of members \\[dimensionless\\]."
    n::N
    "Radius of the member-centre locus \\[m\\]."
    r::R
    "Starting angle \\[rad\\]."
    φ0::T
    "Counter-clockwise angular span \\[rad\\]."
    span::T
    "Fractional clearance between adjacent members \\[dimensionless\\]."
    gap_frac::T

    function Ring{N, T, R}(
            n::N, r::R, φ0::T, span::T, gap_frac::T
    ) where {N, T <: Real, R}
        n isa _DeferredCardinality ||
            (n isa Int && n > 0) ||
            throw(ArgumentError("ring cardinality must be positive or capacity()"))
        r === nothing ||
            (r isa Real && isfinite(r) && r >= zero(r)) ||
            throw(DomainError(r, "ring radius must be nonnegative, finite, or contextual"))
        isfinite(φ0) || throw(DomainError(φ0, "ring start angle must be finite"))
        isfinite(span) && zero(span) < span <= oftype(span, 2π) ||
            throw(DomainError(span, "ring span must lie in (0, 2π]"))
        isfinite(gap_frac) && gap_frac >= zero(gap_frac) || throw(DomainError(
            gap_frac, "ring gap fraction must be nonnegative and finite"
        ))
        return new{N, T, R}(n, r, φ0, span, gap_frac)
    end
end

"Maximum circular-wire count on a ring with fractional adjacent clearance."
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

function _ring(n, r, φ0, span, gap_frac)
    n isa Union{Integer, _DeferredCardinality} || throw(ArgumentError(
        "ring cardinality must be a positive integer or capacity()"
    ))
    numeric = r === nothing ? promote(φ0, span, gap_frac) :
              promote(r, φ0, span, gap_frac)
    values = map(float, numeric)
    if r === nothing
        angle, angular_span, gap = values
        return Ring{typeof(n), typeof(angle), Nothing}(
            n isa Integer ? Int(n) : n,
            nothing,
            angle,
            angular_span,
            gap
        )
    end
    radius, angle, angular_span, gap = values
    return Ring{typeof(n isa Integer ? Int(n) : n), typeof(radius), typeof(radius)}(
        n isa Integer ? Int(n) : n,
        radius,
        angle,
        angular_span,
        gap
    )
end

function Ring(
        n;
        r = nothing,
        φ0 = 0,
        span = 2π,
        gap_frac = 0,
        combine::Symbol = :product
)
    return parameterize(
        Ring, _ring, (n, r, φ0, span, gap_frac); combine
    )
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
        nr > 0 || throw(ArgumentError("polar radial count must be positive"))
        nφ > 0 || throw(ArgumentError("polar angular count must be positive"))
        isfinite(r0) && r0 >= zero(r0) ||
            throw(DomainError(r0, "polar initial radius must be nonnegative and finite"))
        isfinite(dr) && dr >= zero(dr) ||
            throw(DomainError(dr, "polar radial increment must be nonnegative and finite"))
        isfinite(φ0) || throw(DomainError(φ0, "polar start angle must be finite"))
        isfinite(span) && zero(span) < span <= oftype(span, 2π) ||
            throw(DomainError(span, "polar span must lie in (0, 2π]"))
        return new{T}(nr, nφ, r0, dr, φ0, span)
    end
end

function _polar(nr, nφ, r0, dr, φ0, span)
    nr isa Integer && !(nr isa Bool) || throw(ArgumentError(
        "polar radial count must be an integer"
    ))
    nφ isa Integer && !(nφ isa Bool) || throw(ArgumentError(
        "polar angular count must be an integer"
    ))
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
    return parameterize(Polar, _polar, (nr, nφ, r0, dr, φ0, span); combine)
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

    function Fill{T}(r::T, φ::T, φ0::T, span::T) where {T <: Real}
        isfinite(r) && r > zero(r) || throw(DomainError(
            r, "fill radius must be positive and finite"
        ))
        isfinite(φ) || throw(DomainError(φ, "fill course offset must be finite"))
        isfinite(φ0) || throw(DomainError(φ0, "fill start angle must be finite"))
        isfinite(span) && zero(span) < span <= oftype(span, 2π) || throw(
            DomainError(span, "fill span must lie in (0, 2π]")
        )
        return new{T}(r, φ, φ0, span)
    end
end

function _fill(r, φ, φ0, span)
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
    return parameterize(Fill, _fill, (r, φ, φ0, span); combine)
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

    function Lattice{T}(nx::Int, ny::Int, dx::T, dy::T) where {T <: Real}
        nx > 0 || throw(ArgumentError("lattice column count must be positive"))
        ny > 0 || throw(ArgumentError("lattice row count must be positive"))
        isfinite(dx) && (nx == 1 ? dx >= zero(dx) : dx > zero(dx)) ||
            throw(DomainError(dx, "lattice column spacing is invalid"))
        isfinite(dy) && (ny == 1 ? dy >= zero(dy) : dy > zero(dy)) ||
            throw(DomainError(dy, "lattice row spacing is invalid"))
        return new{T}(nx, ny, dx, dy)
    end
end

function _lattice(nx, ny, dx, dy)
    nx isa Integer && !(nx isa Bool) || throw(ArgumentError(
        "lattice column count must be an integer"
    ))
    ny isa Integer && !(ny isa Bool) || throw(ArgumentError(
        "lattice row count must be an integer"
    ))
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
    return parameterize(Lattice, _lattice, (nx, ny, dx, dy); combine)
end

"Return local member poses prescribed by a placement pattern."
function placements end

struct _ResolvedPlacement{P <: Pose2, D <: AbstractShape}
    pose::P
    primitive::D
end

_placement_pose(value::Pose2) = value
_placement_pose(value::_ResolvedPlacement) = value.pose
_placement_definition(value::Pose2, definition::AbstractPrimitive) = definition
_placement_definition(value::_ResolvedPlacement, ::AbstractPrimitive) = value.primitive

function _ring_step(pattern::Ring, count::Int)
    return count == 1 ? zero(pattern.span) :
           isapprox(pattern.span, oftype(pattern.span, 2π)) ?
           pattern.span / count : pattern.span / (count - 1)
end

function _ring_poses(pattern::Ring, count::Int, radius::Real)
    step = _ring_step(pattern, count)
    return Pose2[Pose2(
                     radius * cos(pattern.φ0 + index * step),
                     radius * sin(pattern.φ0 + index * step),
                     pattern.φ0 + index * step
                 )
                 for index in 0:(count - 1)]
end

placements(::Nothing, item, ::Nothing) = Pose2[Pose2(0, 0, 0)]

function _check_ring_clearance(pattern::Ring, tangential_width::Real)
    pattern.n == 1 && return nothing
    interval_count = isapprox(pattern.span, oftype(pattern.span, 2π)) ?
                     pattern.n : pattern.n - 1
    chord = 2pattern.r * sin(pattern.span / interval_count / 2)
    required = tangential_width * (one(pattern.gap_frac) + pattern.gap_frac)
    tolerance = sqrt(eps(float(chord))) * max(one(chord), chord)
    chord + tolerance >= required || throw(DomainError(
        pattern.n,
        "ring members overlap or violate the requested gap"
    ))
    return nothing
end

function _tangential_width(definition::AbstractPrimitive)
    primitive = resolve(EmptyBoundary(), definition)
    return support(primitive, pi / 2) + support(primitive, -pi / 2)
end

function placements(
        pattern::Ring,
        item::AbstractPrimitive,
        ::Nothing
)
    pattern.n isa Int || throw(ArgumentError(
        "capacity() requires contextual placement"
    ))
    pattern.r === nothing && throw(ArgumentError(
        "a contextual ring radius requires placement inside a physical tree"
    ))
    _check_ring_clearance(pattern, _tangential_width(item))
    return _ring_poses(pattern, pattern.n, pattern.r)
end

function placements(pattern::Ring, item::Disk, ::Nothing)
    pattern.n isa Int || throw(ArgumentError(
        "capacity() requires contextual placement"
    ))
    pattern.r === nothing && throw(ArgumentError(
        "a contextual ring radius requires placement inside a physical tree"
    ))
    _check_ring_clearance(pattern, 2item.r)
    return _ring_poses(pattern, pattern.n, pattern.r)
end

function placements(pattern::Ring, item::Rectangle, ::Nothing)
    pattern.n isa Int || throw(ArgumentError(
        "capacity() requires contextual placement"
    ))
    pattern.r === nothing && throw(ArgumentError(
        "a contextual ring radius requires placement inside a physical tree"
    ))
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
    definition = BentStrip(inner, outer, angular_width)
    step = pattern.n == 1 ? zero(pattern.span) : pattern.span / pattern.n
    return _ResolvedPlacement[
        _ResolvedPlacement(
            Pose2(0, 0, pattern.φ0 + index * step),
            definition
        ) for index in 0:(pattern.n - 1)
    ]
end

function placements(pattern::Ring, item::Sector, ::Nothing)
    pattern.n isa Int || throw(ArgumentError(
        "capacity() requires contextual placement"
    ))
    pattern.r === nothing && throw(ArgumentError(
        "a contextual ring radius requires placement inside a physical tree"
    ))
    if iszero(pattern.r)
        return _origin_ring_poses(
            pattern,
            resolve(EmptyBoundary(), item)
        )
    end
    _check_ring_clearance(pattern, _tangential_width(item))
    return _ring_poses(pattern, pattern.n, pattern.r)
end

function placements(pattern::Ring, item::CableGeometry, ::Nothing)
    pattern.n isa Int || throw(ArgumentError(
        "capacity() requires contextual placement"
    ))
    pattern.r === nothing && throw(ArgumentError(
        "a contextual ring radius requires placement inside a physical tree"
    ))
    iszero(pattern.r) && return _origin_ring_poses(pattern, item)
    _check_ring_clearance(pattern, 2support(boundary(item)))
    return _ring_poses(pattern, pattern.n, pattern.r)
end

function _sector_ring_poses(pattern::Ring, shape::SectorShape, allow_contact::Bool)
    pattern.n == 1 && return _ring_poses(pattern, pattern.n, pattern.r)
    iszero(pattern.gap_frac) || throw(DomainError(
        pattern.gap_frac,
        "origin-centred sector spacing is set by its physical side clearance, " *
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
        "an origin-centred sector span must equal its angular pitch $pitch so " *
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
        "bare origin-centred sectors require a positive physical side clearance"
    ))
    return _ring_poses(pattern, pattern.n, pattern.r)
end

_origin_ring_poses(pattern::Ring, shape::SectorShape) =
    _sector_ring_poses(pattern, shape, false)

function _origin_ring_poses(pattern::Ring, item::CableGeometry)
    shape = boundary(item)
    shape isa SectorShape || throw(ArgumentError(
        "only sector-shaped members can share an origin"
    ))
    outer_region = findlast(
        region -> boundary(region.primitive) === shape,
        item.regions
    )
    allow_contact = outer_region !== nothing &&
                    item.regions[outer_region].source.material.kind === :insulator
    return _sector_ring_poses(pattern, shape, allow_contact)
end

function placements(pattern::Polar, item, ::Nothing)
    poses = Pose2[]
    for radial_index in 0:(pattern.nr - 1)
        radius = pattern.r0 + radial_index * pattern.dr
        if iszero(radius)
            push!(poses, Pose2(0, 0, pattern.φ0))
            continue
        end
        append!(poses,
            placements(
                Ring(pattern.nφ; r = radius, φ0 = pattern.φ0, span = pattern.span),
                item,
                nothing
            ))
    end
    return poses
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
    return Pose2[Pose2(x0 + ix * pattern.dx, y0 + iy * pattern.dy, 0)
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

    poses = Pose2[Pose2(0, 0, pattern.φ0)]
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
        append!(poses, placements(
            Ring(count; r = ring.r, φ0 = ring.φ0, span = ring.span),
            item,
            compact
        ))
        course += 1
        radius = 2course * extent
    end
    return poses
end

placements(pattern::Fill, item, ::Nothing) = _fill_placements(pattern, item, nothing)
capacity(pattern::Fill, item, ::Nothing) = length(_fill_placements(pattern, item, nothing))

function placements(poses::AbstractVector{<:Pose2}, item, ::Nothing)
    isempty(poses) && throw(ArgumentError("explicit placements cannot be empty"))
    return copy(poses)
end
