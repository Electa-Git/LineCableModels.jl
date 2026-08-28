"""
$(TYPEDEF)

Place `n` members on one circular locus.

$(TYPEDFIELDS)
"""
struct Ring{T <: Real}
    "Number of members \\[dimensionless\\]."
    n::Int
    "Radius of the member-centre locus \\[m\\]."
    r::T
    "Starting angle \\[rad\\]."
    φ0::T
    "Counter-clockwise angular span \\[rad\\]."
    span::T

    function Ring{T}(n::Int, r::T, φ0::T, span::T) where {T <: Real}
        n > 0 || throw(ArgumentError("ring cardinality must be positive"))
        isfinite(r) && r >= zero(r) ||
            throw(DomainError(r, "ring radius must be nonnegative and finite"))
        isfinite(φ0) || throw(DomainError(φ0, "ring start angle must be finite"))
        isfinite(span) && zero(span) < span <= oftype(span, 2π) ||
            throw(DomainError(span, "ring span must lie in (0, 2π]"))
        return new{T}(n, r, φ0, span)
    end
end

function Ring(n::Integer; r, φ0 = 0, span = 2π)
    values = map(float, promote(r, φ0, span))
    return Ring{typeof(first(values))}(Int(n), values...)
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

function Polar(; nr::Integer, nφ::Integer, r0 = 0, dr, φ0 = 0, span = 2π)
    values = map(float, promote(r0, dr, φ0, span))
    return Polar{typeof(first(values))}(Int(nr), Int(nφ), values...)
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

function Lattice(; nx::Integer, ny::Integer, dx, dy)
    values = map(float, promote(dx, dy))
    return Lattice{typeof(first(values))}(Int(nx), Int(ny), values...)
end

"Return local member poses prescribed by a placement pattern."
function placements end

placements(::Nothing, item, ::Nothing) = Pose2[Pose2(0, 0, 0)]

function placements(pattern::Ring, item, ::Nothing)
    step = pattern.n == 1 ? zero(pattern.span) :
           isapprox(pattern.span, oftype(pattern.span, 2π)) ?
           pattern.span / pattern.n : pattern.span / (pattern.n - 1)
    return Pose2[
        Pose2(
            pattern.r * cos(pattern.φ0 + index * step),
            pattern.r * sin(pattern.φ0 + index * step),
            pattern.φ0 + index * step
        )
        for index in 0:(pattern.n - 1)
    ]
end

function placements(pattern::Polar, item, ::Nothing)
    poses = Pose2[]
    for radial_index in 0:(pattern.nr - 1)
        radius = pattern.r0 + radial_index * pattern.dr
        if iszero(radius)
            push!(poses, Pose2(0, 0, pattern.φ0))
            continue
        end
        append!(poses, placements(
            Ring(pattern.nφ; r = radius, φ0 = pattern.φ0, span = pattern.span),
            item,
            nothing
        ))
    end
    return poses
end

function placements(pattern::Lattice, item, ::Nothing)
    x0 = -(pattern.nx - 1) * pattern.dx / 2
    y0 = -(pattern.ny - 1) * pattern.dy / 2
    return Pose2[
        Pose2(x0 + ix * pattern.dx, y0 + iy * pattern.dy, 0)
        for iy in 0:(pattern.ny - 1) for ix in 0:(pattern.nx - 1)
    ]
end

function placements(poses::AbstractVector{<:Pose2}, item, ::Nothing)
    isempty(poses) && throw(ArgumentError("explicit placements cannot be empty"))
    return copy(poses)
end
