"Minimum surface gap between cable exteriors, and from a cable exterior to the air-earth interface \\[m\\]."
const _CABLE_CLEARANCE = 1.0e-6

"""
$(TYPEDEF)

Cable exteriors that overlap, or a cable exterior that crosses the air-earth interface. A
line-system build throws it. `MonteCarlo` discards a draw that throws it and draws again.

$(TYPEDFIELDS)
"""
struct OverlapError <: Exception
    "Indices of the cables: two for a pair, one at the air-earth interface."
    cables::Tuple{Vararg{Int}}
    "Nominal surface gap \\[m\\]. A negative gap is an overlap."
    gap::Float64
    "Explanation of the failure."
    message::String
end

function Base.showerror(io::IO, error::OverlapError)
    print(io, "OverlapError: ", error.message, " (", length(error.cables) == 1 ? "cable " :
        "cables ", join(error.cables, ", "), "; surface gap ", error.gap, " m)")
end

# The exterior circle of a placed cable. A circular boundary keeps its eccentric local @at
# coordinates in the center, not inside a norm about the declaration origin. Other
# geometric boundaries use their containing circle about the pose.
function _clearance_exterior(shape::Disk, design, pose)
    return (radius = shape.r, center = pose * shape.at)
end
_clearance_exterior(shape, design, pose) = (radius = outer_radius(design), center = pose)

"""
Floating-point allowance for subtraction of placement coordinates \\[m\\].
"""
function _clearance_roundoff(radii, poses)
    scale = maximum(abs ∘ nominal, radii)
    for pose in poses
        scale = max(scale, abs(nominal(pose.x)), abs(nominal(pose.y)))
    end
    return 8eps(float(scale))
end

"""
$(TYPEDSIGNATURES)

Raise every exterior gap below 1 μm to 1 μm, without changing cable dimensions or
rotations. A pair gap from zero, within the roundoff allowance, up to 1 μm is raised by
expanding the relative placements about their centroid, which retains shared translations
and formation symmetry. With `interface`, each half-space group is then translated away
from the air-earth interface until its smallest interface gap is 1 μm. A gap below the
negative roundoff allowance is an overlap and throws [`OverlapError`](@ref). Gaps use
nominal values: uncertainty never moves the nominal geometry.

Circular boundaries use their composed centers and physical radii. Other geometric
boundaries use their containing circles.

# Arguments

- `designs`: completed cable designs in placement order.
- `poses`: cable translations \\[m\\] and rotations \\[rad\\] in the system frame.

# Keywords

- `pairs=true`: resolve the gaps between cables.
- `interface=false`: resolve the gaps to the air-earth interface.
- `adjust=true`: resolve placements. `false` only checks that every gap is at least 1 μm.

# Returns

- Resolved poses and the largest translation magnitude \\[m\\]. The input poses are
  returned unchanged when nothing moves.
"""
function clearance_geometry(designs, poses;
        pairs::Bool = true, interface::Bool = false, adjust::Bool = true)
    length(poses) == length(designs) || throw(DimensionMismatch(
        "clearance requires one pose per cable"))
    exteriors = [_clearance_exterior(boundary(design.geometry), design, pose)
                 for (design, pose) in zip(designs, poses)]
    radii = [exterior.radius for exterior in exteriors]
    centers = [exterior.center for exterior in exteriors]
    T = eltype(eltype(poses))
    roundoff = _clearance_roundoff(radii, centers)
    minimum_gap = oftype(nominal(zero(T)), _CABLE_CLEARANCE)
    minimum_gap < _CABLE_CLEARANCE && (minimum_gap = nextfloat(minimum_gap))
    # The expansion of the relative placements about their centroid that raises every pair
    # gap to `minimum_gap`. With `adjust=false`, every pair is checked instead.
    factor = one(T)
    for i in eachindex(centers), j in 1:(pairs ? i - 1 : 0)
        distance = hypot(centers[i].x - centers[j].x, centers[i].y - centers[j].y)
        gap = distance - radii[i] - radii[j]
        nominal(gap) < -roundoff && throw(OverlapError((j, i), nominal(gap),
            "cable exteriors overlap"))
        if !adjust
            nominal(gap) > 0 && nominal(gap) + roundoff >= minimum_gap ||
                throw(OverlapError((j, i), nominal(gap), "cable exteriors lie closer than 1 μm"))
        elseif nominal(gap) + roundoff < minimum_gap
            required_scale = (radii[i] + radii[j] + minimum_gap + 2roundoff) / distance
            nominal(required_scale) > nominal(factor) && (factor = required_scale)
        end
    end
    expanded = nominal(factor) > 1 ? _expand(centers, factor) : centers
    below, above = interface ?
                   _interface_shift(expanded, centers, radii, minimum_gap, roundoff, adjust, T) :
                   (zero(T), zero(T))
    shifted = nominal(below) > 0 || nominal(above) > 0
    !shifted && expanded === centers && return poses, zero(nominal(zero(T)))
    result = shifted ? [Pose2(expanded[i].x,
                              expanded[i].y + (centers[i].y < 0 ? -below : above),
                              expanded[i].φ) for i in eachindex(expanded)] : expanded
    displacement = maximum(eachindex(result)) do i
        hypot(nominal(result[i].x - centers[i].x), nominal(result[i].y - centers[i].y))
    end
    resolved = [Pose2(poses[i].x + (result[i].x - centers[i].x),
                    poses[i].y + (result[i].y - centers[i].y), poses[i].φ)
                for i in eachindex(poses)]
    return resolved, displacement
end

# Expand the placements by `factor` about their centroid.
function _expand(placements, factor)
    n = length(placements)
    xcenter = sum(pose -> pose.x, placements) / n
    ycenter = sum(pose -> pose.y, placements) / n
    return [Pose2(xcenter + factor * (pose.x - xcenter),
                ycenter + factor * (pose.y - ycenter), pose.φ) for pose in placements]
end

# The translations, away from the air-earth interface, of the cables below it and of those
# above it, that raise each group's smallest interface gap to `minimum_gap`. A common
# translation leaves distances unchanged within a group and increases the separation from
# the other group. `declared` gives each cable's half-space: its center before any
# expansion. With `adjust=false`, every gap is checked instead.
function _interface_shift(placements, declared, radii, minimum_gap, roundoff, adjust::Bool,
        ::Type{T}) where {T}
    below = zero(T)
    above = zero(T)
    for i in eachindex(placements)
        side = sign(nominal(declared[i].y))
        iszero(side) && throw(OverlapError((i,), -nominal(radii[i]),
            "a cable center lies on the air-earth interface"))
        declared_gap = side * declared[i].y - radii[i]
        nominal(declared_gap) < -roundoff && throw(OverlapError((i,), nominal(declared_gap),
            "cable exterior crosses the air-earth interface"))
        gap = side * placements[i].y - radii[i]
        if !adjust
            nominal(gap) > 0 && nominal(gap) + roundoff >= minimum_gap || throw(OverlapError(
                (i,), nominal(gap), "cable exterior lies closer than 1 μm to the air-earth interface"))
        elseif nominal(gap) + roundoff < minimum_gap
            required_shift = minimum_gap + 2roundoff - gap
            if side < 0
                nominal(required_shift) > nominal(below) && (below = required_shift)
            else
                nominal(required_shift) > nominal(above) && (above = required_shift)
            end
        end
    end
    return below, above
end
