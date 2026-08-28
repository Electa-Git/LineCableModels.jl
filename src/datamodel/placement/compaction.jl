"""
$(TYPEDEF)

Scale nominal member-centre offsets by a measured diameter ratio.

$(TYPEDFIELDS)
"""
struct DiameterFactor{T <: Real}
    "Compacted diameter divided by nominal diameter \\[dimensionless\\]."
    k::T
    function DiameterFactor{T}(k::T) where {T <: Real}
        isfinite(k) && zero(k) < k <= one(k) ||
            throw(DomainError(k, "diameter factor must lie in (0, 1]"))
        return new{T}(k)
    end
end
DiameterFactor(k::Real) = DiameterFactor{typeof(float(k))}(float(k))

function placements(pattern, item, factor::DiameterFactor)
    nominal_poses = placements(pattern, item, nothing)
    return Pose2[
        Pose2(factor.k * pose.x, factor.k * pose.y, pose.φ)
        for pose in nominal_poses
    ]
end
