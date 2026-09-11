"Minimum engineering separation of cable exteriors \\[m\\]."
const _CABLE_CLEARANCE = 1.0e-6

# A scoped value, not process-global mutable state: nested user builders can
# participate without acquiring an extra public argument. Each sampling task
# owns its records and counters.
const _CLEARANCE_CONTEXT = Base.ScopedValues.ScopedValue{Any}(nothing)

function _clearance_exterior(design, pose)
    shape = boundary(design.geometry)
    if shape isa Disk
        # Keep eccentric/uncertain local @at coordinates in the centre, not
        # inside a norm about the declaration origin (undefined derivative at 0).
        reserve = uncertainty(outer_radius(design))
        if !isfinite(reserve)
            # At the origin, norm(x,y) is not differentiable. Its RMS offset
            # bounds the standard deviation; retain that conservative reserve
            # without manufacturing an independent Measurement.
            reserve = uncertainty(shape.r) + hypot(uncertainty(shape.at.x), uncertainty(shape.at.y))
        end
        reserve = max(reserve, uncertainty(shape.r))
        return (radius = shape.r, centre = pose * shape.at, reserve = reserve)
    end
    radius = outer_radius(design)
    return (radius = radius, centre = pose, reserve = uncertainty(radius))
end

"Floating-point allowance for subtraction of placement coordinates \\[m\\]."
function _clearance_roundoff(radii, poses)
    scale = maximum(abs ∘ nominal, radii)
    for pose in poses
        scale = max(scale, abs(nominal(pose.x)), abs(nominal(pose.y)))
    end
    return 8eps(float(scale))
end

"""
$(TYPEDSIGNATURES)

Resolve exterior-circle clearance by expanding relative cable placements,
without changing cable dimensions or rotations. Shared translations and
formation symmetry are retained. If `interface` is true, each half-space
group is translated away from the air-earth interface when necessary.

The required pair gap is:

```math
c_{ij} = 10^{-6} + \\max\\left(u_{R,\\max}, u(d_{ij}-r_i-r_j)\\right)
```

where all lengths and standard uncertainties are in \\[m\\]. Diagonal entries
retain the corresponding interface clearance. Supplied `required` values
retain the original uncertainty budget after stochastic sampling.

Circular boundaries use their composed centres and physical radii; other
boundaries use their containing circles. A nominally zero uncertain circular
offset uses its RMS displacement as a conservative radius-uncertainty reserve.

Nominal overlaps are rejected; `sampling=true` permits construction of a
feasible realization from an already checked declaration. With `adjust=false`,
only validate the supplied geometry and clearance requirements.

# Arguments

- `designs`: Completed cable designs in placement order.
- `poses`: Cable translations \\[m\\] and rotations \\[rad\\] in the system frame.

# Keywords

- `required=nothing`: Retained pairwise and interface clearance matrix \\[m\\].
- `reference=poses`: Declared poses used to distinguish nominal overlaps from sampled ones.
- `interface=false`: Enforce clearance from the air-earth interface.
- `sampling=false`: Permit correction of overlaps produced by a sampled declaration.
- `adjust=true`: Resolve placements; `false` performs validation only.
- `reference_centres=nothing`: Composed exterior-centre poses retained before sampling.

# Returns

- Resolved poses, retained clearance matrix \\[m\\], and maximum translation
  magnitude \\[m\\]. Input arrays are not modified.
"""
function clearance_geometry(designs, poses;
        required = nothing, reference = poses, interface::Bool = false,
        sampling::Bool = false, adjust::Bool = true, reference_centres = nothing)
    exteriors = _clearance_exterior.(designs, poses)
    radii = getproperty.(exteriors, :radius)
    centres = getproperty.(exteriors, :centre)
    reference_centres === nothing &&
        (reference_centres = getproperty.(_clearance_exterior.(designs, reference), :centre))
    n = length(radii)
    length(poses) == length(reference) == n || throw(DimensionMismatch(
        "clearance requires one pose per cable"))
    T = eltype(eltype(poses))
    roundoff = _clearance_roundoff(radii, centres)
    if required === nothing
        minimum_gap = oftype(nominal(zero(T)), _CABLE_CLEARANCE)
        minimum_gap < _CABLE_CLEARANCE && (minimum_gap = nextfloat(minimum_gap))
        radius_uncertainty = maximum(exterior -> exterior.reserve, exteriors)
        required = zeros(T, n, n)
        for i in 1:n
            required[i, i] = minimum_gap + max(radius_uncertainty,
                uncertainty(abs(centres[i].y) - radii[i]))
            for j in 1:(i - 1)
                distance = hypot(centres[i].x - centres[j].x, centres[i].y - centres[j].y)
                clearance = minimum_gap + max(radius_uncertainty,
                    uncertainty(distance - radii[i] - radii[j]))
                required[i, j] = required[j, i] = clearance
            end
        end
    end
    size(required) == (n, n) || throw(DimensionMismatch(
        "clearance matrix must match the cable count"))
    all(value -> isfinite(value) && nominal(value) >= oftype(nominal(value), _CABLE_CLEARANCE) &&
                 iszero(uncertainty(value)), required) || throw(ArgumentError(
        "retained clearances must be finite deterministic lengths of at least 1 μm"))
    required == transpose(required) || throw(ArgumentError(
        "retained cable clearances must be symmetric"))

    result = copy(centres)
    if sampling && any(((i, j),) ->
            iszero(nominal(result[i].x - result[j].x)) &&
            iszero(nominal(result[i].y - result[j].y)),
            ((i, j) for i in 1:n for j in 1:(i - 1)))
        # A coincident draw has no direction. Reuse the declared relative
        # layout at the sampled centre, then apply the same clearance rule.
        cx = sum(p -> p.x, result) / n
        cy = sum(p -> p.y, result) / n
        rx = sum(p -> nominal(p.x), reference_centres) / n
        ry = sum(p -> nominal(p.y), reference_centres) / n
        result = [Pose2(cx + nominal(reference_centres[i].x) - rx,
                       cy + nominal(reference_centres[i].y) - ry, centres[i].φ) for i in 1:n]
    end
    factor = one(T)
    for i in 1:n, j in 1:(i - 1)
        distance = hypot(result[i].x - result[j].x, result[i].y - result[j].y)
        gap = distance - radii[i] - radii[j]
        if !adjust
            nominal(gap) > 0 && nominal(gap) + roundoff >= nominal(required[i, j]) ||
                throw(DomainError((i, j), "cable exterior clearance is below its retained minimum"))
            continue
        end
        !sampling && nominal(gap) < -roundoff && throw(DomainError((i, j),
            "cable cross-sections overlap; only touching or insufficient positive clearance can be adjusted"))
        nominal(distance) > 0 || throw(DomainError((i, j),
            "coincident cable centres do not define a separation direction"))
        if nominal(gap) + roundoff < nominal(required[i, j])
            candidate = (radii[i] + radii[j] + required[i, j] + 2roundoff) / distance
            nominal(candidate) > nominal(factor) && (factor = candidate)
        end
    end
    if nominal(factor) > 1
        xcentre = sum(pose -> pose.x, result) / n
        ycentre = sum(pose -> pose.y, result) / n
        result = [Pose2(xcentre + factor * (pose.x - xcentre),
                       ycentre + factor * (pose.y - ycentre), pose.φ) for pose in result]
    end

    if interface
        # A common translation within each half-space preserves all distances
        # within that group and increases separation from the other group.
        for side in (-1, 1)
            indices = findall(pose -> sign(nominal(pose.y)) == side, reference_centres)
            shift = zero(T)
            for i in indices
                original_gap = side * reference_centres[i].y - radii[i]
                gap = side * result[i].y - radii[i]
                adjust && !sampling && nominal(gap) + roundoff < nominal(required[i, i]) &&
                    nominal(original_gap) < -roundoff && throw(DomainError(i,
                    "cable cross-section crosses the air-earth interface"))
                if !adjust
                    nominal(gap) > 0 && nominal(gap) + roundoff >= nominal(required[i, i]) ||
                        throw(DomainError(i, "cable clearance from the air-earth interface is below its retained minimum"))
                elseif nominal(gap) + roundoff < nominal(required[i, i])
                    candidate = required[i, i] + 2roundoff - gap
                    nominal(candidate) > nominal(shift) && (shift = candidate)
                end
            end
            if nominal(shift) > 0
                for i in indices
                    pose = result[i]
                    result[i] = Pose2(pose.x, pose.y + side * shift, pose.φ)
                end
            end
        end
        any(pose -> iszero(nominal(pose.y)), reference_centres) && throw(DomainError(
            reference, "a declared cable centre cannot lie on the air-earth interface"))
        if adjust
            # Extreme draws can exchange half-spaces before projection. If
            # restoring their declared sides shortened a cross-interface
            # distance, expansion about y=0 preserves the restored sides.
            factor = one(T)
            for i in 1:n, j in 1:(i - 1)
                distance = hypot(result[i].x - result[j].x, result[i].y - result[j].y)
                limit = radii[i] + radii[j] + required[i, j]
                if nominal(distance) + roundoff < nominal(limit)
                    candidate = (limit + 2roundoff) / distance
                    nominal(candidate) > nominal(factor) && (factor = candidate)
                end
            end
            if nominal(factor) > 1
                xcentre = sum(p -> p.x, result) / n
                result = [Pose2(xcentre + factor * (pose.x - xcentre),
                    factor * pose.y, pose.φ) for pose in result]
            end
        end
    end
    displacement = maximum(eachindex(result)) do i
        hypot(nominal(result[i].x - centres[i].x), nominal(result[i].y - centres[i].y))
    end
    resolved = iszero(displacement) ? copy(poses) :
        [Pose2(poses[i].x + (result[i].x - centres[i].x),
               poses[i].y + (result[i].y - centres[i].y), poses[i].φ) for i in 1:n]
    if adjust && displacement > 0
        clearance_geometry(designs, resolved; required, reference, reference_centres,
            interface, sampling, adjust = false)
    end
    return resolved, Matrix{T}(required), displacement
end

function _record_clearance_adjustment(system_id, displacement)
    displacement > 0 || return nothing
    context = _CLEARANCE_CONTEXT[]
    if context === nothing
        @warn "Cable placements adjusted to preserve exterior clearance" system_id max_displacement_m=displacement
    else
        context.adjustments[] += 1
        context.max_displacement[] = max(context.max_displacement[], displacement)
    end
    return nothing
end

"""
$(TYPEDSIGNATURES)

Materialize one unresolved line-system point to retain its uncertainty-bearing
clearance requirements before sampling. Return task-scoped construction records
and adjustment counters; do not consume the sampling RNG.
"""
function prepare_clearance(point)
    context = (records = Any[], references = IdDict{Any, Any}(),
        cursor = Ref(0), sampling = Ref(false),
        adjustments = Ref(0), max_displacement = Ref(0.0))
    Base.ScopedValues.with(_CLEARANCE_CONTEXT => context) do
        if Base.get_extension(parentmodule(@__MODULE__), :LineCableModelsMeasurementsExt) === nothing
            # Monte Carlo remains usable without Measurements. The geometric
            # constraint still applies to every draw; propagated uncertainty
            # reserves are available when Measurements is loaded.
            realize(point, realize_arguments(Random.Xoshiro(0), point,
                (_rng, mean, _sigma) -> mean))
        else
            materialize(point)
        end
    end
    context.sampling[] = true
    context.adjustments[] = 0
    context.max_displacement[] = 0.0
    return context
end

"""
$(TYPEDSIGNATURES)

Execute `f` while reconstructing one realization with prepared clearance
requirements. Return the result of `f` and restore the previous scope on exit.
"""
function with_clearance(f, context)
    context === nothing && return f()
    return Base.ScopedValues.with(_CLEARANCE_CONTEXT => context) do
        context.cursor[] = 0
        empty!(context.references)
        value = f()
        context.cursor[] == length(context.records) || throw(ArgumentError(
            "sampled construction produced fewer cable systems than its declaration"))
        return value
    end
end

"""
$(TYPEDSIGNATURES)

Return the number of adjusted builds and maximum translation \\[m\\] recorded
in a sampling context. A missing context returns zero counts and displacement.
"""
function clearance_summary(context)
    context === nothing && return (adjustments = 0, max_displacement_m = 0.0)
    return (adjustments = context.adjustments[],
        max_displacement_m = context.max_displacement[])
end

"""
$(TYPEDSIGNATURES)

Emit one aggregate warning if any sampled placements required adjustment.
Return `nothing`; no warning is emitted for an unchanged geometry.
"""
function warn_clearance_summary(context)
    summary = clearance_summary(context)
    summary.adjustments > 0 &&
        @warn "Sampled cable placements adjusted to preserve exterior clearance" summary...
    return nothing
end

"""
$(TYPEDSIGNATURES)

Sample and build one line-system point with its retained clearance requirements.
Reuse an enclosing sampling context or prepare one for a standalone draw.
Return the completed system or problem.
"""
function realize_clearance(rng, point, distribution)
    _CLEARANCE_CONTEXT[] === nothing ||
        return realize(point, realize_arguments(rng, point, distribution))
    context = prepare_clearance(point)
    try
        return with_clearance(context) do
            realize(point, realize_arguments(rng, point, distribution))
        end
    finally
        warn_clearance_summary(context)
    end
end
