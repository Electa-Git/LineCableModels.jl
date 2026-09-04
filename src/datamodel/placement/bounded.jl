"Record the resolved boundary shared by every member of one bounded formation."
struct BoundedPlacement{B <: AbstractShape}
    boundary::B
end

function resolve(at::Pose2, region::PlacedRegion)
    patterns = map(region.placement.patterns) do entry
        pattern = entry.pattern isa BoundedPlacement ?
                  BoundedPlacement(resolve(at, entry.pattern.boundary)) : entry.pattern
        (pattern = pattern, member = entry.member, pose = entry.pose)
    end
    return PlacedRegion(
        region.source,
        resolve(at, region.primitive),
        region.terminal,
        (patterns = patterns,),
        region.paths
    )
end

"""
$(TYPEDEF)

Store the exact annular deformation of one rectangular strip.

`BentStrip` is resolved geometry, not a modelling declaration. The source
remains a [`Rectangle`](@ref); the strip preserves its source area while its
radial faces follow the containing circular course.

For angular coverage ``\\Delta\\phi`` and radial limits ``r_i`` and ``r_o``,

```math
A = \\frac{\\Delta\\phi}{2}\\left(r_o^2-r_i^2\\right) = w h,
```

where ``w h`` is the area of the source rectangle.

$(TYPEDFIELDS)
"""
struct BentStrip{T <: Real, P <: Pose2{T}} <: AbstractShape{T}
    "Inner course radius \\[m\\]."
    ri::T
    "Outer course radius \\[m\\]."
    ro::T
    "Angular coverage \\[rad\\]."
    span::T
    "Resolved pose of the strip symmetry axis."
    at::P

    function BentStrip{T, P}(ri::T, ro::T, span::T, at::P) where {
            T <: Real, P <: Pose2{T}
    }
        isfinite(ri) && ri >= zero(ri) || throw(DomainError(
            ri, "bent-strip inner radius must be nonnegative and finite"
        ))
        isfinite(ro) && ro > ri || throw(DomainError(
            ro, "bent-strip outer radius must exceed its inner radius"
        ))
        isfinite(span) && zero(span) < span <= oftype(span, 2pi) ||
            throw(DomainError(span, "bent-strip span must lie in (0, 2pi]"))
        return new{T, P}(ri, ro, span, at)
    end
end

function BentStrip(ri::Real, ro::Real, span::Real, at::Pose2 = Pose2(0, 0, 0))
    values = map(float, promote(ri, ro, span, at.x, at.y, at.φ))
    inner, outer, angle, x, y, rotation = values
    pose = Pose2(x, y, rotation)
    return BentStrip{typeof(inner), typeof(pose)}(inner, outer, angle, pose)
end

area(shape::BentStrip) = shape.span * (shape.ro^2 - shape.ri^2) / 2
perimeter(shape::BentStrip) =
    shape.span * (shape.ro + shape.ri) + 2(shape.ro - shape.ri)
r_in(shape::BentStrip) = shape.ri
r_ex(shape::BentStrip) = shape.ro
thickness(shape::BentStrip) = shape.ro - shape.ri
boundary(shape::BentStrip) = shape

function centroid(shape::BentStrip)
    radius = isapprox(shape.span, oftype(shape.span, 2pi)) ? zero(shape.span) :
             4sin(shape.span / 2) * (shape.ro^3 - shape.ri^3) /
             (3shape.span * (shape.ro^2 - shape.ri^2))
    return (
        shape.at.x + radius * cos(shape.at.φ),
        shape.at.y + radius * sin(shape.at.φ)
    )
end

function support(shape::BentStrip, φ::Real)
    direction = φ - shape.at.φ
    half = shape.span / 2
    candidates = (
        shape.ro * cos(direction - half),
        shape.ro * cos(direction + half),
        shape.ri * cos(direction - half),
        shape.ri * cos(direction + half)
    )
    period = oftype(direction, 2pi)
    wrapped = mod(direction + pi, period) - pi
    radial = abs(wrapped) <= half ? shape.ro : maximum(candidates)
    return shape.at.x * cos(φ) + shape.at.y * sin(φ) + radial
end

support(shape::BentStrip) = hypot(shape.at.x, shape.at.y) + shape.ro

function resolve(at::Pose2, shape::BentStrip)
    return BentStrip(shape.ri, shape.ro, shape.span, at * shape.at)
end

function tessellate(shape::BentStrip; points_per_arc::Integer = 32)
    points_per_arc >= 2 || throw(ArgumentError(
        "points_per_arc must be at least two"
    ))
    outer = range(-shape.span / 2, shape.span / 2; length = Int(points_per_arc))
    inner = reverse(outer)
    local_points = [
        (shape.ro * cos(angle), shape.ro * sin(angle)) for angle in outer
    ]
    append!(local_points,
        [(shape.ri * cos(angle), shape.ri * sin(angle)) for angle in inner])
    cosine = cos(shape.at.φ)
    sine = sin(shape.at.φ)
    return [
        (
            shape.at.x + cosine * point[1] - sine * point[2],
            shape.at.y + sine * point[1] + cosine * point[2]
        ) for point in local_points
    ]
end

"Return a counter-clockwise polygonal approximation of a resolved boundary."
function boundary_polygon(shape::Disk; points::Integer = 512)
    points >= 16 || throw(ArgumentError(
        "a disk boundary polygon requires at least 16 points"
    ))
    centre = (shape.at.x, shape.at.y)
    angles = range(zero(shape.r), oftype(shape.r, 2π); length = Int(points) + 1)
    polygon = [(centre[1] + shape.r * cos(angle), centre[2] + shape.r * sin(angle))
               for angle in Iterators.drop(angles, 1)]
    return polygon
end

function boundary_polygon(shape::SectorShape; points::Integer = 64)
    polygon = collect(tessellate(shape; points_per_arc = points))
    length(polygon) > 1 &&
        isapprox(first(polygon)[1], last(polygon)[1]) &&
        isapprox(first(polygon)[2], last(polygon)[2]) && pop!(polygon)
    return signed_polygon_area(polygon) < 0 ? reverse(polygon) : polygon
end

function signed_polygon_area(points)
    return sum(eachindex(points); init = zero(first(points)[1])) do index
        next = mod1(index + 1, length(points))
        points[index][1] * points[next][2] - points[next][1] * points[index][2]
    end / 2
end

function polygon_centroid(points)
    signed_area = signed_polygon_area(points)
    abs(signed_area) > eps(float(max(abs(signed_area), one(signed_area)))) ||
        throw(ArgumentError("a compaction polygon must have positive area"))
    xmoment = zero(signed_area)
    ymoment = zero(signed_area)
    for index in eachindex(points)
        next = mod1(index + 1, length(points))
        cross = points[index][1] * points[next][2] -
                points[next][1] * points[index][2]
        xmoment += (points[index][1] + points[next][1]) * cross
        ymoment += (points[index][2] + points[next][2]) * cross
    end
    return (xmoment / (6signed_area), ymoment / (6signed_area))
end

function normalize_polygon_area(points, target_area::Real, centre)
    polygon = signed_polygon_area(points) < 0 ? reverse(points) : points
    current_area = signed_polygon_area(polygon)
    current_area > zero(current_area) || throw(ArgumentError(
        "a bounded formation requires a counter-clockwise positive-area boundary"
    ))
    factor = sqrt(target_area / current_area)
    return [
        (
            centre[1] + factor * (point[1] - centre[1]),
            centre[2] + factor * (point[2] - centre[2])
        ) for point in polygon
    ]
end

function clip_halfplane(points, normal, offset)
    isempty(points) && return points
    scale = max(abs(offset), one(offset))
    tolerance = 128eps(float(scale)) * scale
    output = eltype(points)[]
    previous = last(points)
    previous_value = normal[1] * previous[1] + normal[2] * previous[2] - offset
    previous_inside = previous_value <= tolerance
    for current in points
        current_value = normal[1] * current[1] + normal[2] * current[2] - offset
        current_inside = current_value <= tolerance
        if current_inside != previous_inside
            denominator = previous_value - current_value
            fraction = iszero(denominator) ? zero(denominator) :
                       previous_value / denominator
            push!(output,
                (
                    previous[1] + fraction * (current[1] - previous[1]),
                    previous[2] + fraction * (current[2] - previous[2])
                ))
        end
        current_inside && push!(output, current)
        previous = current
        previous_value = current_value
        previous_inside = current_inside
    end
    return output
end

function power_cell(boundary, sites, weights, index::Int)
    cell = copy(boundary)
    site = sites[index]
    for neighbour in eachindex(sites)
        neighbour == index && continue
        other = sites[neighbour]
        dx = other[1] - site[1]
        dy = other[2] - site[2]
        iszero(dx) && iszero(dy) && throw(ArgumentError(
            "bounded compaction requires distinct initial sites"
        ))
        offset = other[1]^2 + other[2]^2 - site[1]^2 - site[2]^2 +
                 weights[index] - weights[neighbour]
        cell = clip_halfplane(cell, (2dx, 2dy), offset)
        isempty(cell) && break
    end
    return cell
end

power_cells(boundary, sites, weights) =
    [power_cell(boundary, sites, weights, index) for index in eachindex(sites)]

function shared_face_length(cell, left, right, left_weight, right_weight, scale)
    dx = right[1] - left[1]
    dy = right[2] - left[2]
    distance = hypot(dx, dy)
    iszero(distance) && return zero(distance)
    offset = right[1]^2 + right[2]^2 - left[1]^2 - left[2]^2 +
             left_weight - right_weight
    tolerance = 2048eps(float(scale)) * max(scale, one(scale))
    projections = typeof(float(scale))[]
    for point in cell
        residual = abs(2dx * point[1] + 2dy * point[2] - offset) / (2distance)
        residual <= tolerance || continue
        push!(projections, (-dy * point[1] + dx * point[2]) / distance)
    end
    length(projections) >= 2 || return zero(distance)
    return maximum(projections) - minimum(projections)
end

function weight_hessian(cells, sites, weights, scale)
    count = length(sites)
    matrix = zeros(typeof(float(scale)), count, count)
    for left in 1:(count - 1), right in (left + 1):count
        distance = hypot(
            sites[right][1] - sites[left][1],
            sites[right][2] - sites[left][2]
        )
        face = shared_face_length(
            cells[left], sites[left], sites[right],
            weights[left], weights[right], scale
        )
        iszero(face) && continue
        derivative = face / (2distance)
        matrix[left, left] += derivative
        matrix[right, right] += derivative
        matrix[left, right] -= derivative
        matrix[right, left] -= derivative
    end
    return matrix
end

function balance_power_cells(boundary, sites, targets; maxiter::Integer = 50)
    count = length(sites)
    count == length(targets) || throw(DimensionMismatch(
        "power-cell sites and target areas must have equal lengths"
    ))
    count == 1 && return ([copy(boundary)], zeros(typeof(first(targets)), 1))
    scale = maximum(point -> hypot(point...), boundary)
    weights = zeros(typeof(float(first(targets))), count)
    best_error = oftype(first(targets), Inf)
    for _ in 1:Int(maxiter)
        cells = power_cells(boundary, sites, weights)
        any(cell -> length(cell) < 3, cells) && throw(ArgumentError(
            "bounded compaction produced an empty power cell"
        ))
        areas = signed_polygon_area.(cells)
        residual = targets .- areas
        relative_error = maximum(abs.(residual) ./ targets)
        relative_error <= 2.0e-6 && return (cells, weights)
        if relative_error < best_error
            best_error = relative_error
        end

        hessian = weight_hessian(cells, sites, weights, scale)
        reduced = hessian[1:(end - 1), 1:(end - 1)]
        regularization = eps(float(scale^2)) / max(scale^2, eps(float(scale^2)))
        for index in axes(reduced, 1)
            reduced[index, index] += regularization
        end
        step = try
            LinearAlgebra.qr(reduced) \ residual[1:(end - 1)]
        catch
            (scale^2 / sum(targets)) .* residual[1:(end - 1)]
        end

        accepted = false
        damping = one(eltype(weights))
        for _ in 1:10
            candidate = copy(weights)
            candidate[1:(end - 1)] .+= damping .* step
            candidate .-= candidate[end]
            candidate_cells = power_cells(boundary, sites, candidate)
            if all(cell -> length(cell) >= 3, candidate_cells)
                candidate_areas = signed_polygon_area.(candidate_cells)
                candidate_error = maximum(abs.((targets .- candidate_areas) ./ targets))
                if candidate_error < relative_error
                    weights = candidate
                    accepted = true
                    break
                end
            end
            damping /= 2
        end
        if !accepted
            weights .+= (scale^2 / sum(targets)) .* residual ./ 4
            weights .-= weights[end]
        end
    end
    throw(ArgumentError(
        "bounded compaction did not converge; maximum relative cell-area error " *
        "was $(best_error), target is 2e-6"
    ))
end

function compact_power_cells(boundary, sites, targets; relaxations::Integer = 1)
    current_sites = copy(sites)
    cells = Vector{eltype(boundary)}[]
    for pass in 1:Int(relaxations)
        cells, _ = balance_power_cells(boundary, current_sites, targets)
        pass == relaxations && break
        current_sites = polygon_centroid.(cells)
    end
    return cells
end

function axial_lattice(shells::Integer, spacing::Real)
    shells >= 0 || throw(ArgumentError("axial lattice shell count must be nonnegative"))
    T = typeof(float(spacing))
    result = NamedTuple{(:site, :course), Tuple{Tuple{T, T}, Int}}[]
    push!(result, (site = (zero(T), zero(T)), course = 0))
    for course in 1:Int(shells)
        coordinates = Tuple{Int, Int}[]
        for q in (-course):course, r in (-course):course
            max(abs(q), abs(r), abs(q + r)) == course || continue
            push!(coordinates, (q, r))
        end
        sort!(coordinates; by = coordinate -> atan(
            sqrt(T(3)) * (coordinate[2] + coordinate[1] / 2),
            T(3) * coordinate[1] / 2
        ))
        for (q, r) in coordinates
            push!(result, (
                site = (
                    T(spacing) * (T(q) + T(r) / 2),
                    T(spacing) * sqrt(T(3)) * T(r) / 2
                ),
                course
            ))
        end
    end
    return result
end

function disk_wire_sites(shape::Disk, centre::Disk, wire::Disk)
    available = shape.r - wire.r
    available >= zero(available) || throw(DomainError(
        wire.r, "one circular wire does not fit inside the disk boundary"
    ))
    spacing = 2wire.r
    shell_limit = ceil(Int, nominal(available / spacing)) + 2
    tolerance = 256 * _geometry_tolerance(max(shape.r, centre.r, wire.r))
    sites = NamedTuple[]
    for candidate in Iterators.drop(axial_lattice(shell_limit, spacing), 1)
        x = shape.at.x + candidate.site[1]
        y = shape.at.y + candidate.site[2]
        radius = hypot(candidate.site...)
        radius + wire.r <= shape.r + tolerance || continue
        radius + tolerance >= centre.r + wire.r || continue
        push!(sites, (site = (x, y), course = candidate.course))
    end
    isempty(sites) && throw(DomainError(
        wire.r, "the disk boundary admits no strand outside its centre wire"
    ))
    return sites
end

function compact_disk_sites(shape::Disk, count::Integer)
    count > 0 || throw(ArgumentError("compacted disk inventory must be positive"))
    count == 1 && return [(site = (shape.at.x, shape.at.y), course = 0)]
    shell = 1
    lattice = axial_lattice(shell, one(shape.r))
    while length(lattice) < count
        shell += 1
        lattice = axial_lattice(shell, one(shape.r))
    end
    selected = lattice[1:Int(count)]
    outer = maximum(candidate -> hypot(candidate.site...), selected)
    scale = iszero(outer) ? zero(shape.r) : 0.82shape.r / outer
    return [(
                site = (
                    shape.at.x + scale * candidate.site[1],
                    shape.at.y + scale * candidate.site[2]
                ),
                candidate.course
            ) for candidate in selected]
end

function point_inside_convex(point, polygon)
    tolerance = 256 * _geometry_tolerance(maximum(value -> hypot(value...), polygon))
    return all(eachindex(polygon)) do index
        next = mod1(index + 1, length(polygon))
        edge = (
            polygon[next][1] - polygon[index][1],
            polygon[next][2] - polygon[index][2]
        )
        offset = (
            point[1] - polygon[index][1],
            point[2] - polygon[index][2]
        )
        edge[1] * offset[2] - edge[2] * offset[1] >= -tolerance
    end
end

function sector_compact_sites(shape::SectorShape, count::Integer)
    count > 0 || throw(ArgumentError("compacted sector inventory must be positive"))
    polygon = boundary_polygon(shape; points = 64)
    xs = first.(polygon)
    ys = last.(polygon)
    resolution = max(12, ceil(Int, sqrt(12count)))
    candidates = Tuple{typeof(shape.at.x), typeof(shape.at.y)}[]
    while length(candidates) < 4count
        empty!(candidates)
        for y in range(minimum(ys), maximum(ys); length = resolution)
            for x in range(minimum(xs), maximum(xs); length = resolution)
                point_inside_convex((x, y), polygon) && push!(candidates, (x, y))
            end
        end
        resolution *= 2
        resolution <= 2048 || throw(ArgumentError(
            "could not seed the compacted sector inventory"
        ))
    end
    centre = centroid(shape)
    first_index = findmin(point -> hypot(
        point[1] - centre[1], point[2] - centre[2]
    ), candidates)[2]
    selected = [candidates[first_index]]
    deleteat!(candidates, first_index)
    distances = [
        (point[1] - selected[1][1])^2 + (point[2] - selected[1][2])^2
        for point in candidates
    ]
    while length(selected) < count
        index = findmax(distances)[2]
        point = candidates[index]
        push!(selected, point)
        deleteat!(candidates, index)
        deleteat!(distances, index)
        for candidate in eachindex(candidates)
            distance = (candidates[candidate][1] - point[1])^2 +
                       (candidates[candidate][2] - point[2])^2
            distances[candidate] = min(distances[candidate], distance)
        end
    end
    return selected
end

function rectangular_strands(shape::Disk, centre::Disk, strand::Rectangle)
    remaining_area = area(shape) - area(centre)
    remaining_area > zero(remaining_area) || throw(DomainError(
        centre.r, "the centre wire leaves no area for rectangular strands"
    ))
    count = floor(Int, nominal(remaining_area / area(strand)) +
                       64eps(float(nominal(remaining_area / area(strand)))))
    count > 0 || throw(DomainError(
        area(strand), "one rectangular strand does not fit outside the centre wire"
    ))
    result = NamedTuple[]
    inner = centre.r
    remaining = count
    course = 0
    while remaining > 0
        course += 1
        mean_radius = inner + strand.h / 2
        nominal_count = max(1, floor(Int, nominal(2pi * mean_radius / strand.w)))
        member_count = min(remaining, nominal_count)
        outer = sqrt(inner^2 + member_count * area(strand) / pi)
        if outer > shape.r
            member_count = min(
                remaining,
                floor(Int, nominal(pi * (shape.r^2 - inner^2) / area(strand)))
            )
            member_count > 0 || break
            outer = sqrt(inner^2 + member_count * area(strand) / pi)
        end
        span = 2pi / member_count
        for member in 1:member_count
            angle = (member - 1) * span
            primitive = BentStrip(
                inner,
                outer,
                span,
                Pose2(shape.at.x, shape.at.y, shape.at.φ + angle)
            )
            push!(result, (
                primitive,
                site = centroid(primitive),
                course,
                member
            ))
        end
        inner = outer
        remaining -= member_count
    end
    return result
end

function clearance(shape::SectorShape, centre)
    cosine = cos(shape.at.φ)
    sine = sin(shape.at.φ)
    dx = centre[1] - shape.at.x
    dy = centre[2] - shape.at.y
    point = (cosine * dx + sine * dy, -sine * dx + cosine * dy)
    distance = oftype(shape.primitive.r_back, Inf)

    for segment in values(shape.contacts.segments)
        first_point, last_point = segment
        edge = (
            last_point[1] - first_point[1],
            last_point[2] - first_point[2]
        )
        edge_length = hypot(edge...)
        outward = (edge[2] / edge_length, -edge[1] / edge_length)
        clearance =
            first_point[1] * outward[1] + first_point[2] * outward[2] -
            point[1] * outward[1] - point[2] * outward[2]
        distance = min(distance, clearance)
    end

    for arc in values(shape.contacts.arcs)
        iszero(arc.radius) && continue
        offset = (point[1] - arc.center[1], point[2] - arc.center[2])
        direction = atan(offset[2], offset[1])
        projection = if _angle_in_arc(direction, arc)
            hypot(offset...)
        else
            max(
                offset[1] * cos(arc.start) + offset[2] * sin(arc.start),
                offset[1] * cos(arc.stop) + offset[2] * sin(arc.stop)
            )
        end
        distance = min(distance, arc.radius - projection)
    end
    return distance
end

function accommodates(shape::SectorShape, centre, radius::Real)
    tolerance = 256 * _geometry_tolerance(
        max(shape.primitive.r_back, radius)
    )
    return _geometry_scalar(clearance(shape, centre) + tolerance - radius) >= 0
end

function sector_point(shape::SectorShape, radius, angle)
    direction = shape.at.φ + angle
    return (
        shape.at.x + radius * cos(direction),
        shape.at.y + radius * sin(direction)
    )
end

function sector_limit(shape::SectorShape, radius, strand_radius)
    accommodates(shape, sector_point(shape, radius, zero(radius)), strand_radius) ||
        return nothing
    lower = zero(shape.primitive.span)
    upper = shape.primitive.span / 2
    accommodates(shape, sector_point(shape, radius, upper), strand_radius) &&
        return upper
    for _ in 1:64
        middle = (lower + upper) / 2
        if accommodates(shape, sector_point(shape, radius, middle), strand_radius)
            lower = middle
        else
            upper = middle
        end
    end
    return lower
end

function shell_clear(left, right, diameter)
    (isempty(left) || isempty(right)) && return true
    scale = max(
        diameter,
        maximum(point -> hypot(point...), left),
        maximum(point -> hypot(point...), right)
    )
    tolerance = 512 * _geometry_tolerance(scale)
    threshold = _geometry_scalar(diameter - tolerance)
    return all(left_point -> all(right_point ->
            _geometry_scalar(hypot(
                left_point[1] - right_point[1],
                left_point[2] - right_point[2]
            )) >= threshold, right), left)
end

function shell_gap(left, right)
    return minimum(left_point -> minimum(right_point -> hypot(
            left_point[1] - right_point[1],
            left_point[2] - right_point[2]
        ), right), left)
end

function shell_options(shape::SectorShape, radius, wire::Disk, shell::Int)
    limit = sector_limit(shape, radius, wire.r)
    limit === nothing && return NamedTuple[]
    ratio = wire.r / radius
    _geometry_scalar(ratio) <= 1 || return NamedTuple[]
    minimum_step = 2asin(ratio)
    width = 2limit
    maximum_count = floor(Int,
        _geometry_scalar(width / minimum_step) + 1 +
        64eps(float(_geometry_scalar(width / minimum_step))))
    maximum_count = max(maximum_count, 1)
    phase_order = isodd(shell) ?
                  (0.0, 1.0, 0.5, 0.25, 0.75, 0.125, 0.875, 0.375, 0.625) :
                  (1.0, 0.0, 0.5, 0.75, 0.25, 0.875, 0.125, 0.625, 0.375)
    options = NamedTuple[]
    for count in maximum_count:-1:1
        if count == 1
            angles = (zero(limit),)
            points = [sector_point(shape, radius, only(angles))]
            push!(options, (radius, angles, points, count, arc = zero(radius)))
            continue
        end
        maximum_step = width / (count - 1)
        _geometry_scalar(maximum_step - minimum_step) >=
            -256 * _geometry_tolerance(maximum_step) || continue
        steps = typeof(minimum_step)[]
        for fraction in phase_order
            step = minimum_step + fraction * (maximum_step - minimum_step)
            any(candidate -> isapprox(
                    _geometry_scalar(candidate), _geometry_scalar(step);
                    rtol = 64eps(float(_geometry_scalar(step))), atol = 0
                ), steps) || push!(steps, step)
        end
        for step in steps
            angles = ntuple(
                index -> (index - (count + 1) / 2) * step,
                count
            )
            points = [sector_point(shape, radius, angle) for angle in angles]
            all(point -> accommodates(shape, point, wire.r), points) || continue
            shell_clear(points, 2wire.r) || continue
            push!(options, (
                radius,
                angles,
                points,
                count,
                arc = radius * (last(angles) - first(angles))
            ))
        end
    end
    return options
end

function shell_clear(points, diameter)
    length(points) < 2 && return true
    scale = max(diameter, maximum(point -> hypot(point...), points))
    tolerance = 512 * _geometry_tolerance(scale)
    threshold = _geometry_scalar(diameter - tolerance)
    return all(
        _geometry_scalar(hypot(
            points[left][1] - points[right][1],
            points[left][2] - points[right][2]
        )) >= threshold
        for left in 1:(length(points) - 1)
        for right in (left + 1):length(points)
    )
end

function phase_layout(shape::SectorShape, wire::Disk, radii, odd_first::Bool)
    states = NamedTuple[]
    for (shell, radius) in enumerate(radii)
        options = shell_options(shape, radius, wire, shell)
        isempty(options) && return nothing
        odd_shell = isodd(shell) == odd_first
        candidates = filter(option -> isodd(option.count) == odd_shell, options)
        isempty(candidates) && return nothing

        if shell == 1
            states = [(
                layout = [candidate],
                score = (
                    candidate.count,
                    zero(_geometry_scalar(candidate.arc)),
                    _geometry_scalar(candidate.arc)
                )
            ) for candidate in candidates]
            continue
        end

        next_states = NamedTuple[]
        for candidate in candidates
            best_state = nothing
            best_score = nothing
            for state in states
                previous = last(state.layout)
                candidate.count >= previous.count || continue
                shell_clear(previous.points, candidate.points, 2wire.r) || continue
                score = (
                    state.score[1] + candidate.count,
                    state.score[2] - abs(_geometry_scalar(
                        shell_gap(previous.points, candidate.points) - 2wire.r
                    )),
                    state.score[3] + _geometry_scalar(candidate.arc)
                )
                if best_state === nothing || score > best_score
                    best_state = state
                    best_score = score
                end
            end
            best_state === nothing && continue
            push!(next_states, (
                layout = [best_state.layout; candidate],
                score = best_score
            ))
        end
        isempty(next_states) && return nothing
        states = next_states
    end
    best = first(states)
    for state in Iterators.drop(states, 1)
        state.score > best.score && (best = state)
    end
    return best.layout
end

function sector_wire_sites(shape::SectorShape, wire::Disk)
    minimum_pitch = sqrt(oftype(wire.r, 3)) * wire.r
    maximum_pitch = 2wire.r
    minimum_radius = shape.primitive.r_base + wire.r
    maximum_radius = shape.primitive.r_back - wire.r
    _geometry_scalar(maximum_radius - minimum_radius) >= 0 || throw(DomainError(
        wire.r,
        "one circular wire does not fit inside the sector boundary"
    ))
    best = nothing
    best_score = nothing
    tolerance = 256 * _geometry_tolerance(maximum_radius)
    for pitch_fraction in range(0.0, 1.0; length = 9)
        radial_pitch = minimum_pitch +
                       pitch_fraction * (maximum_pitch - minimum_pitch)
        ratio = _geometry_scalar(
            (maximum_radius - minimum_radius) / radial_pitch
        )
        outer_phase = ratio - floor(ratio)
        phases = collect(range(0.0, 1.0; length = 9))[1:(end - 1)]
        any(phase -> isapprox(
                phase, outer_phase; atol = 64eps(Float64), rtol = 0
            ), phases) || push!(phases, outer_phase)
        for phase in phases, odd_first in (true, false)
            radii = typeof(minimum_radius)[]
            radius = minimum_radius + phase * radial_pitch
            while _geometry_scalar(radius - maximum_radius) <= tolerance
                sector_limit(shape, radius, wire.r) === nothing || push!(radii, radius)
                radius += radial_pitch
            end
            isempty(radii) && continue
            layout = phase_layout(shape, wire, radii, odd_first)
            layout === nothing && continue
            count = sum(option -> option.count, layout)
            arc = sum(option -> option.arc, layout)
            radial_extent = last(layout).radius - first(layout).radius
            contact = sum(2:length(layout); init = zero(radial_pitch)) do shell
                -abs(shell_gap(layout[shell - 1].points, layout[shell].points) - 2wire.r)
            end
            score = (
                count,
                length(layout),
                _geometry_scalar(radial_extent),
                _geometry_scalar(arc),
                _geometry_scalar(contact),
                -pitch_fraction,
                -phase,
                odd_first
            )
            if best === nothing || score > best_score
                best = layout
                best_score = score
            end
        end
    end
    best === nothing && throw(DomainError(
        wire.r,
        "one circular wire does not fit inside the sector boundary"
    ))
    return reduce(vcat, getproperty.(best, :points))
end

function boundary_anchor(cell, boundary, centre)
    scale = maximum(point -> hypot(point...), boundary)
    tolerance = 4096eps(float(scale)) * max(scale, one(scale))
    candidates = eltype(cell)[]
    for point in cell
        for index in eachindex(boundary)
            next = mod1(index + 1, length(boundary))
            edge = (
                boundary[next][1] - boundary[index][1],
                boundary[next][2] - boundary[index][2]
            )
            length_squared = edge[1]^2 + edge[2]^2
            iszero(length_squared) && continue
            offset = (point[1] - boundary[index][1], point[2] - boundary[index][2])
            fraction = clamp(
                (offset[1] * edge[1] + offset[2] * edge[2]) / length_squared,
                zero(length_squared), one(length_squared)
            )
            projection = (
                boundary[index][1] + fraction * edge[1],
                boundary[index][2] + fraction * edge[2]
            )
            hypot(point[1] - projection[1], point[2] - projection[2]) <= tolerance ||
                continue
            push!(candidates, point)
            break
        end
    end
    isempty(candidates) && return polygon_centroid(cell)
    index = findmax(
        point -> (point[1] - centre[1])^2 + (point[2] - centre[2])^2,
        candidates
    )[2]
    return candidates[index]
end

function area_preserving_strand(cell, source_area, boundary, boundary_centre)
    cell_area = signed_polygon_area(cell)
    source_area <= cell_area * (one(cell_area) + 4.0e-6) || throw(DomainError(
        source_area, "a source strand does not fit inside its assigned power cell"
    ))
    factor = min(one(cell_area), sqrt(source_area / cell_area))
    anchor = boundary_anchor(cell, boundary, boundary_centre)
    points = [
        (
            anchor[1] + factor * (point[1] - anchor[1]),
            anchor[2] + factor * (point[2] - anchor[2])
        ) for point in cell
    ]
    return normalize_polygon_area(points, source_area, anchor)
end

function resolved_polygon(points, angle)
    centre = polygon_centroid(points)
    cosine = cos(angle)
    sine = sin(angle)
    local_points = map(points) do point
        dx = point[1] - centre[1]
        dy = point[2] - centre[2]
        return (cosine * dx + sine * dy, -sine * dx + cosine * dy)
    end
    return _polygon(local_points, Pose2(centre[1], centre[2], angle))
end

function deform_disk_members(boundary_shape, members)
    all(member -> member.source.primitive isa Disk, members) || throw(ArgumentError(
        "boundary-constrained compaction currently requires circular source strands"
    ))
    source_areas = [area(member.source.primitive) for member in members]
    total_source_area = sum(source_areas)
    nominal_source_areas = nominal.(source_areas)
    nominal_total_source_area = sum(nominal_source_areas)
    nominal_boundary_area = nominal(area(boundary_shape))
    nominal_total_source_area <= nominal_boundary_area * (1 + 2.0e-6) ||
        throw(DomainError(
            total_source_area,
            "the compacted strand inventory exceeds its authoritative boundary area"
        ))
    points = boundary_shape isa Disk ? 256 : 64
    polygon = [nominal.(point)
               for point in boundary_polygon(boundary_shape; points)]
    while nominal_total_source_area > signed_polygon_area(polygon) * (1 + 2.0e-6)
        points *= 2
        points <= 16_384 || throw(ArgumentError(
            "bounded compaction could not resolve the declared fill factor " *
            "within its geometric tolerance"
        ))
        polygon = [nominal.(point)
                   for point in boundary_polygon(boundary_shape; points)]
    end
    targets = signed_polygon_area(polygon) .* nominal_source_areas ./
              nominal_total_source_area
    sites = [nominal.(member.site) for member in members]
    cells = compact_power_cells(polygon, sites, targets)
    centre = nominal.(centroid(boundary_shape))
    return [
        resolved_polygon(
            area_preserving_strand(
                cells[index], source_areas[index], polygon, centre
            ),
            members[index].angle
        ) for index in eachindex(members)
    ]
end
