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

function compact_power_cells(boundary, sites, targets; relaxations::Integer = 4)
    current_sites = copy(sites)
    cells = Vector{eltype(boundary)}[]
    for pass in 1:Int(relaxations)
        cells, _ = balance_power_cells(boundary, current_sites, targets)
        pass == relaxations && break
        current_sites = polygon_centroid.(cells)
    end
    return cells
end

function sector_course_sites(shape::SectorShape, members, polygon)
    centre = centroid(shape)
    maximum_course = maximum(member -> member.course, members)
    maximum_course > 0 || throw(ArgumentError(
        "a sector stranded formation requires at least one strand course"
    ))
    scale = maximum(point -> hypot(point[1] - centre[1], point[2] - centre[2]), polygon)
    tolerance = 256eps(float(scale)) * max(scale, one(scale))
    return map(members) do member
        member.course > 0 || throw(ArgumentError(
            "a sector stranded formation cannot contain a central member"
        ))
        angle = shape.at.φ + member.angle
        direction = (cos(angle), sin(angle))
        extent = oftype(scale, Inf)
        for index in eachindex(polygon)
            next = mod1(index + 1, length(polygon))
            edge = (
                polygon[next][1] - polygon[index][1],
                polygon[next][2] - polygon[index][2]
            )
            denominator = edge[1] * direction[2] - edge[2] * direction[1]
            denominator < -tolerance || continue
            offset = (
                centre[1] - polygon[index][1],
                centre[2] - polygon[index][2]
            )
            numerator = edge[1] * offset[2] - edge[2] * offset[1]
            candidate = numerator / -denominator
            candidate >= -tolerance || continue
            extent = min(extent, max(candidate, zero(candidate)))
        end
        isfinite(extent) && extent > tolerance || throw(ArgumentError(
            "a sector strand course could not be mapped inside its boundary"
        ))
        fraction = member.course / (maximum_course + 0.35)
        return (
            centre[1] + fraction * extent * direction[1],
            centre[2] + fraction * extent * direction[2]
        )
    end
end

function disk_course_sites(shape::Disk, members)
    centre = (shape.at.x, shape.at.y)
    maximum_course = maximum(member -> member.course, members)
    maximum_course == 0 && return [centre]
    sites = map(members) do member
        member.course == 0 && return centre
        radial_fraction = member.course / (maximum_course + 0.35)
        radius = 0.92shape.r * radial_fraction
        angle = member.angle
        return (centre[1] + radius * cos(angle), centre[2] + radius * sin(angle))
    end
    return sites
end

function compaction_factor(compact::FillFactor, source_area, boundary_area)
    derived = source_area / boundary_area
    positive = derived > zero(derived)
    admissible = derived <= one(derived) ||
                 isapprox(derived, one(derived); rtol = 2.0e-6, atol = zero(derived))
    positive && admissible || throw(DomainError(
        derived,
        "bounded strand area must be positive and cannot exceed its boundary area"
    ))
    derived = min(derived, one(derived))
    isapprox(compact.η, derived; rtol = 2.0e-6, atol = zero(derived)) ||
        throw(DomainError(
            compact.η,
            "the declared fill factor conflicts with strand area / boundary area = $derived"
        ))
    return derived
end

function natural_disk_members(shape::Disk, members)
    centre = (shape.at.x, shape.at.y)
    resolved = Vector{AbstractShape}(undef, length(members))
    outer_radius = zero(shape.r)
    for course in 0:maximum(member -> member.course, members)
        indices = findall(member -> member.course == course, members)
        isempty(indices) && throw(ArgumentError(
            "a natural circular formation must contain every declared course"
        ))
        radii = [members[index].source.primitive.r for index in indices]
        radius = first(radii)
        all(value -> isapprox(value, radius), radii) || throw(ArgumentError(
            "one natural circular strand course requires equal wire radii"
        ))
        if iszero(course)
            length(indices) == 1 || throw(ArgumentError(
                "a natural circular formation requires exactly one central wire"
            ))
            resolved[only(indices)] = Disk(radius, Pose2(centre[1], centre[2], shape.at.φ))
            outer_radius = radius
            continue
        end
        lay_radius = outer_radius + radius
        count = length(indices)
        chord = count == 1 ? oftype(lay_radius, Inf) :
                2lay_radius * sin(π / count)
        tolerance = sqrt(eps(float(chord))) * max(one(chord), chord)
        chord + tolerance >= 2radius || throw(DomainError(
            count,
            "natural circular strands overlap in course $course"
        ))
        for index in indices
            angle = shape.at.φ + members[index].angle
            pose = Pose2(
                centre[1] + lay_radius * cos(angle),
                centre[2] + lay_radius * sin(angle),
                angle
            )
            resolved[index] = resolve(pose, members[index].source.primitive)
        end
        outer_radius = lay_radius + radius
    end
    tolerance = sqrt(eps(float(shape.r))) * max(one(shape.r), shape.r)
    outer_radius <= shape.r + tolerance || throw(DomainError(
        shape.r,
        "the natural circular strands require radius $outer_radius"
    ))
    return resolved
end

function circle_inside_polygon(point, radius, polygon, tolerance)
    for index in eachindex(polygon)
        next = mod1(index + 1, length(polygon))
        edge = (
            polygon[next][1] - polygon[index][1],
            polygon[next][2] - polygon[index][2]
        )
        offset = (
            point[1] - polygon[index][1],
            point[2] - polygon[index][2]
        )
        clearance = edge[1] * offset[2] - edge[2] * offset[1]
        clearance + tolerance >= radius * hypot(edge...) || return false
    end
    return true
end

function order_sector_sites(sites, centre)
    return sort(sites; by = point -> (
        atan(point[2] - centre[2], point[1] - centre[1]),
        hypot(point[1] - centre[1], point[2] - centre[2])
    ))
end

function select_sector_sites(sites, count::Integer, centre)
    count > 0 || throw(ArgumentError(
        "a sector stranded formation requires at least one wire"
    ))
    count <= length(sites) || throw(DomainError(
        count,
        "the sector admits at most $(length(sites)) natural circular wires"
    ))
    count == length(sites) && return order_sector_sites(sites, centre)
    available = collect(eachindex(sites))
    first_index = findmin(index -> hypot(
        sites[index][1] - centre[1], sites[index][2] - centre[2]
    ), available)[2]
    selected = Int[available[first_index]]
    deleteat!(available, first_index)
    while length(selected) < count
        position = findmax(available) do candidate
            minimum(selected) do accepted
                (sites[candidate][1] - sites[accepted][1])^2 +
                (sites[candidate][2] - sites[accepted][2])^2
            end
        end[2]
        push!(selected, available[position])
        deleteat!(available, position)
    end
    return order_sector_sites(sites[selected], centre)
end

function sector_wire_sites(shape::SectorShape, wire::Disk)
    polygon = boundary_polygon(shape; points = 128)
    centre = centroid(shape)
    spacing = 2wire.r
    scale = maximum(point -> hypot(
        point[1] - centre[1], point[2] - centre[2]
    ), polygon)
    tolerance = 512eps(float(scale)) * max(scale, one(scale))
    limit = ceil(Int, 2nominal(scale / spacing)) + 3
    best = Tuple{typeof(centre[1]), typeof(centre[2])}[]
    best_offset = oftype(float(nominal(scale)), Inf)
    for orientation in range(shape.at.φ, shape.at.φ + π / 6; length = 4)
        first_basis = (spacing * cos(orientation), spacing * sin(orientation))
        second_basis = (
            spacing * cos(orientation + π / 3),
            spacing * sin(orientation + π / 3)
        )
        for first_phase in 0:4, second_phase in 0:4
            first_fraction = first_phase / 5
            second_fraction = second_phase / 5
            sites = Tuple{typeof(centre[1]), typeof(centre[2])}[]
            for first_index in -limit:limit, second_index in -limit:limit
                first_coordinate = first_index + first_fraction
                second_coordinate = second_index + second_fraction
                point = (
                    centre[1] + first_coordinate * first_basis[1] +
                    second_coordinate * second_basis[1],
                    centre[2] + first_coordinate * first_basis[2] +
                    second_coordinate * second_basis[2]
                )
                circle_inside_polygon(point, wire.r, polygon, tolerance) &&
                    push!(sites, point)
            end
            isempty(sites) && continue
            site_centre = (
                sum(first, sites) / length(sites),
                sum(last, sites) / length(sites)
            )
            offset = hypot(
                site_centre[1] - centre[1], site_centre[2] - centre[2]
            )
            if length(sites) > length(best) ||
               (length(sites) == length(best) && offset < best_offset)
                best = sites
                best_offset = offset
            end
        end
    end
    isempty(best) && throw(DomainError(
        wire.r,
        "one circular wire does not fit inside the sector boundary"
    ))
    return order_sector_sites(best, centre)
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
    return _polygon(Tuple(local_points), Pose2(centre[1], centre[2], angle))
end

function compact_disk_members(boundary_shape, members, compact)
    all(member -> member.source.primitive isa Disk, members) || throw(ArgumentError(
        "boundary-constrained compaction currently requires circular source strands"
    ))
    source_areas = [area(member.source.primitive) for member in members]
    total_source_area = sum(source_areas)
    compaction_factor(compact, total_source_area, area(boundary_shape))
    points = boundary_shape isa Disk ? 256 : 64
    coordinate_type = foldl(
        promote_type,
        (typeof(value) for value in source_areas);
        init = eltype(boundary_shape)
    )
    polygon = [convert.(coordinate_type, point)
               for point in boundary_polygon(boundary_shape; points)]
    while total_source_area > signed_polygon_area(polygon) * (1 + 2.0e-6)
        points *= 2
        points <= 16_384 || throw(ArgumentError(
            "bounded compaction could not resolve the declared fill factor " *
            "within its geometric tolerance"
        ))
        polygon = [convert.(coordinate_type, point)
                   for point in boundary_polygon(boundary_shape; points)]
    end
    targets = signed_polygon_area(polygon) .* source_areas ./ total_source_area
    sites = if boundary_shape isa Disk
        [convert.(coordinate_type, point)
         for point in disk_course_sites(boundary_shape, members)]
    else
        [convert.(coordinate_type, point)
         for point in sector_course_sites(boundary_shape, members, polygon)]
    end
    cells = compact_power_cells(polygon, sites, targets)
    centre = centroid(boundary_shape)
    return [
        resolved_polygon(
            area_preserving_strand(
                cells[index], source_areas[index], polygon, centre
            ),
            members[index].angle
        ) for index in eachindex(members)
    ]
end

function compact_bounded_members(boundary, members, compact, at::Pose2)
    boundary_shape = resolve(at, boundary)
    source_shapes = map(member -> member.source.primitive, members)
    all(shape -> shape isa Disk, source_shapes) || throw(ArgumentError(
        "stranded formations require circular source wires"
    ))
    resolved = if compact === nothing && boundary_shape isa Disk
        natural_disk_members(boundary_shape, members)
    elseif compact === nothing && boundary_shape isa SectorShape
        [resolve(Pose2(member.site..., member.angle), member.source.primitive)
         for member in members]
    else
        compact isa FillFactor || throw(ArgumentError(
            "bounded strand deformation requires one FillFactor"
        ))
        compact_disk_members(boundary_shape, members, compact)
    end
    return resolved, boundary_shape
end
