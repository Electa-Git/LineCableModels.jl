"""
$(TYPEDEF)

Record the bounded-formation provenance of one resolved member.

The boundary identifies the complete formation while `course` records the
inferred radial course. Course zero identifies the centre strand of a circular
bundle or of one sector segment.

$(TYPEDFIELDS)
"""
struct BoundedPlacement{B <: AbstractShape}
    "Resolved boundary shared by every member of the formation."
    boundary::B
    "Inferred radial course; zero denotes a centre member."
    course::Int

    function BoundedPlacement(boundary::B, course::Integer) where {
            B <: AbstractShape
    }
        course >= 0 || throw(DomainError(
            course, "bounded-placement course must be nonnegative"
        ))
        return new{B}(boundary, Int(course))
    end
end

function resolve(at::Pose2, region::PlacedRegion)
    patterns = map(region.placement.patterns) do entry
        pattern = entry.pattern isa BoundedPlacement ?
                  BoundedPlacement(
            resolve(at, entry.pattern.boundary), entry.pattern.course
        ) : entry.pattern
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
    return clip_halfplane!(eltype(points)[], points, normal, offset)
end

function clip_halfplane!(output, points, normal, offset)
    empty!(output)
    isempty(points) && return output
    scale = max(abs(offset), one(offset))
    tolerance = 128eps(float(scale)) * scale
    previous = last(points)
    previous_value = normal[1] * previous[1] + normal[2] * previous[2] - offset
    previous_inside = previous_value <= tolerance
    for current in points
        current_value = normal[1] * current[1] + normal[2] * current[2] - offset
        current_inside = current_value <= tolerance
        if current_inside != previous_inside
            denominator = previous_value - current_value
            fraction = iszero(denominator) ? zero(denominator) :
                       clamp(previous_value / denominator,
                zero(denominator), one(denominator))
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

function balance_power_cells(boundary, sites, targets; maxiter::Integer = 50,
        rtol::Real = 1.0e-11)
    count = length(sites)
    count == length(targets) || throw(DimensionMismatch(
        "power-cell sites and target areas must have equal lengths"
    ))
    count == 1 && return ([copy(boundary)], zeros(typeof(first(targets)), 1))
    origin = polygon_centroid(boundary)
    length_scale = sqrt(sum(targets))
    boundary = [((point[1] - origin[1]) / length_scale,
                 (point[2] - origin[2]) / length_scale) for point in boundary]
    sites = [((point[1] - origin[1]) / length_scale,
              (point[2] - origin[2]) / length_scale) for point in sites]
    targets = targets ./ length_scale^2
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
        if relative_error <= rtol
            resolved = [[(origin[1] + length_scale * point[1],
                          origin[2] + length_scale * point[2]) for point in cell]
                        for cell in cells]
            return (resolved, weights .* length_scale^2)
        end
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
        "was $(best_error), target is $rtol"
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

function course_count(available_area, strand_area)
    ratio = nominal(available_area / strand_area)
    capacity = ratio + 64eps(float(ratio))
    capacity >= 6 || throw(DomainError(
        available_area,
        "a stranded formation requires enough area for its first six-wire course"
    ))
    courses = floor(Int, (sqrt(1 + 4ratio / 3) - 1) / 2)
    # The inverse can round below an integer at exact full occupancy. Check
    # the actual inventory in both directions, using the same area tolerance.
    while 3(courses + 1) * (courses + 2) <= capacity
        courses += 1
    end
    while 3courses * (courses + 1) > capacity
        courses -= 1
    end
    return courses
end

function circular_courses(shape::Disk, centre::Disk, wire::Disk, compact::Bool)
    courses = if compact
        course_count(area(shape) - area(centre), area(wire))
    else
        ratio = nominal((shape.r - centre.r) / (2wire.r))
        floor(Int, ratio + 64eps(float(ratio)))
    end
    courses >= 1 || throw(DomainError(
        wire.r, "the disk boundary admits no complete six-wire course"
    ))
    source_area = area(centre) + 3courses * (courses + 1) * area(wire)
    sites = NamedTuple[]
    for course in 1:courses
        count = 6course
        radius = if compact
            area_before = area(centre) + 3(course - 1) * course * area(wire)
            area_middle = area_before + 3course * area(wire)
            shape.r * sqrt(area_middle / source_area)
        else
            centre.r + (2course - 1) * wire.r
        end
        phase = isodd(course) ? zero(radius) : π / count
        for member in 1:count
            angle = phase + 2π * (member - 1) / count
            push!(sites, (
                site = (
                    shape.at.x + radius * cos(angle),
                    shape.at.y + radius * sin(angle)
                ),
                course,
                member,
                angle
            ))
        end
    end
    return sites
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

function sector_polygon(shape::SectorShape; points::Integer = 17)
    polygon = boundary_polygon(shape; points)
    cosine = cos(shape.at.φ)
    sine = sin(shape.at.φ)
    first_index = findmin(polygon) do point
        dx = point[1] - shape.at.x
        dy = point[2] - shape.at.y
        cosine * dx + sine * dy
    end[2]
    return [polygon[mod1(first_index + index - 1, length(polygon))]
            for index in eachindex(polygon)]
end

function fan_measure(polygon, centre)
    areas = map(eachindex(polygon)) do index
        next = mod1(index + 1, length(polygon))
        first_point = polygon[index]
        last_point = polygon[next]
        abs(
            (first_point[1] - centre[1]) * (last_point[2] - centre[2]) -
            (last_point[1] - centre[1]) * (first_point[2] - centre[2])
        ) / 2
    end
    total = sum(areas)
    total > zero(total) || throw(ArgumentError(
        "a sector course requires a positive-area boundary"
    ))
    return [zero(total); cumsum(areas) ./ total]
end

function fan_point(polygon, measure, fraction)
    wrapped = mod(fraction, one(fraction))
    index = min(searchsortedlast(measure, wrapped), length(polygon))
    width = measure[index + 1] - measure[index]
    local_fraction = iszero(width) ? zero(width) :
                     (wrapped - measure[index]) / width
    next = mod1(index + 1, length(polygon))
    return (
        polygon[index][1] + local_fraction *
        (polygon[next][1] - polygon[index][1]),
        polygon[index][2] + local_fraction *
        (polygon[next][2] - polygon[index][2])
    )
end
"""
$(TYPEDSIGNATURES)

Map one centre strand and complete circular `6k` courses into a sector.
The mapped sites retain their course identities while a prescribed-area power
diagram allocates space; clipped disks supply the actual conductor shapes.

For `L` courses and `N = 1 + 3L(L+1)` strands, the course sites are

```math
p_0=c,\\qquad
p_{k,m}=c+\\sqrt{\\frac{1+3k^2}{N}}\\,[q(t_{k,m})-c],
\\qquad t_{k,m}=\\frac{m+\\delta_k}{6k}.
```

Here `c` is the polygonal sector centroid and `q` traverses the boundary
by normalized swept area. Sites remain fixed during power-cell balancing.

# Arguments

- `shape`: Resolved sector boundary, with coordinates in metres.
- `wire`: Circular source strand; its area is preserved in every member.

# Returns

- Mapped members, resolved conductor polygons, and the number of outer courses.
"""
function sector_courses(shape::SectorShape, wire::Disk)
    source_area = area(wire)
    courses = course_count(area(shape) - source_area, source_area)
    total_count = 1 + 3courses * (courses + 1)
    points = 17
    polygon = sector_polygon(shape; points)
    while total_count * nominal(source_area) >
          signed_polygon_area([nominal.(point) for point in polygon]) * (1 + 2.0e-6)
        points = 2points - 1
        points <= 16_385 || throw(ArgumentError(
            "sector course tessellation could not preserve its source inventory"
        ))
        polygon = sector_polygon(shape; points)
    end
    polygon = [nominal.(point) for point in polygon]
    centre = polygon_centroid(polygon)
    measure = fan_measure(polygon, centre)
    members = [(site = centre, course = 0, member = 1, angle = zero(shape.at.φ))]
    for course in 1:courses
        count = 6course
        scale = sqrt((1 + 3course^2) / total_count)
        phase = isodd(course) ? 0.0 : 0.5
        for member in 1:count
            point = fan_point(polygon, measure, (member - 1 + phase) / count)
            site = (
                centre[1] + scale * (point[1] - centre[1]),
                centre[2] + scale * (point[2] - centre[2])
            )
            push!(members, (
                site,
                course,
                member,
                angle = atan(site[2] - centre[2], site[1] - centre[1])
            ))
        end
    end
    targets = fill(signed_polygon_area(polygon) / total_count, total_count)
    cells, _ = balance_power_cells(polygon, getproperty.(members, :site), targets)
    primitives = [
        resolved_polygon(
            area_preserving_strand(cell, source_area; angle = shape.at.φ),
            zero(shape.at.φ)
        ) for cell in cells
    ]
    return members, primitives, courses
end

function clipped_disk!(points, scratch, cell, radius, directions)
    empty!(points)
    for point in directions
        push!(points, (radius * point[1], radius * point[2]))
    end
    for index in eachindex(cell)
        next = mod1(index + 1, length(cell))
        first_point, last_point = cell[index], cell[next]
        normal = (last_point[2] - first_point[2], first_point[1] - last_point[1])
        offset = normal[1] * first_point[1] + normal[2] * first_point[2]
        clip_halfplane!(scratch, points, normal, offset)
        points, scratch = scratch, points
    end
    return points
end

"""
$(TYPEDSIGNATURES)

Grow a disk about a convex cell's centroid until its intersection with the cell
has the prescribed strand area:

```math
W_i=P_i\\cap D(c_i,\\rho_i),\\qquad |W_i|=A_i,
\\qquad c_i=\\operatorname{centroid}(P_i).
```

The disk is approximated by a regular polygon. Bisection solves its radius
against the area of the clipped polygon, so the area tolerance is independent
of the arc tessellation. Cell faces bound contact regions; free boundaries
follow the disk. Full occupancy returns the cell within the boundary's
geometric tolerance.

# Arguments

- `cell`: Counter-clockwise vertices of an allocated convex cell \\[m\\].
- `source_area`: Prescribed transverse conductor area \\[m²\\].

# Keywords

- `angle=0`: Orientation of the disk tessellation \\[rad\\].
- `points=128`: Vertices in the disk approximation.

# Returns

- Vertices of the strand polygon, contained in its cell to numerical tolerance.
"""
function area_preserving_strand(cell, source_area; angle = 0, points::Integer = 128)
    points >= 16 || throw(ArgumentError("a clipped disk requires at least 16 vertices"))
    target = nominal(source_area)
    cell_area = signed_polygon_area(cell)
    target > 0 || throw(DomainError(source_area, "strand area must be positive"))
    target <= cell_area * (1 + 4.0e-6) || throw(DomainError(
        source_area, "a source strand does not fit inside its allocated cell"
    ))
    centre = polygon_centroid(cell)
    if target >= cell_area * (1 - 64eps(float(cell_area)))
        return normalize_polygon_area(cell, source_area, centre)
    end

    # Solve at unit area to avoid coordinate-scale tolerances and cancellation.
    scale = sqrt(target)
    local_cell = [((point[1] - centre[1]) / scale,
                   (point[2] - centre[2]) / scale) for point in cell]
    directions = [(cos(angle + 2pi * index / points),
                   sin(angle + 2pi * index / points)) for index in 0:(points - 1)]
    lower = inv(sqrt(pi))
    upper = maximum(point -> hypot(point...), local_cell) / cos(pi / points)
    radius = lower
    buffer = similar(directions, 0)
    scratch = similar(buffer)
    sizehint!(buffer, points + length(cell))
    sizehint!(scratch, points + length(cell))
    strand = buffer
    for _ in 1:64
        radius = (lower + upper) / 2
        strand = clipped_disk!(buffer, scratch, local_cell, radius, directions)
        value = signed_polygon_area(strand)
        abs(value - 1) <= 64eps(float(value)) && break
        value < 1 ? (lower = radius) : (upper = radius)
    end

    # Differentiate the scalar area constraint for uncertainty-bearing inputs;
    # the nominal cell topology is fixed, as in the existing compaction model.
    if !iszero(source_area - target)
        step = cbrt(eps(float(radius))) * radius
        area_plus = signed_polygon_area(
            clipped_disk!(buffer, scratch, local_cell, radius + step, directions)
        )
        area_minus = signed_polygon_area(
            clipped_disk!(buffer, scratch, local_cell, radius - step, directions)
        )
        derivative = (area_plus - area_minus) / (2step)
        radius += (source_area / target - 1) / derivative
        point_type = typeof((radius * first(directions)[1], radius * first(directions)[2]))
        strand = clipped_disk!(point_type[], point_type[], local_cell, radius, directions)
    end
    return [(centre[1] + scale * point[1], centre[2] + scale * point[2])
            for point in strand]
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
    return [
        resolved_polygon(
            area_preserving_strand(
                cells[index], source_areas[index]; angle = boundary_shape.at.φ
            ),
            members[index].angle
        ) for index in eachindex(members)
    ]
end
