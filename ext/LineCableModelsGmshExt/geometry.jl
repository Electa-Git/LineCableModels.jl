struct FEMLoop
    ccw::Int
    cw::Int
    curves::Vector{Int}
    oriented::Vector{Int}
end

const FEM_BOUNDARY_ANGLE_TOLERANCE = 1.0e-12

mutable struct FEMLoopRegistry
    points::Dict{Tuple{Float64, Float64}, Int}
    point_buckets::Dict{Tuple{Int, Int}, Vector{Tuple{Float64, Float64}}}
    point_sizes::Dict{Int, Float64}
    lines::Dict{Tuple{Int, Int}, Int}
    curve_points::Dict{Int, Tuple{Int, Int}}
    curve_samples::Dict{Int, Vector{Tuple{Float64, Float64}}}
    circle_arcs::Dict{Any, Tuple{Int, Int, Int}}
    circle_breaks::Dict{Any, Set{Float64}}
    circle_break_points::Dict{Any, Dict{Float64, Tuple{Float64, Float64}}}
    loops::Dict{Any, FEMLoop}
    mesh_size::Float64
end

function FEMLoopRegistry(mesh_size)
    FEMLoopRegistry(
        Dict{Tuple{Float64, Float64}, Int}(),
        Dict{Tuple{Int, Int}, Vector{Tuple{Float64, Float64}}}(),
        Dict{Int, Float64}(),
        Dict{Tuple{Int, Int}, Int}(),
        Dict{Int, Tuple{Int, Int}}(),
        Dict{Int, Vector{Tuple{Float64, Float64}}}(),
        Dict{Any, Tuple{Int, Int, Int}}(),
        Dict{Any, Set{Float64}}(),
        Dict{Any, Dict{Float64, Tuple{Float64, Float64}}}(),
        Dict{Any, FEMLoop}(),
        Float64(mesh_size)
    )
end

_coordinate_key(value) = round(Float64(value); sigdigits = 15)

function _matching_point_key(registry::FEMLoopRegistry, point)
    points = registry.points
    key = (_coordinate_key(point[1]), _coordinate_key(point[2]))
    haskey(points, key) && return key
    scale = max(abs(Float64(point[1])), abs(Float64(point[2])), 1.0)
    tolerance = 64eps(scale)
    for x in floor(Int, (key[1] - tolerance) / registry.mesh_size):floor(Int, (key[1] + tolerance) / registry.mesh_size),
        y in floor(Int, (key[2] - tolerance) / registry.mesh_size):floor(Int, (key[2] + tolerance) / registry.mesh_size)
        for candidate in get(registry.point_buckets, (x, y), ())
            abs(candidate[1] - key[1]) <= tolerance &&
                abs(candidate[2] - key[2]) <= tolerance && return candidate
        end
    end
    return key
end

function _circle_key(centre, radius)
    return (
        _coordinate_key(centre[1]),
        _coordinate_key(centre[2]),
        _coordinate_key(radius)
    )
end

function _angle_key(angle)
    value = mod(Float64(angle), 2π)
    isapprox(value, 2π; rtol = 0, atol = 128eps(Float64)) && (value = 0.0)
    return value
end

function _matching_circle_break(breaks, angle)
    for existing in breaks
        difference = abs(existing - angle)
        min(difference, 2π - difference) <= FEM_BOUNDARY_ANGLE_TOLERANCE &&
            return existing
    end
    return nothing
end

function _register_circle_break!(
        registry::FEMLoopRegistry,
        centre,
        radius,
        angle;
        point = nothing
)
    key = _circle_key(centre, radius)
    breaks = get!(registry.circle_breaks, key) do
        Set{Float64}()
    end
    value = _angle_key(angle)
    stored = something(_matching_circle_break(breaks, value), value)
    push!(breaks, stored)
    if point !== nothing
        points = get!(registry.circle_break_points, key) do
            Dict{Float64, Tuple{Float64, Float64}}()
        end
        get!(points, stored) do
            (Float64(point[1]), Float64(point[2]))
        end
    end
    return stored
end

function _register_full_circle_breaks!(registry::FEMLoopRegistry, centre, radius)
    foreach(
        angle -> _register_circle_break!(registry, centre, radius, angle),
        (0.0, π / 2, π, 3π / 2)
    )
    return nothing
end

function _tangent_circular_fill_plan(shape)
    shape isa DataModel.DifferenceShape || return nothing
    outer = shape.outer
    outer isa DataModel.Annulus || return nothing
    isempty(shape.holes) && return nothing
    iszero(outer.ri) && return nothing

    centre = (Float64(outer.at.x), Float64(outer.at.y))
    disks = filter(hole -> hole isa DataModel.Disk, shape.holes)
    cutouts = filter(hole -> !(hole isa DataModel.Disk), shape.holes)
    isempty(disks) && return nothing
    all(hole -> hole isa DataModel.BentStrip, cutouts) || return nothing
    entries = Tuple{Float64, Any}[]
    outer_contacts = Float64[]
    for hole in disks
        dx = Float64(hole.at.x) - centre[1]
        dy = Float64(hole.at.y) - centre[2]
        distance = hypot(dx, dy)
        radius = Float64(hole.r)
        scale = max(Float64(outer.ro), distance, radius, 1.0)
        tolerance = 1.0e-12 * scale
        distance > tolerance || return nothing
        isapprox(
            distance - radius,
            Float64(outer.ri);
            rtol = 0,
            atol = tolerance
        ) || return nothing
        push!(outer_contacts, distance + radius)
        push!(entries, (_angle_key(atan(dy, dx)), hole))
    end
    computed_outer_radius = first(outer_contacts)
    tolerance = 1.0e-12 * max(
        Float64(outer.ro), computed_outer_radius, 1.0
    )
    all(
        radius -> isapprox(
            radius, computed_outer_radius; rtol = 0, atol = tolerance
        ), outer_contacts) || return nothing
    computed_outer_radius <= Float64(outer.ro) + tolerance || return nothing

    # Prefer the radius declared by the enclosing material boundary.  The
    # contact radius obtained from distance + wire radius can differ by an ULP;
    # using it as a separate circle would make coincident material interfaces
    # overlap instead of sharing the same Gmsh curves.
    tangent_outer_radius = if !isempty(cutouts)
        Float64(first(cutouts).ri)
    elseif isapprox(
        computed_outer_radius, Float64(outer.ro); rtol = 0, atol = tolerance
    )
        Float64(outer.ro)
    else
        computed_outer_radius
    end
    isapprox(
        tangent_outer_radius, computed_outer_radius;
        rtol = 0,
        atol = tolerance
    ) || return nothing

    for cutout in cutouts
        isapprox(Float64(cutout.at.x), centre[1]; rtol = 0, atol = tolerance) ||
            return nothing
        isapprox(Float64(cutout.at.y), centre[2]; rtol = 0, atol = tolerance) ||
            return nothing
        isapprox(
            Float64(cutout.ri), tangent_outer_radius;
            rtol = 0,
            atol = tolerance
        ) || return nothing
        isapprox(
            Float64(cutout.ro), Float64(outer.ro);
            rtol = 0,
            atol = tolerance
        ) || return nothing
        FEM_BOUNDARY_ANGLE_TOLERANCE < Float64(cutout.span) <
        2π - FEM_BOUNDARY_ANGLE_TOLERANCE || return nothing
    end
    has_outer_shell = Float64(outer.ro) - tangent_outer_radius > tolerance
    has_outer_shell || isempty(cutouts) || return nothing
    outer_spans = has_outer_shell ? _complementary_sector_spans(cutouts) :
                  Tuple{Float64, Float64}[]
    outer_spans === nothing && return nothing

    sort!(entries; by = first)
    if length(entries) > 1
        angles = first.(entries)
        separations = [mod(angles[mod1(index + 1, end)] - angles[index], 2π)
                       for index in eachindex(angles)]
        minimum(separations) > FEM_BOUNDARY_ANGLE_TOLERANCE || return nothing
    end
    return (
        centre = centre,
        inner_radius = Float64(outer.ri),
        outer_radius = tangent_outer_radius,
        shell_outer_radius = Float64(outer.ro),
        angles = first.(entries),
        holes = last.(entries),
        outer_spans
    )
end

function _complementary_sector_spans(cutouts)
    isempty(cutouts) && return Tuple{Float64, Float64}[(0.0, 2π)]
    occupied = Tuple{Float64, Float64}[]
    for cutout in cutouts
        start = _angle_key(Float64(cutout.at.φ - cutout.span / 2))
        stop = start + Float64(cutout.span)
        if stop <= 2π + FEM_BOUNDARY_ANGLE_TOLERANCE
            push!(occupied, (start, min(stop, 2π)))
        else
            push!(occupied, (start, 2π))
            push!(occupied, (0.0, stop - 2π))
        end
    end
    sort!(occupied; by = first)
    merged = Tuple{Float64, Float64}[]
    for interval in occupied
        if isempty(merged) || interval[1] > last(merged)[2] +
                                            FEM_BOUNDARY_ANGLE_TOLERANCE
            push!(merged, interval)
            continue
        end
        interval[1] < last(merged)[2] - FEM_BOUNDARY_ANGLE_TOLERANCE &&
            return nothing
        merged[end] = (last(merged)[1], max(last(merged)[2], interval[2]))
    end
    spans = Tuple{Float64, Float64}[]
    cursor = 0.0
    for interval in merged
        interval[1] - cursor > FEM_BOUNDARY_ANGLE_TOLERANCE &&
            push!(spans, (cursor, interval[1] - cursor))
        cursor = max(cursor, interval[2])
    end
    2π - cursor > FEM_BOUNDARY_ANGLE_TOLERANCE &&
        push!(spans, (cursor, 2π - cursor))
    return spans
end

function _radial_point(centre, radius, angle)
    return (
        centre[1] + radius * cos(angle),
        centre[2] + radius * sin(angle)
    )
end

function _annular_sector_surface!(
        registry::FEMLoopRegistry,
        centre,
        inner_radius,
        outer_radius,
        start_angle,
        span,
        mesh_size
)
    if isapprox(span, 2π; rtol = 0, atol = FEM_BOUNDARY_ANGLE_TOLERANCE)
        outer = _circle_loop!(
            registry, centre, outer_radius; mesh_size
        )
        inner = _circle_loop!(
            registry, centre, inner_radius; mesh_size
        )
        return gmsh.model.geo.add_plane_surface([outer.ccw, inner.cw])
    end

    outer_start_point = _radial_point(centre, outer_radius, start_angle)
    outer_stop_point = _radial_point(
        centre, outer_radius, start_angle + span
    )
    inner_start_point = _radial_point(centre, inner_radius, start_angle)
    inner_stop_point = _radial_point(
        centre, inner_radius, start_angle + span
    )
    outer_start = _point!(registry, outer_start_point; mesh_size)
    outer_stop = _point!(registry, outer_stop_point; mesh_size)
    inner_start = _point!(registry, inner_start_point; mesh_size)
    inner_stop = _point!(registry, inner_stop_point; mesh_size)
    curves = _circle_arc_path!(
        registry,
        centre,
        outer_radius,
        start_angle,
        span;
        mesh_size,
        first_point = outer_start_point,
        last_point = outer_stop_point
    )
    append!(curves, _line_path!(registry, outer_stop_point, inner_stop_point; mesh_size))
    append!(curves,
        _circle_arc_path!(
            registry,
            centre,
            inner_radius,
            start_angle + span,
            -span;
            mesh_size,
            first_point = inner_stop_point,
            last_point = inner_start_point
        ))
    append!(curves, _line_path!(registry, inner_start_point, outer_start_point; mesh_size))
    loop = gmsh.model.geo.add_curve_loop(curves)
    return gmsh.model.geo.add_plane_surface([loop])
end

function _register_tangent_fill_contacts!(registry::FEMLoopRegistry, shape)
    plan = _tangent_circular_fill_plan(shape)
    plan === nothing && return nothing
    for (angle, hole) in zip(plan.angles, plan.holes)
        outer_point = _radial_point(plan.centre, plan.outer_radius, angle)
        inner_point = _radial_point(plan.centre, plan.inner_radius, angle)
        hole_centre = (Float64(hole.at.x), Float64(hole.at.y))
        hole_radius = Float64(hole.r)
        _register_circle_break!(
            registry,
            plan.centre,
            plan.outer_radius,
            angle;
            point = outer_point
        )
        _register_circle_break!(
            registry,
            plan.centre,
            plan.inner_radius,
            angle;
            point = inner_point
        )
        _register_circle_break!(
            registry,
            hole_centre,
            hole_radius,
            angle;
            point = outer_point
        )
        _register_circle_break!(
            registry,
            hole_centre,
            hole_radius,
            angle + π;
            point = inner_point
        )
        _register_circle_break!(
            registry,
            hole_centre,
            hole_radius,
            angle + π / 2
        )
        _register_circle_break!(
            registry,
            hole_centre,
            hole_radius,
            angle - π / 2
        )
    end
    return plan
end

_register_shape_breaks!(::FEMLoopRegistry, shape) = nothing

function _register_shape_breaks!(registry::FEMLoopRegistry, shape::DataModel.Polygon)
    # Register topology first; the consuming material supplies the mesh size.
    foreach(point -> _point!(registry, point; mesh_size = nothing), _shape_points(shape))
    return nothing
end

function _register_shape_breaks!(registry::FEMLoopRegistry, shape::DataModel.Disk)
    _register_full_circle_breaks!(registry, (shape.at.x, shape.at.y), shape.r)
    return nothing
end

function _register_shape_breaks!(registry::FEMLoopRegistry, shape::DataModel.Annulus)
    centre = (shape.at.x, shape.at.y)
    _register_full_circle_breaks!(registry, centre, shape.ri)
    _register_full_circle_breaks!(registry, centre, shape.ro)
    return nothing
end

function _register_shape_breaks!(registry::FEMLoopRegistry, shape::DataModel.BentStrip)
    centre = (shape.at.x, shape.at.y)
    start_angle = shape.at.φ - shape.span / 2
    stop_angle = shape.at.φ + shape.span / 2
    _register_circle_break!(registry, centre, shape.ro, start_angle)
    _register_circle_break!(registry, centre, shape.ro, stop_angle)
    if !iszero(shape.ri)
        _register_circle_break!(registry, centre, shape.ri, start_angle)
        _register_circle_break!(registry, centre, shape.ri, stop_angle)
    end
    return nothing
end

function _register_shape_breaks!(
        registry::FEMLoopRegistry,
        shape::DataModel.SectorShape
)
    for arc in (
        shape.contacts.arcs.base,
        shape.contacts.arcs.lower,
        shape.contacts.arcs.back,
        shape.contacts.arcs.upper
    )
        iszero(arc.radius) && continue
        centre = _transform_point(arc.center, shape.at)
        _register_circle_break!(registry, centre, arc.radius, shape.at.φ + arc.start)
        _register_circle_break!(registry, centre, arc.radius, shape.at.φ + arc.stop)
    end
    return nothing
end

function _register_shape_breaks!(registry::FEMLoopRegistry, shape::DataModel.ShellShape)
    _register_shape_breaks!(registry, shape.inner)
    _register_shape_breaks!(registry, shape.outer)
    return nothing
end

function _register_shape_breaks!(
        registry::FEMLoopRegistry,
        shape::DataModel.DifferenceShape
)
    _register_shape_breaks!(registry, shape.outer)
    foreach(hole -> _register_shape_breaks!(registry, hole), shape.holes)
    _register_tangent_fill_contacts!(registry, shape)
    return nothing
end

function _register_shape_breaks!(
        registry::FEMLoopRegistry,
        shape::DataModel.AssemblyShape
)
    foreach(member -> _register_shape_breaks!(registry, member), shape.members)
    return nothing
end

function _register_circle_contacts!(registry::FEMLoopRegistry)
    circles = Tuple{Float64, Float64, Float64}[key for key in keys(registry.circle_breaks)]
    # A polygonal strand can end on an analytic enclosing arc. Both consumers
    # must share that vertex before either curve loop is constructed.
    for circle in circles, point in keys(registry.points)
        dx, dy = point[1] - circle[1], point[2] - circle[2]
        abs(hypot(dx, dy) - circle[3]) <= 64eps(max(circle[3], 1.0)) || continue
        _register_circle_break!(registry, circle[1:2], circle[3], atan(dy, dx); point)
    end
    for first_index in eachindex(circles)
        first = circles[first_index]
        first_centre = (first[1], first[2])
        first_radius = first[3]
        for second_index in (first_index + 1):length(circles)
            second = circles[second_index]
            second_centre = (second[1], second[2])
            second_radius = second[3]
            dx = second_centre[1] - first_centre[1]
            dy = second_centre[2] - first_centre[2]
            distance = hypot(dx, dy)
            scale = max(first_radius, second_radius, distance, 1.0)
            tolerance = 1e-12 * scale
            distance <= tolerance && continue
            distance > first_radius + second_radius + tolerance && continue
            distance < abs(first_radius - second_radius) - tolerance && continue

            along = (
                first_radius^2 - second_radius^2 + distance^2
            ) / (2distance)
            height_squared = first_radius^2 - along^2
            height_squared < -tolerance * scale && continue
            tangent = isapprox(
                distance,
                first_radius + second_radius;
                rtol = 0,
                atol = tolerance
            ) || isapprox(
                distance,
                abs(first_radius - second_radius);
                rtol = 0,
                atol = tolerance
            )
            height = tangent ? 0.0 : sqrt(max(height_squared, 0.0))
            unit_x = dx / distance
            unit_y = dy / distance
            base_x = first_centre[1] + along * unit_x
            base_y = first_centre[2] + along * unit_y
            contacts = height <= tolerance ?
                       ((base_x, base_y),) :
                       (
                (base_x - height * unit_y, base_y + height * unit_x),
                (base_x + height * unit_y, base_y - height * unit_x)
            )
            for contact in contacts
                _register_circle_break!(
                    registry,
                    first_centre,
                    first_radius,
                    atan(
                        contact[2] - first_centre[2],
                        contact[1] - first_centre[1]
                    );
                    point = contact
                )
                _register_circle_break!(
                    registry,
                    second_centre,
                    second_radius,
                    atan(
                        contact[2] - second_centre[2],
                        contact[1] - second_centre[1]
                    );
                    point = contact
                )
            end
        end
    end
    return nothing
end

function _transform_point(point, at)
    cosine = cos(at.φ)
    sine = sin(at.φ)
    return (
        at.x + cosine * point[1] - sine * point[2],
        at.y + sine * point[1] + cosine * point[2]
    )
end

function _circle_point(centre, radius, angle)
    cosine = cos(angle)
    sine = sin(angle)
    isapprox(cosine, 0; rtol = 0, atol = 128eps(Float64)) && (cosine = 0.0)
    isapprox(sine, 0; rtol = 0, atol = 128eps(Float64)) && (sine = 0.0)
    return (
        centre[1] + radius * cosine,
        centre[2] + radius * sine
    )
end

function _point!(registry::FEMLoopRegistry, point; mesh_size = registry.mesh_size)
    # Equivalent boundary constructions can differ by one or two Float64 ULPs.
    # Reuse the existing topological point without applying a physical-scale
    # tolerance to every coordinate in the model.
    key = _matching_point_key(registry, point)
    tag = get!(registry.points, key) do
        bucket = (floor(Int, key[1] / registry.mesh_size),
                  floor(Int, key[2] / registry.mesh_size))
        push!(get!(Vector{Tuple{Float64, Float64}}, registry.point_buckets, bucket), key)
        gmsh.model.geo.add_point(
            key[1], key[2], 0.0, mesh_size === nothing ? 0.0 : Float64(mesh_size)
        )
    end
    if mesh_size !== nothing
        registry.point_sizes[tag] = min(
            get(registry.point_sizes, tag, Inf), Float64(mesh_size)
        )
    end
    return tag
end

function _line!(registry::FEMLoopRegistry, first_point::Int, last_point::Int)
    key = minmax(first_point, last_point)
    tag = get!(registry.lines, key) do
        value = gmsh.model.geo.add_line(key[1], key[2])
        registry.curve_points[value] = key
        value
    end
    return first_point == registry.curve_points[tag][1] ? tag : -tag
end

function _line_path!(registry::FEMLoopRegistry, first, last;
        mesh_size = registry.mesh_size)
    first_point = _point!(registry, first; mesh_size)
    last_point = _point!(registry, last; mesh_size)
    dx, dy = last[1] - first[1], last[2] - first[2]
    distance = hypot(dx, dy)
    distance > 0 || return Int[]
    tolerance = 64eps(max(abs(first[1]), abs(first[2]), abs(last[1]), abs(last[2]), 1.0))
    points = Tuple{Float64, Int}[(0.0, first_point), (1.0, last_point)]
    # Visit only buckets intersecting this edge's bounding box. Scanning every
    # vertex for every strand edge is quadratic for large bounded formations.
    x_bins = range(floor(Int, (min(first[1], last[1]) - tolerance) / registry.mesh_size),
                   floor(Int, (max(first[1], last[1]) + tolerance) / registry.mesh_size))
    y_bins = range(floor(Int, (min(first[2], last[2]) - tolerance) / registry.mesh_size),
                   floor(Int, (max(first[2], last[2]) + tolerance) / registry.mesh_size))
    bins = length(x_bins) <= length(registry.point_buckets) ÷ length(y_bins) ?
           Iterators.product(x_bins, y_bins) : keys(registry.point_buckets)
    for (x, y) in bins
        x in x_bins && y in y_bins || continue
        for point in get(registry.point_buckets, (x, y), ())
            tag = registry.points[point]
            tag in (first_point, last_point) && continue
            px, py = point[1] - first[1], point[2] - first[2]
            position = (px * dx + py * dy) / distance^2
            tolerance / distance < position < 1 - tolerance / distance || continue
            abs(px * dy - py * dx) <= tolerance * distance || continue
            push!(points, (position, tag))
        end
    end
    sort!(points)
    return [_line!(registry, points[index][2], points[index + 1][2])
            for index in 1:(length(points) - 1)]
end

function _refine_loop_mesh_size!(
        registry::FEMLoopRegistry, loop::FEMLoop, mesh_size
)
    value = Float64(mesh_size)
    for curve in loop.curves
        for point in registry.curve_points[abs(curve)]
            registry.point_sizes[point] = min(
                get(registry.point_sizes, point, Inf), value
            )
        end
    end
    return loop
end

function _apply_point_mesh_sizes!(registry::FEMLoopRegistry)
    for (point, mesh_size) in registry.point_sizes
        gmsh.model.geo.mesh.set_size([(0, point)], mesh_size)
    end
    return nothing
end

function _signed_area(points)
    origin = first(points)
    return sum(eachindex(points)) do index
        next = mod1(index + 1, length(points))
        (points[index][1] - origin[1]) * (points[next][2] - origin[2]) -
        (points[next][1] - origin[1]) * (points[index][2] - origin[2])
    end / 2
end

function _canonical_points(points)
    values = [(Float64(point[1]), Float64(point[2])) for point in points]
    if length(values) > 1
        first_point = first(values)
        last_point = last(values)
        hypot(first_point[1] - last_point[1], first_point[2] - last_point[2]) <=
        1e-14 && pop!(values)
    end
    length(values) >= 3 || throw(ArgumentError("a FEM boundary requires three points"))
    _signed_area(values) < 0 && reverse!(values)
    keys = [(_coordinate_key(point[1]), _coordinate_key(point[2])) for point in values]
    first_index = argmin(keys)
    values = [values[mod1(first_index + offset, length(values))]
              for offset in 0:(length(values) - 1)]
    return values
end

function _polygon_loop!(registry::FEMLoopRegistry, points; mesh_size = registry.mesh_size)
    values = _canonical_points(points)
    key = (:polygon,
        Tuple(
            (_coordinate_key(point[1]), _coordinate_key(point[2])) for point in values
        ))
    loop = get!(registry.loops, key) do
        curves = reduce(vcat,
            [_line_path!(registry, values[index], values[mod1(index + 1, end)]; mesh_size)
             for index in eachindex(values)])
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        FEMLoop(ccw, cw, abs.(curves), curves)
    end
    return _refine_loop_mesh_size!(registry, loop, mesh_size)
end

function _circle_loop!(
        registry::FEMLoopRegistry,
        centre,
        radius;
        mesh_size = registry.mesh_size
)
    key = (:circle, _circle_key(centre, radius))
    loop = get!(registry.loops, key) do
        _register_full_circle_breaks!(registry, centre, radius)
        first_point = (centre[1] + radius, centre[2])
        curves = _circle_arc_path!(
            registry,
            centre,
            radius,
            0.0,
            2π;
            mesh_size,
            first_point,
            last_point = first_point
        )
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        FEMLoop(ccw, cw, abs.(curves), curves)
    end
    return _refine_loop_mesh_size!(registry, loop, mesh_size)
end

function _ellipse_loop!(
        registry::FEMLoopRegistry,
        shape::DataModel.Ellipse;
        mesh_size = registry.mesh_size
)
    key = (
        :ellipse,
        _coordinate_key(shape.at.x),
        _coordinate_key(shape.at.y),
        _coordinate_key(shape.at.φ),
        _coordinate_key(shape.a),
        _coordinate_key(shape.b)
    )
    loop = get!(registry.loops, key) do
        centre = (shape.at.x, shape.at.y)
        centre_tag = _point!(registry, centre; mesh_size)
        local_points = (
            (shape.a, 0.0),
            (0.0, shape.b),
            (-shape.a, 0.0),
            (0.0, -shape.b)
        )
        angles = Float64[0, π/2, π, 3π/2, 2π]
        sort!(unique!(angles))
        point_tags = [_point!(registry,
            _transform_point((shape.a*cos(a),shape.b*sin(a)),shape.at); mesh_size)
            for a in angles[1:end-1]]
        major_point = shape.a >= shape.b ? local_points[1] : local_points[2]
        major_tag = _point!(
            registry, _transform_point(major_point, shape.at); mesh_size
        )
        curves = Int[]
        for index in eachindex(point_tags)
            next_index = mod1(index + 1, length(point_tags))
            curve = gmsh.model.geo.add_ellipse_arc(
                point_tags[index], centre_tag, major_tag, point_tags[next_index]
            )
            registry.curve_points[curve] = (
                point_tags[index], point_tags[next_index]
            )
            push!(curves, curve)
        end
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        for (index, curve) in enumerate(curves)
            registry.curve_samples[curve] = [
                _transform_point((shape.a * cos(angle), shape.b * sin(angle)), shape.at)
                for angle in range(angles[index], angles[index+1]; length = 17)
            ]
        end
        FEMLoop(ccw, cw, curves, curves)
    end
    return _refine_loop_mesh_size!(registry, loop, mesh_size)
end

function _circle_break_point(registry::FEMLoopRegistry, circle_key, angle)
    points = get(registry.circle_break_points, circle_key, nothing)
    points === nothing && return nothing
    stored = _matching_circle_break(keys(points), _angle_key(angle))
    stored === nothing && return nothing
    return points[stored]
end

function _circle_arc_path!(
        registry::FEMLoopRegistry,
        centre,
        radius,
        start_angle,
        span;
        mesh_size = registry.mesh_size,
        first_point = nothing,
        last_point = nothing
)
    value_span = Float64(span)
    iszero(value_span) && return Int[]
    direction = sign(value_span)
    distance = abs(value_span)
    distance <= 2π + 128eps(Float64) || throw(ArgumentError(
        "a FEM circle-arc path cannot span more than one revolution"
    ))
    # Every consumer of a circular interface must use the same subdivision.
    # Span-relative subdivisions create overlapping curves when another
    # material traverses the same interface over a different angular range.
    _register_full_circle_breaks!(registry, centre, radius)
    _register_circle_break!(registry, centre, radius, start_angle)
    _register_circle_break!(registry, centre, radius, start_angle + value_span)
    circle_key = _circle_key(centre, radius)
    tolerance = FEM_BOUNDARY_ANGLE_TOLERANCE
    distances = Float64[0.0, distance]
    for angle in registry.circle_breaks[circle_key]
        offset = mod(direction * (angle - Float64(start_angle)), 2π)
        if tolerance < offset < distance - tolerance &&
           all(existing -> abs(existing - offset) > tolerance, distances)
            push!(distances, offset)
        end
    end
    sort!(distances)
    angles = Float64(start_angle) .+ direction .* distances
    points = [something(
                  _circle_break_point(registry, circle_key, angle),
                  _circle_point(centre, radius, angle)
              ) for angle in angles]
    first_point === nothing || (points[1] = first_point)
    last_point === nothing || (points[end] = last_point)
    centre_tag = _point!(registry, centre; mesh_size)
    point_tags = [_point!(registry, point; mesh_size) for point in points]
    return [begin
                first_tag = point_tags[index]
                last_tag = point_tags[index + 1]
                key = (:circle_arc, circle_key, minmax(first_tag, last_tag))
                tag, stored_first,
                stored_last = get!(registry.circle_arcs, key) do
                    (
                        gmsh.model.geo.add_circle_arc(
                            first_tag, centre_tag, last_tag
                        ),
                        first_tag,
                        last_tag
                    )
                end
                registry.curve_points[tag] = (stored_first, stored_last)
                get!(registry.curve_samples, tag) do
                    a = first_tag == stored_first ? angles[index] : angles[index + 1]
                    b = first_tag == stored_first ? angles[index + 1] : angles[index]
                    [_circle_point(centre, radius, a + fraction * (b - a))
                     for fraction in (0.0, 1e-6, ((1:15) ./ 16)..., 1 - 1e-6, 1.0)]
                end
                first_tag == stored_first && last_tag == stored_last ? tag : -tag
            end
            for index in 1:(length(point_tags) - 1)
            if point_tags[index] != point_tags[index + 1]]
end

function _sector_loop!(
        registry::FEMLoopRegistry,
        shape::DataModel.SectorShape;
        mesh_size = registry.mesh_size
)
    contacts = shape.contacts
    points = contacts.points
    transformed = map(
        point -> _transform_point(point, shape.at),
        (
            points.base_upper,
            points.base_lower,
            points.side_lower,
            points.back_lower,
            points.back_upper,
            points.side_upper
        )
    )
    key = (:sector,
        Tuple(
            (_coordinate_key(point[1]), _coordinate_key(point[2]))
        for point in transformed
        ))
    loop = get!(registry.loops, key) do
        point_tags = [_point!(registry, point; mesh_size) for point in transformed]
        arcs = contacts.arcs
        arc_specs = (
            (arcs.base, 1, 2),
            (arcs.lower, 3, 4),
            (arcs.back, 4, 5),
            (arcs.upper, 5, 6)
        )
        arc_curves = Dict{Symbol, Vector{Int}}()
        for (name, (arc, first_index, last_index)) in zip(
            (:base, :lower, :back, :upper), arc_specs
        )
            if iszero(arc.radius)
                arc_curves[name] = point_tags[first_index] == point_tags[last_index] ?
                                   Int[] :
                                   _line_path!(registry,
                    transformed[first_index], transformed[last_index]; mesh_size)
                continue
            end
            centre = _transform_point(arc.center, shape.at)
            arc_curves[name] = _circle_arc_path!(
                registry,
                centre,
                arc.radius,
                shape.at.φ + arc.start,
                arc.stop - arc.start;
                mesh_size,
                first_point = transformed[first_index],
                last_point = transformed[last_index]
            )
        end
        curves = Int[]
        append!(curves, arc_curves[:base])
        append!(curves, _line_path!(registry, transformed[2], transformed[3]; mesh_size))
        append!(curves, arc_curves[:lower])
        append!(curves, arc_curves[:back])
        append!(curves, arc_curves[:upper])
        append!(curves, _line_path!(registry, transformed[6], transformed[1]; mesh_size))
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        FEMLoop(ccw, cw, abs.(curves), curves)
    end
    return _refine_loop_mesh_size!(registry, loop, mesh_size)
end

function _bent_strip_loop!(
        registry::FEMLoopRegistry,
        shape::DataModel.BentStrip;
        mesh_size = registry.mesh_size
)
    centre = (shape.at.x, shape.at.y)
    if iszero(shape.ri) && isapprox(
            shape.span, 2π; rtol = 0, atol = FEM_BOUNDARY_ANGLE_TOLERANCE
    )
        return _circle_loop!(registry, centre, shape.ro; mesh_size)
    end
    start_angle = shape.at.φ - shape.span / 2
    stop_angle = shape.at.φ + shape.span / 2
    key = (
        :bent_strip,
        _coordinate_key(centre[1]),
        _coordinate_key(centre[2]),
        _coordinate_key(start_angle),
        _coordinate_key(shape.span),
        _coordinate_key(shape.ri),
        _coordinate_key(shape.ro)
    )
    loop = get!(registry.loops, key) do
        outer_start_point = _circle_point(centre, shape.ro, start_angle)
        outer_stop_point = _circle_point(centre, shape.ro, stop_angle)
        outer_start = _point!(registry, outer_start_point; mesh_size)
        outer_stop = _point!(registry, outer_stop_point; mesh_size)
        curves = _circle_arc_path!(
            registry,
            centre,
            shape.ro,
            start_angle,
            shape.span;
            mesh_size,
            first_point = outer_start_point,
            last_point = outer_stop_point
        )
        if iszero(shape.ri)
            centre_point = _point!(registry, centre; mesh_size)
            append!(curves, _line_path!(registry, outer_stop_point, centre; mesh_size))
            append!(curves, _line_path!(registry, centre, outer_start_point; mesh_size))
        else
            inner_start_point = _circle_point(centre, shape.ri, start_angle)
            inner_stop_point = _circle_point(centre, shape.ri, stop_angle)
            inner_start = _point!(registry, inner_start_point; mesh_size)
            inner_stop = _point!(registry, inner_stop_point; mesh_size)
            append!(curves, _line_path!(registry, outer_stop_point, inner_stop_point; mesh_size))
            append!(curves,
                _circle_arc_path!(
                    registry,
                    centre,
                    shape.ri,
                    stop_angle,
                    -shape.span;
                    mesh_size,
                    first_point = inner_stop_point,
                    last_point = inner_start_point
                ))
            append!(curves, _line_path!(registry, inner_start_point, outer_start_point; mesh_size))
        end
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        FEMLoop(ccw, cw, abs.(curves), curves)
    end
    return _refine_loop_mesh_size!(registry, loop, mesh_size)
end

function _shape_points(shape::DataModel.Rectangle)
    half_width = shape.w / 2
    half_height = shape.h / 2
    return [_transform_point(point, shape.at)
            for point in (
        (-half_width, -half_height),
        (half_width, -half_height),
        (half_width, half_height),
        (-half_width, half_height)
    )]
end

function _shape_points(shape::DataModel.Polygon)
    return [_transform_point(point, shape.at) for point in shape.points]
end

function _shape_points(shape::DataModel.SectorShape)
    return DataModel.tessellate(shape; points_per_arc = 24)
end

function _shape_points(shape::DataModel.EllipseOffset)
    return DataModel.tessellate(shape; points_per_arc = 96)
end

function _shape_points(shape::DataModel.BentStrip)
    return DataModel.tessellate(shape; points_per_arc = 24)
end

function _shape_points(shape)
    _fem_error(
        :unsupported,
        string(typeof(shape)),
        :primitive,
        "the resolved shape has no built-in geo-kernel boundary adaptation"
    )
end

function _boundary_loop!(
        registry::FEMLoopRegistry,
        shape;
        mesh_size = registry.mesh_size
)
    if shape isa DataModel.Disk
        return _circle_loop!(
            registry,
            (shape.at.x, shape.at.y),
            shape.r;
            mesh_size
        )
    end
    shape isa DataModel.Ellipse && return _ellipse_loop!(
        registry, shape; mesh_size
    )
    shape isa DataModel.SectorShape && return _sector_loop!(
        registry, shape; mesh_size
    )
    shape isa DataModel.BentStrip && return _bent_strip_loop!(
        registry, shape; mesh_size
    )
    return _polygon_loop!(registry, _shape_points(shape); mesh_size)
end

function _surface_boundaries(shape::DataModel.Annulus)
    outer = DataModel.Disk(shape.ro, shape.at)
    inner = DataModel.Disk(shape.ri, shape.at)
    return outer, Any[inner]
end

_surface_boundaries(shape::DataModel.ShellShape) = (shape.outer, Any[shape.inner])

function _surface_boundaries(shape::DataModel.DifferenceShape)
    outer, outer_holes = _surface_boundaries(shape.outer)
    return outer, Any[outer_holes..., shape.holes...]
end

_surface_boundaries(shape) = (shape, Any[])

function _surface_faces!(registry::FEMLoopRegistry, shape, mesh_size)
    outer, holes = _surface_boundaries(shape)
    outer_loop = _boundary_loop!(registry, outer; mesh_size)
    hole_loops = [_boundary_loop!(registry, hole; mesh_size)
                  for hole in holes]
    isempty(holes) && return Int[gmsh.model.geo.add_plane_surface([outer_loop.ccw])]
    seen = Set{Int}()
    touching = false
    for loop in Iterators.flatten(((outer_loop,), hole_loops))
        points = Set(point for curve in loop.curves for point in registry.curve_points[curve])
        touching |= any(point -> point in seen, points)
        union!(seen, points)
    end
    if !touching
        return Int[gmsh.model.geo.add_plane_surface([outer_loop.ccw; getproperty.(hole_loops, :cw)])]
    end
    return _material_faces!(registry,
        (outer_loop.oriented, (-loop.oriented for loop in hole_loops)...))
end

function _material_faces!(registry::FEMLoopRegistry, boundaries)
    # Every oriented edge has material on its left. Cancel shared hole edges:
    # a metal-metal seam cannot also bound a third (filler) material.
    counts = Dict{Int, Int}()
    for curves in boundaries
        for curve in curves
            counts[abs(curve)] = get(counts, abs(curve), 0) + sign(curve)
        end
    end
    all(value -> abs(value) <= 1, Base.values(counts)) || throw(ArgumentError(
        "overlapping boundaries in a FEM material face"))
    edges = sort!([tag * count for (tag, count) in counts if !iszero(count)]; by=abs)
    positions = Dict(tag => point for (point, tag) in registry.points)
    samples = Dict{Int, Vector{Tuple{Float64, Float64}}}()
    outgoing = Dict{Int, Vector{Int}}()
    endpoints = Dict{Int, Tuple{Int, Int}}()
    for edge in edges
        first_point, last_point = registry.curve_points[abs(edge)]
        values = copy(get(registry.curve_samples, abs(edge),
            [positions[first_point], positions[last_point]]))
        values[1], values[end] = positions[first_point], positions[last_point]
        if edge < 0
            first_point, last_point = last_point, first_point
            reverse!(values)
        end
        endpoints[edge] = (first_point, last_point)
        samples[edge] = values
        push!(get!(Vector{Int}, outgoing, first_point), edge)
    end
    # Compare incident curves at a common physical distance from the contact.
    # Equal parameter fractions on tangent circles of different radii can give
    # identical directions and the wrong face pairing.
    steps = Dict{Int, Float64}()
    for (first_point, last_point) in values(endpoints)
        a, b = positions[first_point], positions[last_point]
        step = hypot(b[1]-a[1], b[2]-a[2]) / 64
        for point in (first_point, last_point)
            steps[point] = min(get(steps, point, Inf), step)
        end
    end
    circles = Dict{Int, Tuple{Float64, Float64, Float64, Int}}()
    for (key, (tag, first_point, last_point)) in registry.circle_arcs
        cx, cy, radius = key[2]
        a, b = positions[first_point], positions[last_point]
        orientation = sign((a[1]-cx)*(b[2]-cy)-(a[2]-cy)*(b[1]-cx))
        circles[tag] = (cx, cy, radius, Int(orientation))
    end
    function direction(edge, step)
        if haskey(circles, abs(edge))
            cx, cy, radius, orientation = circles[abs(edge)]
            first_point, last_point = registry.curve_points[abs(edge)]
            point = positions[edge > 0 ? first_point : last_point]
            turn = sign(edge) * orientation
            return atan(turn*(point[1]-cx), -turn*(point[2]-cy)) + turn*step/(2radius)
        end
        points = haskey(samples, edge) ? samples[edge] : reverse(samples[-edge])
        return atan(points[2][2]-points[1][2], points[2][1]-points[1][1])
    end
    successors = Dict{Int, Int}()
    for edge in edges
        step = steps[endpoints[edge][2]]
        backwards = direction(-edge, step)
        candidates = get(outgoing, endpoints[edge][2], Int[])
        isempty(candidates) && throw(ArgumentError("open FEM material boundary"))
        successors[edge] = candidates[argmin(map(candidates) do candidate
            mod(backwards - direction(candidate, step), 2π)
        end)]
    end
    length(unique(Base.values(successors))) == length(edges) || throw(ArgumentError(
        "ambiguous FEM material boundary at a contact"))

    cycles = Vector{Int}[]
    polygons = Vector{Tuple{Float64, Float64}}[]
    remaining = Set(edges)
    for first_edge in edges
        first_edge in remaining || continue
        cycle = Int[]
        polygon = Tuple{Float64, Float64}[]
        edge = first_edge
        while edge in remaining
            delete!(remaining, edge)
            push!(cycle, edge)
            append!(polygon, samples[edge][1:end-1])
            edge = successors[edge]
        end
        edge == first_edge || throw(ArgumentError("nonclosing FEM material face"))
        push!(cycles, cycle)
        push!(polygons, polygon)
    end
    areas = _signed_area.(polygons)
    outers = findall(>(0), areas)
    isempty(outers) && throw(ArgumentError("FEM material has no positive-area face"))
    interiors = Dict(index => Int[] for index in outers)
    function contains(polygon, point)
        inside = false
        for index in eachindex(polygon)
            a, b = polygon[index], polygon[mod1(index + 1, end)]
            (a[2] > point[2]) == (b[2] > point[2]) && continue
            point[1] < a[1] + (point[2] - a[2]) * (b[1] - a[1]) / (b[2] - a[2]) &&
                (inside = !inside)
        end
        return inside
    end
    for index in findall(<(0), areas)
        parents = filter(parent -> contains(polygons[parent], first(polygons[index])), outers)
        isempty(parents) && throw(ArgumentError("uncontained FEM material hole"))
        parent = parents[argmin(areas[parents])]
        push!(interiors[parent], index)
    end
    loops = [gmsh.model.geo.add_curve_loop(cycle) for cycle in cycles]
    return [gmsh.model.geo.add_plane_surface([loops[index]; loops[interiors[index]]])
            for index in outers]
end

function _tangent_fill_surfaces!(registry::FEMLoopRegistry, shape, mesh_size)
    plan = _tangent_circular_fill_plan(shape)
    plan === nothing && return nothing
    surfaces = Int[]
    count = length(plan.holes)
    for index in eachindex(plan.holes)
        next_index = mod1(index + 1, count)
        angle = plan.angles[index]
        next_angle = plan.angles[next_index]
        span = mod(next_angle - angle, 2π)
        span <= FEM_BOUNDARY_ANGLE_TOLERANCE && (span = 2π)
        hole = plan.holes[index]
        next_hole = plan.holes[next_index]
        hole_centre = (Float64(hole.at.x), Float64(hole.at.y))
        next_hole_centre = (
            Float64(next_hole.at.x), Float64(next_hole.at.y)
        )
        outer_point = _radial_point(plan.centre, plan.outer_radius, angle)
        next_outer_point = _radial_point(
            plan.centre, plan.outer_radius, next_angle
        )
        inner_point = _radial_point(plan.centre, plan.inner_radius, angle)
        next_inner_point = _radial_point(
            plan.centre, plan.inner_radius, next_angle
        )
        curves = Int[]
        append!(curves,
            _circle_arc_path!(
                registry,
                plan.centre,
                plan.outer_radius,
                angle,
                span;
                mesh_size,
                first_point = outer_point,
                last_point = next_outer_point
            ))
        append!(curves,
            _circle_arc_path!(
                registry,
                next_hole_centre,
                Float64(next_hole.r),
                next_angle,
                -π;
                mesh_size,
                first_point = next_outer_point,
                last_point = next_inner_point
            ))
        append!(curves,
            _circle_arc_path!(
                registry,
                plan.centre,
                plan.inner_radius,
                next_angle,
                -span;
                mesh_size,
                first_point = next_inner_point,
                last_point = inner_point
            ))
        append!(curves,
            _circle_arc_path!(
                registry,
                hole_centre,
                Float64(hole.r),
                angle + π,
                -π;
                mesh_size,
                first_point = inner_point,
                last_point = outer_point
            ))
        # Neighbouring wires may also touch each other. The annular sector is
        # then a pinched walk, not one face: trace each material-left component
        # with the same contact rule used by every other enclosure boundary.
        append!(surfaces, _material_faces!(registry, (curves,)))
    end
    for (start_angle, span) in plan.outer_spans
        push!(surfaces,
            _annular_sector_surface!(
                registry,
                plan.centre,
                plan.outer_radius,
                plan.shell_outer_radius,
                start_angle,
                span,
                mesh_size
            ))
    end
    return surfaces
end

function _surfaces!(registry::FEMLoopRegistry, shape, mesh_size)
    compartments = _tangent_fill_surfaces!(registry, shape, mesh_size)
    compartments === nothing || return compartments
    return _surface_faces!(registry, shape, mesh_size)
end

_boundary_components(shape::DataModel.AssemblyShape) = collect(shape.members)
_boundary_components(shape) = Any[shape]

function _entity_boundary(surfaces)
    isempty(surfaces) && return Int[]
    dimtags = [(2, surface) for surface in surfaces]
    return sort!(unique(Int(tag)
    for (dim, tag) in gmsh.model.get_boundary(dimtags, true, false, false) if dim == 1))
end

_lowest_shape_point(shape) = argmin(p -> (p[2],p[1]), _shape_points(shape))
_lowest_shape_point(shape::DataModel.Disk) = (shape.at.x,shape.at.y-shape.r)
_lowest_shape_point(shape::DataModel.Annulus) = (shape.at.x,shape.at.y-shape.ro)
_lowest_shape_point(shape::Union{DataModel.ShellShape,DataModel.DifferenceShape}) =
    _lowest_shape_point(shape.outer)
_lowest_shape_point(shape::DataModel.AssemblyShape) =
    argmin(p -> (p[2],p[1]), _lowest_shape_point.(shape.members))
function _lowest_shape_point(shape::DataModel.Ellipse)
    return argmin(p -> (p[2],p[1]), [_transform_point(p,shape.at)
        for p in ((shape.a,0.),(0.,shape.b),(-shape.a,0.),(0.,-shape.b))])
end
function _lowest_shape_point(shape::DataModel.SectorShape)
    points = [_transform_point(p,shape.at) for p in values(shape.contacts.points)]
    for arc in values(shape.contacts.arcs)
        iszero(arc.radius) && continue
        start, span = shape.at.φ+arc.start, arc.stop-arc.start
        if mod(3π/2-start,2π) <= span
            centre = _transform_point(arc.center,shape.at)
            push!(points,(centre[1],centre[2]-arc.radius))
        end
    end
    return argmin(p -> (p[2],p[1]),points)
end
function _lowest_shape_point(shape::DataModel.BentStrip)
    centre = (shape.at.x,shape.at.y)
    start, stop = shape.at.φ-shape.span/2, shape.at.φ+shape.span/2
    angles = mod(3π/2-start,2π) <= shape.span ? (start,stop,3π/2) : (start,stop)
    return argmin(p -> (p[2],p[1]), [_circle_point(centre,r,a)
        for r in (shape.ri,shape.ro) for a in angles])
end

# Relative circle-area error is bounded by 2pi^2/(3N^2). Dyadic multiples
# of twelve reproduce the qualified 96/192-point levels without measuring
# a mesh or refining in response to its result.
function _conductor_curve_geometry(shape, curve)
    lower, upper = gmsh.model.get_parametrization_bounds(1, curve)
    a = gmsh.model.get_value(1, curve, lower)[1:2] .- (shape.at.x, shape.at.y)
    b = gmsh.model.get_value(1, curve, upper)[1:2] .- (shape.at.x, shape.at.y)
    angle = atan(abs(a[1]*b[2]-a[2]*b[1]), a[1]*b[1]+a[2]*b[2])
    return (; fraction=angle/(2π), length=angle*hypot(a...))
end

function _sector_partition!(registry, shape::DataModel.SectorShape, mesh_size)
    outer = _boundary_loop!(registry, shape; mesh_size)
    centre = DataModel.centroid(shape)
    # Retain the physical contour. A 0.1-scale internal copy gives one native
    # four-sided strip per contour segment and an unstructured central patch.
    coordinates = Dict(tag => point for (point,tag) in registry.points)
    inner_point = Dict{Int,Int}()
    spokes = Dict{Int,Int}()
    depth = 0.0
    for curve in outer.curves, p in registry.curve_points[curve]
        haskey(inner_point,p) && continue
        point = coordinates[p]
        target = centre .+ 0.1 .* (point .- centre)
        q = _point!(registry,target;mesh_size)
        inner_point[p] = q
        spokes[p] = _line!(registry,p,q)
        depth = max(depth,0.9hypot((point .- centre)...))
    end
    circles = Dict(tag => key[2] for (key,(tag,_,_)) in registry.circle_arcs)
    inner_curves = Int[]
    curve_pairs = Tuple{Int,Int,Float64,Float64}[]
    patches = Dict{Int,NTuple{4,Int}}()
    for curve in outer.oriented
        a,b = registry.curve_points[abs(curve)]
        p,q = curve > 0 ? (a,b) : (b,a)
        ip,iq = inner_point[p],inner_point[q]
        if haskey(circles,abs(curve))
            x,y,radius = circles[abs(curve)]
            copied_centre = centre .+ 0.1 .* ((x,y) .- centre)
            c = _point!(registry,copied_centre;mesh_size)
            inner = gmsh.model.geo.add_circle_arc(ip,c,iq)
            u,v = coordinates[p] .- (x,y), coordinates[q] .- (x,y)
            turn = atan(abs(u[1]*v[2]-u[2]*v[1]),u[1]*v[1]+u[2]*v[2])
            length = radius*turn
        else
            inner = _line!(registry,ip,iq)
            length,turn = hypot((coordinates[q] .- coordinates[p])...),0.0
        end
        push!(inner_curves,inner)
        push!(curve_pairs,(abs(curve),abs(inner),length,turn))
        loop = gmsh.model.geo.add_curve_loop([curve,spokes[q],-inner,-spokes[p]])
        surface = gmsh.model.geo.add_plane_surface([loop])
        patches[surface] = (p,q,iq,ip)
        gmsh.model.geo.mesh.set_transfinite_surface(surface,"AlternateLeft",[p,q,iq,ip])
    end
    core = gmsh.model.geo.add_plane_surface([gmsh.model.geo.add_curve_loop(inner_curves)])
    partition = FEMSectorPartition(core,curve_pairs,collect(values(spokes)),depth,patches)
    return [sort!(collect(keys(patches)));core],partition
end

# Build the frequency-independent physical CAD once. Native exterior and mesh
# sources consume this topology; material and terminal ownership stays with the CAD.
function _build_physical_geometry!(
        model::FEMResolvedModel,
        model_name::String
)
    gmsh.model.add(model_name)
    registry = FEMLoopRegistry(model.cad_scale)
    for region in model.region_plans
        _register_shape_breaks!(registry, region.shape)
    end
    foreach(
        boundary -> _register_shape_breaks!(registry, boundary),
        model.cable_boundaries
    )
    _register_circle_contacts!(registry)
    material_surfaces = [Int[] for _ in model.material_plans]
    region_surfaces = Vector{Int}[]
    sector_partitions = Dict{Int,FEMSectorPartition}()
    terminal_surfaces = [Int[] for _ in model.terminal_ids]
    for (index,region) in enumerate(model.region_plans)
        surfaces = if region.shape isa DataModel.SectorShape &&
                      model.material_plans[region.material_index].kind === :conductor
            members,partition = _sector_partition!(registry,region.shape,model.cad_scale)
            sector_partitions[index] = partition
            members
        else
            _surfaces!(registry, region.shape, model.cad_scale)
        end
        push!(region_surfaces, surfaces)
        append!(material_surfaces[region.material_index], surfaces)
        region.terminal_index > 0 && append!(
            terminal_surfaces[region.terminal_index], surfaces
        )
    end
    cable_curves = [Int[] for _ in model.cable_boundaries]
    cable_loops = [Int[] for _ in model.cable_boundaries]
    cable_loop_curves = [Vector{Int}[] for _ in model.cable_boundaries]
    for cable_index in eachindex(model.cable_boundaries)
        for component in _boundary_components(model.cable_boundaries[cable_index])
            loop = _boundary_loop!(registry, component;
                mesh_size = model.cad_scale)
            push!(cable_loops[cable_index], loop.cw)
            push!(cable_loop_curves[cable_index], -reverse(loop.oriented))
            append!(cable_curves[cable_index], loop.curves)
        end
    end
    _apply_point_mesh_sizes!(registry)
    gmsh.model.geo.synchronize()
    # Interior anchors are geometry facts, not reference vertices. Retain a
    # horizontal clearance of 1% of the metal dimension for the native shift.
    anchors = Tuple{Float64,Float64}[]
    metal_dimensions = Float64[]
    for terminal in eachindex(model.terminal_ids)
        indices = findall(r -> r.terminal_index == terminal,model.region_plans)
        dimensions = map(indices) do i
            shape = model.region_plans[i].shape
            section = _conductor_section(shape)
            section === nothing ? minimum(gmsh.model.get_bounding_box(2,first(region_surfaces[i]))[k+3]-gmsh.model.get_bounding_box(2,first(region_surfaces[i]))[k] for k in 1:2) :
                section.kind == 1 ? 2section.width : section.width
        end
        dimension = minimum(dimensions)
        anchor = nothing
        for (i,width) in zip(indices,dimensions)
            low = _lowest_shape_point(model.region_plans[i].shape)
            for dy in (.125,.25,.5,.75), dx in (0.,.125,-.125,.25,-.25)
                point = (low[1]+dx*width,low[2]+dy*width)
                if any(region_surfaces[i]) do surface
                        all(shift -> gmsh.model.is_inside(2,surface,[point[1]+shift*dimension,point[2],0.])>0,(0.,.01))
                    end
                    anchor = point
                    break
                end
            end
            anchor === nothing || break
        end
        anchor === nothing && _fem_error(:geometry,model.problem.system.system_id,:measurement_line,"cannot locate an interior metal anchor for terminal $terminal")
        push!(anchors,anchor); push!(metal_dimensions,dimension)
    end

    return (; anchors, metal_dimensions, material_surfaces, region_surfaces,
        sector_partitions, terminal_surfaces, cable_curves, cable_loops, cable_loop_curves)
end

# Serialize intrinsic source dimensions and CAD ownership. Native parameters
# apply the public numerical controls; no frequency-dependent target is written.
function _write_native_geometry_data(io, model)
    system = model.problem.system
    println(io,"CableX() = ",_pro_array([p.x for p in system.positions]),";")
    println(io,"CableY() = ",_pro_array([p.y for p in system.positions]),";")
    println(io,"CableRadius() = ",_pro_array(LineCableModels.outer_radius.(system.designs)),";")
    println(io,"NumPhysicalRegions = ",length(model.region_plans),";")
    offset = 0
    sources = Vector{Vector{Any}}(undef,length(model.region_plans))
    for (cable,design) in pairs(system.designs)
        count = length(design.geometry.regions)
        regions = @view system.geometry[offset+1:offset+count]
        terminals = @view system.terminal_map[offset+1:offset+count]
        formations = _formations(regions,terminals,"cable_$cable")
        for (index,plan) in pairs(model.region_plans)
            plan.cable_index == cable || continue
            formation = findfirst(f -> f.complete && first(f.members) == plan.region_index,formations)
            members = formation === nothing ? [plan.region_index] : formations[formation].members
            sources[index] = Any[regions[member] for member in members]
        end
        offset += count
    end
    for (index,members) in pairs(sources)
        dimensions = map(members) do region
            shape = region.source.primitive
            repeated = any(p -> p.owner === DataModel.Group,region.placement.patterns)
            if shape isa DataModel.Disk
                (1,shape.r,0.,repeated)
            elseif shape isa DataModel.Annulus
                (2,shape.ro,shape.ri,repeated)
            elseif shape isa DataModel.Rectangle
                (3,shape.w,shape.h,repeated)
            elseif shape isa DataModel.Ellipse
                (4,shape.a,shape.b,repeated)
            elseif region.primitive isa DataModel.Annulus
                placed = region.primitive
                (2,placed.ro,placed.ri,repeated)
            else
                (0,LineCableModels.area(region.primitive),DataModel.perimeter(region.primitive),repeated)
            end
        end
        for (name,column) in (("RegionDimensionKinds",1),("RegionDimension1",2),("RegionDimension2",3),("RegionRepeated",4))
            println(io,name,"~{",index-1,"}() = ",_pro_array(getindex.(dimensions,column)),";")
        end
    end
    for cable in eachindex(system.designs)
        indices = findall(p -> p.cable_index == cable,model.region_plans)
        sizing = model.cable_boundaries[cable] isa DataModel.AssemblyShape ? indices : [last(indices)]
        println(io,"CableSizingRegions~{",cable-1,"}() = ",_pro_array(sizing.-1),";")
    end
end

function _write_physical_geometry(path, model, physical)
    gmsh.model.geo.synchronize()
    (;material_surfaces,region_surfaces,terminal_surfaces,cable_curves,cable_loops,
        sector_partitions,anchors,metal_dimensions,cable_loop_curves) = physical
    raw = path*"_unrolled"
    gmsh.write(raw)
    open(path,"w") do io
        println(io,"// Frequency-independent physical CAD; numerical sizing is native.")
        for line in eachline(raw)
            startswith(line,"Transfinite") && continue
            startswith(line,"Recombine") && continue
            occursin(r"^cl__\d+ =",line) && continue
            point = match(r"^Point\((\d+)\)",line)
            if point !== nothing
                tag = parse(Int,point[1])
                xyz = gmsh.model.get_value(0,tag,Float64[])
                println(io,"Point(",tag,") = {",join(_pro_number.(xyz),", "),", MeshFine};")
            else
                println(io,line)
            end
        end
        emit(name,values) = println(io,name,"() = ",_pro_array(values),";")
        emit("FEMMaterialSurfaces",reduce(vcat,material_surfaces;init=Int[]))
        emit("FEMCableCurves",unique(reduce(vcat,cable_curves;init=Int[])))
        println(io,"FEMAirHoles() = {}; FEMEarthHoles() = {};")
        for cable in eachindex(cable_curves)
            for curves in cable_loop_curves[cable]
                println(io,"FEMHole = newll; Curve Loop(FEMHole) = {",join(curves,","),"};")
                println(io,model.cable_hosts[cable] === :air ? "FEMAirHoles" : "FEMEarthHoles","() += {FEMHole};")
            end
            emit("FEMCableCurves~{$(cable-1)}",cable_curves[cable])
            passive_regions = [i for (i,r) in pairs(model.region_plans) if r.cable_index == cable && model.material_plans[r.material_index].kind !== :conductor]
            passive_surfaces = reduce(vcat,(region_surfaces[i] for i in passive_regions);init=Int[])
            emit("FEMCablePassiveSurfaces~{$(cable-1)}",passive_surfaces)
            emit("FEMCablePassiveCurves~{$(cable-1)}",_entity_boundary(passive_surfaces))
            _write_native_physical_group(io,1,model.tags.cable_contour_base+cable,
                @sprintf("LCM/cable_contour/%04d/%s",cable,model.problem.system.designs[cable].cable_id),cable_curves[cable])
        end
        for (i,surfaces) in pairs(region_surfaces)
            curves = _entity_boundary(surfaces)
            tags = unique(last.(gmsh.model.get_boundary([(1,c) for c in curves],false,false,false)))
            emit("FEMRegionSurfaces~{$(i-1)}",surfaces)
            emit("FEMRegionCurves~{$(i-1)}",curves)
            emit("FEMRegionPoints~{$(i-1)}",tags)
            region = model.region_plans[i]
            section = _conductor_section(region.shape)
            isconductor = section !== nothing && model.material_plans[region.material_index].kind === :conductor
            println(io,"FEMConductor~{",i-1,"} = ",Int(isconductor),";")
            isconductor || continue
            println(io,"FEMConductorMaterial~{",i-1,"} = ",region.material_index,";")
            println(io,"FEMConductorWidth~{",i-1,"} = ",_pro_number(section.width),";")
            println(io,"FEMConductorKind~{",i-1,"} = ",section.kind,";")
            partition = get(sector_partitions,i,nothing)
            println(io,"FEMConductorSector~{",i-1,"} = ",Int(partition !== nothing),";")
            if partition === nothing
                arcs = [_conductor_curve_geometry(region.shape,curve) for curve in curves]
                emit("FEMConductorArcLengths~{$(i-1)}",getproperty.(arcs,:length))
                emit("FEMConductorArcFractions~{$(i-1)}",getproperty.(arcs,:fraction))
                emit("FEMConductorBulkSurfaces~{$(i-1)}",surfaces)
            else
                println(io,"FEMSectorDepth~{",i-1,"} = ",_pro_number(partition.depth),";")
                emit("FEMSectorOuter~{$(i-1)}",first.(partition.curve_pairs))
                emit("FEMSectorInner~{$(i-1)}",getindex.(partition.curve_pairs,2))
                emit("FEMSectorLengths~{$(i-1)}",getindex.(partition.curve_pairs,3))
                emit("FEMSectorTurns~{$(i-1)}",getindex.(partition.curve_pairs,4))
                emit("FEMSectorSpokes~{$(i-1)}",partition.spokes)
                emit("FEMConductorBulkSurfaces~{$(i-1)}",[partition.core])
            end
        end
        for (i,material) in pairs(model.material_plans)
            _write_native_physical_group(io,2,material.physical_tag,material.physical_name,material_surfaces[i])
        end
        conductor = Int[]; passive = Int[]
        for (i,material) in pairs(model.material_plans)
            append!(material.kind === :conductor ? conductor : passive,material_surfaces[i])
        end
        _write_native_physical_group(io,2,6001,"LCM/domain/conductors",conductor)
        isempty(passive) || _write_native_physical_group(io,2,6002,"LCM/domain/passive_cable_media",passive)
        for (i,surfaces) in pairs(terminal_surfaces)
            _write_native_physical_group(io,2,model.tags.terminal_base+i,model.terminal_names[i],surfaces)
            _write_native_physical_group(io,1,model.tags.terminal_contour_base+i,
                @sprintf("LCM/terminal_contour/%04d",i),_entity_boundary(surfaces))
        end
        emit("FEMReceiverX",first.(anchors))
        emit("FEMReceiverY",last.(anchors))
        emit("FEMReceiverMetalDimension",metal_dimensions)
        for partition in values(sector_partitions), (surface,corners) in partition.patches
            println(io,"Transfinite Surface {",surface,"} = {",join(corners,","),"} AlternateLeft;")
        end
    end
    rm(raw)
    return nothing
end

function _write_native_physical_group(io, dim, tag, name, entities)
    println(io,"Physical ",("Point","Curve","Surface")[dim+1],"(",_pro_string(name),",",tag,") = {",join(entities,","),"};")
end

function _write_native_mesh_entry(path, data, physical, assets)
    write(path,"""
    Include "$data";
    Include "$assets/parameters.pro";
    Include "$physical";
    Include "$assets/geometry.geo";
    Include "$assets/mesh.geo";
    """)
    return path
end
