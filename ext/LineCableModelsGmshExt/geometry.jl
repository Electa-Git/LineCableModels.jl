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
    voltage_abscissae::Vector{Float64}
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
        Float64(mesh_size),
        Float64[]
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
    # Split CAD interfaces at the measurement abscissae before making surfaces.
    # Gmsh can then embed each path interval using the shared interface vertices.
    if !iszero(dx)
        for x in registry.voltage_abscissae
            position = (x - first[1]) / dx
            tolerance / distance < position < 1 - tolerance / distance || continue
            _point!(registry, (x, first[2] + position * dy); mesh_size)
        end
    end
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
        ax, bx = shape.a*cos(shape.at.φ), -shape.b*sin(shape.at.φ)
        phase, amplitude = atan(bx, ax), hypot(ax, bx)
        for x in registry.voltage_abscissae
            offset = (x - shape.at.x)/amplitude
            abs(offset) < 1 || continue
            for angle in mod.((phase-acos(offset),phase+acos(offset)),2π)
                all(a -> abs(a-angle)>FEM_BOUNDARY_ANGLE_TOLERANCE,angles) && push!(angles,angle)
            end
        end
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
    for x in registry.voltage_abscissae
        offset = (x - centre[1]) / radius
        abs(offset) < 1 || continue
        angle = acos(offset)
        for a in (angle, -angle)
            _register_circle_break!(registry, centre, radius, a;
                point=(x, centre[2] + radius * sin(a)))
        end
    end
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

function _physical_group(dim::Int, entities, tag::Int, name::String)
    values = sort!(unique(Int.(entities)))
    isempty(values) && return nothing
    gmsh.model.add_physical_group(dim, values, tag)
    gmsh.model.set_physical_name(dim, tag, name)
    return tag
end

function _entity_boundary(surfaces)
    isempty(surfaces) && return Int[]
    dimtags = [(2, surface) for surface in surfaces]
    return sort!(unique(Int(tag)
    for (dim, tag) in gmsh.model.get_boundary(dimtags, true, false, false) if dim == 1))
end

function _validate_material_interfaces!(model::FEMResolvedModel, material_surfaces)
    surfaces = collect(Iterators.flatten(material_surfaces))
    length(unique(surfaces)) == length(surfaces) || _fem_error(
        :geometry,
        model.problem.system.system_id,
        :material_partition,
        "a FEM surface is assigned to more than one material"
    )
    curves = Int[]
    for surfaces in material_surfaces
        append!(curves, _entity_boundary(surfaces))
    end
    offending = Int[]
    for curve in sort!(unique(curves))
        adjacent, _ = gmsh.model.get_adjacencies(1, curve)
        length(unique(adjacent)) == 2 || push!(offending, curve)
    end
    isempty(offending) || _fem_error(
        :geometry,
        model.problem.system.system_id,
        :material_partition,
        "material interfaces must have exactly two adjacent surfaces; " *
        "invalid curves: $(join(offending, ", "))"
    )
    return nothing
end

function _interface_mesh_sizes(model::FEMResolvedModel, mesh_plan::FEMMeshPlan)
    centre_x = model.centre[1]
    sizes = Dict{Float64, Float64}()
    function register(x, size)
        key = _coordinate_key(x)
        sizes[key] = min(get(sizes, key, Inf), Float64(size))
    end
    register(centre_x, mesh_plan.interface_mesh_size)
    for offset in (-2.0, 2.0)
        register(centre_x + offset, mesh_plan.domain_mesh_size)
    end
    for (index, position) in enumerate(model.problem.system.positions)
        register(position.x, mesh_plan.cable_interface_mesh_sizes[index])
    end
    return sizes
end

# Geometric intervals from the interface towards the far boundary. The
# logarithmic mean chooses a count that bounds both endpoint element sizes.
function _exterior_edge_grading(length, first_size, last_size)
    first_size = min(first_size, last_size)
    first_size == last_size && return (max(2, ceil(Int, length / last_size)), 1.0)
    relative_gap = (last_size - first_size) / first_size
    # Avoid losing the size difference when the endpoint targets are close.
    log_ratio = relative_gap < 0.5 ? log1p(relative_gap) : log(last_size / first_size)
    count = max(2, ceil(Int, length * log_ratio / (last_size - first_size)))
    return count, exp(log_ratio / (count - 1))
end

# Native geometric progression, with representable nodes from interface to wall.
# Evaluate the normalized exponential without overflowing exp(g). The rounded
# native ratio defines the check: sufficiently small exponents round to one.
function _pml_progression(count, grading, inner, outer)
    ratio = count == 1 ? 1.0 : exp(grading / count)
    isfinite(ratio) && isfinite(inner) && isfinite(outer) && inner != outer ||
        throw(ArgumentError("PML progression or endpoints are not representable"))
    exponent = count * log(ratio)
    previous = inner
    direction = sign(outer - inner)
    for i in 1:count-1
        fraction = ratio == 1 ? i / count :
            exp(exponent * (i / count - 1)) *
            (-expm1(-exponent * i / count)) / (-expm1(-exponent))
        coordinate = inner + (outer - inner) * fraction
        isfinite(coordinate) && direction * (coordinate - previous) > 0 &&
            direction * (outer - coordinate) > 0 || throw(ArgumentError(
                "PML spacing is not representable for count=$count, grading=$grading, endpoints=($inner, $outer)"))
        previous = coordinate
    end
    return ratio
end

# The same CAD vertices used by the primitive constructors select the voltage
# endpoint. This is geometry bookkeeping, independent of mesh nodes or fields.
_voltage_endpoint(shape) = argmin(p -> (p[2],p[1]), _shape_points(shape))
_voltage_endpoint(shape::DataModel.Disk) = (shape.at.x,shape.at.y-shape.r)
_voltage_endpoint(shape::DataModel.Annulus) = (shape.at.x,shape.at.y-shape.ro)
_voltage_endpoint(shape::Union{DataModel.ShellShape,DataModel.DifferenceShape}) =
    _voltage_endpoint(shape.outer)
_voltage_endpoint(shape::DataModel.AssemblyShape) =
    argmin(p -> (p[2],p[1]), _voltage_endpoint.(shape.members))
function _voltage_endpoint(shape::DataModel.Ellipse)
    return argmin(p -> (p[2],p[1]), [_transform_point(p,shape.at)
        for p in ((shape.a,0.),(0.,shape.b),(-shape.a,0.),(0.,-shape.b))])
end
function _voltage_endpoint(shape::DataModel.SectorShape)
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
function _voltage_endpoint(shape::DataModel.BentStrip)
    centre = (shape.at.x,shape.at.y)
    start, stop = shape.at.φ-shape.span/2, shape.at.φ+shape.span/2
    angles = mod(3π/2-start,2π) <= shape.span ? (start,stop,3π/2) : (start,stop)
    return argmin(p -> (p[2],p[1]), [_circle_point(centre,r,a)
        for r in (shape.ri,shape.ro) for a in angles])
end

# Relative circle-area error is bounded by 2pi^2/(3N^2). Dyadic multiples
# of twelve reproduce the qualified 96/192-point levels without measuring
# a mesh or refining in response to its result.
_conductor_circle_segments(tolerance) = 12 * 2^max(0,
    ceil(Int, log2(sqrt(2π^2/(3tolerance))/12)))

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

function _build_geometry!(
        model::FEMResolvedModel,
        model_name::String,
        mesh_plan::FEMMeshPlan = last(model.mesh_plans);
        reuse::Union{Nothing, FEMGeometry} = nothing
)
    endpoints = [argmin(p -> (p[2],p[1]),
        [_voltage_endpoint(region.shape) for region in model.region_plans
         if region.terminal_index == i]) for i in eachindex(model.terminal_ids)]
    abscissae = sort!(unique(_coordinate_key(p[1]) for p in endpoints))
    if reuse === nothing
        gmsh.model.add(model_name)
        registry = FEMLoopRegistry(model.fine_mesh_size)
        append!(registry.voltage_abscissae,abscissae)
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
                members,partition = _sector_partition!(registry,region.shape,region.mesh_size)
                sector_partitions[index] = partition
                members
            else
                _surfaces!(registry, region.shape, region.mesh_size)
            end
            push!(region_surfaces, surfaces)
            append!(material_surfaces[region.material_index], surfaces)
            region.terminal_index > 0 && append!(
                terminal_surfaces[region.terminal_index], surfaces
            )
        end
        cable_curves = [Int[] for _ in model.cable_boundaries]
        cable_loops = [Int[] for _ in model.cable_boundaries]
        for cable_index in eachindex(model.cable_boundaries)
            for component in _boundary_components(model.cable_boundaries[cable_index])
                loop = _boundary_loop!(registry, component;
                    mesh_size = model.cable_outer_mesh_sizes[cable_index])
                push!(cable_loops[cable_index], loop.cw)
                append!(cable_curves[cable_index], loop.curves)
            end
        end
        _apply_point_mesh_sizes!(registry)
    else
        gmsh.model.set_current(model_name)
        material_surfaces = reuse.material_surfaces
        region_surfaces = reuse.region_surfaces
        sector_partitions = reuse.sector_partitions
        terminal_surfaces = reuse.terminal_surfaces
        cable_curves = reuse.cable_curves
        cable_loops = reuse.cable_loops
    end

    # Keep exterior bookkeeping separate so frequency changes never recreate
    # or transform the authoritative cable geometry.
    gmsh.model.geo.synchronize()
    material_curves = unique(abs.(last.(gmsh.model.get_boundary(
        [(2,s) for s in Iterators.flatten(material_surfaces)],false,true,false))))
    material_points = unique(last.(gmsh.model.get_boundary(
        [(1,c) for c in material_curves],false,false,false)))
    registry = FEMLoopRegistry(model.fine_mesh_size)
    for tag in material_points
        point = Tuple(_coordinate_key.(gmsh.model.get_value(0,tag,Float64[])[1:2]))
        registry.points[point] = tag
        bucket = (floor(Int,point[1]/registry.mesh_size),floor(Int,point[2]/registry.mesh_size))
        push!(get!(Vector{Tuple{Float64,Float64}},registry.point_buckets,bucket),point)
    end
    for curve in material_curves
        gmsh.model.get_type(1,curve) == "Line" || continue
        lo,hi = gmsh.model.get_parametrization_bounds(1,curve)
        a,b = [registry.points[_matching_point_key(registry,
            gmsh.model.get_value(1,curve,[u]))] for u in (only(lo),only(hi))]
        registry.lines[minmax(a,b)] = curve
        registry.curve_points[curve] = (a,b)
    end
    centre_x, _ = model.centre
    halfwidth = mesh_plan.domain_halfwidth
    side, top, bottom = mesh_plan.pml_thickness
    # Split the Cartesian PML into columns at buried measurement paths. Each
    # path is then a shared transfinite edge, with the bottom PML grading.
    left, right = centre_x-halfwidth, centre_x+halfwidth
    buried_x = [endpoints[i][1] for i in eachindex(endpoints)
        if model.cable_hosts[model.problem.system.terminal_order[i].cable] !== :air]
    prescription = NamedTuple{(:side,:top,:bottom)}(mesh_plan.pml_strips)
    # Structural representability only; this does not validate field accuracy.
    for (strips,origin,distance) in ((prescription.side,right,side),
            (prescription.side,left,-side),(prescription.top,halfwidth,top),
            (prescription.bottom,-halfwidth,-bottom)), strip in strips
        _pml_progression(strip.count,strip.count*log(strip.ratio),
            origin+distance*strip.start,origin+distance*strip.stop)
    end
ns, nt, nb = length(prescription.side), length(prescription.top), length(prescription.bottom)
inner_xs = sort!(unique([left,buried_x...,right]))
xs = [left .- side .* reverse([s.stop for s in prescription.side]);
    inner_xs; right .+ side .* [s.stop for s in prescription.side]]

    inner_left = _point!(registry,(left,0.0);mesh_size=mesh_plan.domain_mesh_size)
    inner_right = _point!(registry,(right,0.0);mesh_size=mesh_plan.domain_mesh_size)
    air_outline = [(right,0.0), [(x,halfwidth) for x in reverse(inner_xs)]..., (left,0.0)]
    earth_outline = [(left,0.0), [(x,-halfwidth) for x in inner_xs]..., (right,0.0)]
    outlines = map((air_outline,earth_outline)) do points
        tags = [_point!(registry,p;mesh_size=mesh_plan.domain_mesh_size) for p in points]
        [_line!(registry,tags[k],tags[k+1]) for k in 1:length(tags)-1]
    end
    inner_curves = [outlines[1];outlines[2]]
    interface_sizes = _interface_mesh_sizes(model, mesh_plan)
    for x in abscissae
        interface_sizes[x] = min(get(interface_sizes,x,Inf),mesh_plan.interface_mesh_size)
    end
    interface_points = Int[inner_left]
    for (x, mesh_size) in sort!(collect(interface_sizes); by = first)
        centre_x - halfwidth < x < centre_x + halfwidth || continue
        push!(interface_points, _point!(registry, (x, 0.0); mesh_size))
    end
    push!(interface_points, inner_right)
    finite_interfaces = [_line!(registry, interface_points[index], interface_points[index + 1])
                         for index in 1:(length(interface_points) - 1)]
    air_holes = Int[]
    earth_holes = Int[]
    for cable_index in eachindex(cable_loops)
        target = model.cable_hosts[cable_index] === :air ? air_holes : earth_holes
        append!(target, cable_loops[cable_index])
    end
    air_loop = gmsh.model.geo.add_curve_loop([finite_interfaces; outlines[1]])
    earth_loop = gmsh.model.geo.add_curve_loop([-reverse(finite_interfaces); outlines[2]])
    air_surface = gmsh.model.geo.add_plane_surface([air_loop; air_holes])
    earth_surface = gmsh.model.geo.add_plane_surface([earth_loop; earth_holes])

    air_pml_surfaces, earth_pml_surfaces = Int[], Int[]
    interface_curves = copy(finite_interfaces)
    air_size, earth_size = mesh_plan.exterior_mesh_sizes
    vertical_grading = map((air_size, earth_size), mesh_plan.wave_mesh_sizes) do remote, wave
        first_size = remote == mesh_plan.domain_mesh_size ? remote : min(mesh_plan.domain_mesh_size, wave)
        _exterior_edge_grading(halfwidth, first_size, remote)
    end
ys = [-halfwidth .- bottom .* reverse([s.stop for s in prescription.bottom]);
    -halfwidth; 0.0; halfwidth; halfwidth .+ top .* [s.stop for s in prescription.top]]
    vertices = [_point!(registry, (x,y); mesh_size=mesh_plan.domain_mesh_size)
        for x in xs, y in ys]
    transfinite_curves = Dict{Int, Tuple{Int, Float64}}()
    transfinite_surfaces = Dict{Int, Tuple{String, NTuple{4, Int}}}()
    for partition in values(sector_partitions), (surface,corners) in partition.patches
        transfinite_surfaces[surface] = ("AlternateLeft",corners)
    end
    function mesh_line(first_point, last_point, count, ratio=1.0)
        curve = _line!(registry, first_point, last_point)
        coefficient = curve > 0 ? ratio : inv(ratio)
        assignment = (count+1, coefficient)
        if haskey(transfinite_curves, abs(curve))
            transfinite_curves[abs(curve)] == assignment || throw(ArgumentError(
                "inconsistent transfinite constraints on shared PML curve $(abs(curve))"))
            return curve
        end
        gmsh.model.geo.mesh.set_transfinite_curve(abs(curve), count+1,
            "Progression", coefficient)
        transfinite_curves[abs(curve)] = assignment
        return curve
    end
    outer_air_curves, outer_earth_curves = Int[], Int[]
    for i in 1:length(xs)-1, j in 1:length(ys)-1
    interior = ns < i < length(xs)-ns
    physical_row = nb < j <= nb+2
    interior && physical_row && continue
    medium = ys[j] >= 0 ? 1 : 2
    nx, rx = if interior
        (max(2,ceil(Int,(xs[i+1]-xs[i])/mesh_plan.exterior_mesh_sizes[medium])),1.0)
    else
        strip = i <= ns ? prescription.side[ns+1-i] : prescription.side[i-(length(xs)-ns-1)]
        (strip.count, i <= ns ? inv(strip.ratio) : strip.ratio)
    end
    ny, ry = if physical_row
        count, ratio = vertical_grading[medium]
        (count, medium == 2 ? inv(ratio) : ratio)
    else
        strip = j <= nb ? prescription.bottom[nb+1-j] : prescription.top[j-nb-2]
        (strip.count, j <= nb ? inv(strip.ratio) : strip.ratio)
    end
        a, b, c, d = vertices[i,j], vertices[i+1,j], vertices[i+1,j+1], vertices[i,j+1]
        lower = mesh_line(a,b,nx,rx)
        right = mesh_line(b,c,ny,ry)
        upper = mesh_line(d,c,nx,rx)
        left = mesh_line(a,d,ny,ry)
        loop = gmsh.model.geo.add_curve_loop([lower,right,-upper,-left])
        surface = gmsh.model.geo.add_plane_surface([loop])
        arrangement, corners = "AlternateLeft", (a,b,c,d)
        gmsh.model.geo.mesh.set_transfinite_surface(surface,arrangement,collect(corners))
        if mesh_plan.pml_element_family === :quadrangle
            gmsh.model.geo.mesh.set_recombine(2, surface)
        end
        transfinite_surfaces[surface] = (arrangement,corners)
        push!(medium == 1 ? air_pml_surfaces : earth_pml_surfaces, surface)
        outer = medium == 1 ? outer_air_curves : outer_earth_curves
        i == 1 && push!(outer,left)
        i == length(xs)-1 && push!(outer,right)
        j == 1 && push!(outer,lower)
        j == length(ys)-1 && push!(outer,upper)
        j == nb+1 && !interior && push!(interface_curves,upper)
    end
    pml_inner_curves = abs.(inner_curves)
    outer_air_curves = sort!(unique(abs.(outer_air_curves)))
    outer_earth_curves = sort!(unique(abs.(outer_earth_curves)))
    outer_curves = sort!([outer_air_curves;outer_earth_curves])
    interface_curves = sort!(unique(abs.(interface_curves)))

    _apply_point_mesh_sizes!(registry)
    gmsh.model.geo.synchronize()

    _validate_material_interfaces!(model, material_surfaces)

    for (index, material) in enumerate(model.material_plans)
        _physical_group(
            2,
            material_surfaces[index],
            material.physical_tag,
            material.physical_name
        )
    end
    terminal_curves = [_entity_boundary(surfaces) for surfaces in terminal_surfaces]
    # Each interval belongs to one CAD surface. Its endpoints share that
    # surface's boundary vertices; embedding makes its line elements actual
    # field edges. Metal intervals are excluded from the electric measurement.
    voltage_paths, voltage_references = Vector{Int}[], Int[]
    embeddings = Dict{Int,Vector{Int}}()
    surfaces = [(s,material.kind !== :conductor) for (i,material) in enumerate(model.material_plans)
        for s in material_surfaces[i]]
    append!(surfaces,[(air_surface,true),(earth_surface,true)])
    for (index, endpoint) in enumerate(endpoints)
        cable = model.problem.system.terminal_order[index].cable
        overhead = model.cable_hosts[cable] === :air
        reference_y = overhead ? 0.0 : -halfwidth-bottom
        reference = _point!(registry, (endpoint[1], reference_y);
            mesh_size=mesh_plan.interface_mesh_size)
        push!(voltage_references, reference)
        curves = Int[]
        physical_reference = if overhead
            reference
        else
            inner = _point!(registry, (endpoint[1], -halfwidth);
                mesh_size=mesh_plan.exterior_mesh_sizes[2])
        previous = reference
        for strip in reverse(prescription.bottom)
            next = _point!(registry, (endpoint[1], -halfwidth-bottom*strip.start);
                mesh_size=mesh_plan.exterior_mesh_sizes[2])
            push!(curves, mesh_line(previous, next, strip.count, inv(strip.ratio)))
            previous = next
        end

            inner
        end
        receiver = registry.points[_matching_point_key(registry,endpoint)]
        physical_y = overhead ? 0.0 : -halfwidth
        points = [(physical_y,physical_reference),(endpoint[2],receiver)]
        append!(points,[(point[2],tag) for (point,tag) in registry.points
            if abs(point[1]-endpoint[1]) <= 64eps(max(abs(endpoint[1]),1.0)) &&
               tag != physical_reference && tag != receiver &&
               physical_y < point[2] < endpoint[2]])
        sort!(points)
        remote_size = overhead ? mesh_plan.interface_mesh_size : mesh_plan.exterior_mesh_sizes[2]
        for k in 1:length(points)-1
            y0,a = points[k]; y1,b = points[k+1]
            y1 > y0 || continue
            middle = [endpoint[1],(y0+y1)/2,0.0]
            host = findfirst(pair -> gmsh.model.is_inside(2,first(pair),middle)>0,surfaces)
            host === nothing && _fem_error(:geometry,model.problem.system.system_id,
                :voltage_path,"measurement interval has no CAD host surface")
            surface, electric = surfaces[host]
            electric || continue
            key = minmax(a,b)
            existing = get(registry.lines,key,nothing)
            if existing === nothing
                # Only air/earth paths grow towards the remote boundary. A
                # cable interval uses local spacing: applying the far-field
                # ratio to a short insulation gap creates nearly coincident
                # nodes and can prevent Gmsh from recovering its field edges.
                last_size = surface in (air_surface,earth_surface) ?
                    remote_size : model.fine_mesh_size
                count, ratio = _exterior_edge_grading(y1-y0,model.fine_mesh_size,last_size)
                curve = mesh_line(a,b,count,inv(ratio))
                push!(get!(Vector{Int},embeddings,surface),abs(curve))
            else
                curve = existing
            end
            push!(curves,curve)
        end
        push!(voltage_paths, abs.(curves))
    end
    gmsh.model.geo.synchronize()
    for (surface,curves) in embeddings
        gmsh.model.mesh.embed(1,curves,2,surface)
    end
    for index in eachindex(terminal_surfaces)
        _physical_group(
            2,
            terminal_surfaces[index],
            model.tags.terminal_base + index,
            model.terminal_names[index]
        )
        _physical_group(
            1,
            terminal_curves[index],
            model.tags.terminal_contour_base + index,
            @sprintf("LCM/terminal_contour/%04d", index)
        )
        _physical_group(1, voltage_paths[index], model.tags.voltage_path_base + index,
            @sprintf("LCM/voltage_path/%04d", index))
        _physical_group(0, [voltage_references[index]], model.tags.voltage_reference_base + index,
            @sprintf("LCM/voltage_reference/%04d", index))
    end
    for index in eachindex(cable_curves)
        _physical_group(
            1,
            cable_curves[index],
            model.tags.cable_contour_base + index,
            @sprintf("LCM/cable_contour/%04d/%s", index,
                model.problem.system.designs[index].cable_id)
        )
    end
    _physical_group(2, [air_surface], model.tags.air, "LCM/domain/air")
    _physical_group(2, [earth_surface], model.tags.earth, "LCM/domain/earth")
    _physical_group(2, air_pml_surfaces, model.tags.air_pml, "LCM/domain/air_pml")
    _physical_group(2, earth_pml_surfaces, model.tags.earth_pml, "LCM/domain/earth_pml")
    _physical_group(2, [air_pml_surfaces; earth_pml_surfaces], model.tags.pml, "LCM/domain/pml")
    _physical_group(1, outer_curves, model.tags.outer_boundary, "LCM/boundary/magnetic_dirichlet")
    _physical_group(1, outer_air_curves, model.tags.outer_air_boundary, "LCM/boundary/electric_reference_air")
    _physical_group(1, outer_earth_curves, model.tags.outer_earth_boundary, "LCM/boundary/electric_reference_earth")
    _physical_group(1, pml_inner_curves, model.tags.pml_inner_boundary, "LCM/boundary/pml_inner")
    _physical_group(1, interface_curves, model.tags.interface, "LCM/interface/air_earth")
    conductor_surfaces = reduce(vcat,
        [material_surfaces[index]
         for (index, material) in enumerate(model.material_plans)
         if material.kind === :conductor];
        init = Int[])
    passive_surfaces = reduce(vcat,
        [material_surfaces[index]
         for (index, material) in enumerate(model.material_plans)
         if material.kind !== :conductor];
        init = Int[])
    _physical_group(2, conductor_surfaces, 6_001, "LCM/domain/conductors")
    _physical_group(2, passive_surfaces, 6_002, "LCM/domain/passive_cable_media")
    _physical_group(
        2,
        [conductor_surfaces;
         passive_surfaces;
         air_surface;
         earth_surface;
         air_pml_surfaces;
         earth_pml_surfaces],
        6_003,
        "LCM/domain/field_maps"
    )

    return FEMGeometry(
        model_name,
        terminal_surfaces,
        terminal_curves,
        material_surfaces,
        region_surfaces,
        cable_curves,
        [air_surface; air_pml_surfaces],
        [earth_surface; earth_pml_surfaces],
        [air_pml_surfaces; earth_pml_surfaces],
        outer_curves,
        outer_air_curves,
        outer_earth_curves,
        pml_inner_curves,
        interface_curves,
        cable_loops,
        sort!(setdiff(collect(values(registry.lines)),material_curves)),
        sort!(setdiff(collect(values(registry.points)),material_points)),
        transfinite_curves,
        transfinite_surfaces,
        Dict{Int, Tuple{Int, Int}}(),
        sector_partitions
    )
end
