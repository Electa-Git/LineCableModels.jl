struct FEMLoop
    ccw::Int
    cw::Int
    curves::Vector{Int}
end

const FEM_BOUNDARY_ANGLE_TOLERANCE = 1.0e-12

mutable struct FEMLoopRegistry
    points::Dict{Tuple{Float64, Float64}, Int}
    point_sizes::Dict{Int, Float64}
    lines::Dict{Tuple{Int, Int}, Int}
    curve_points::Dict{Int, Set{Int}}
    circle_arcs::Dict{Any, Tuple{Int, Int, Int}}
    circle_breaks::Dict{Any, Set{Float64}}
    circle_break_points::Dict{Any, Dict{Float64, Tuple{Float64, Float64}}}
    loops::Dict{Any, FEMLoop}
    mesh_size::Float64
end

function FEMLoopRegistry(mesh_size)
    FEMLoopRegistry(
        Dict{Tuple{Float64, Float64}, Int}(),
        Dict{Int, Float64}(),
        Dict{Tuple{Int, Int}, Int}(),
        Dict{Int, Set{Int}}(),
        Dict{Any, Tuple{Int, Int, Int}}(),
        Dict{Any, Set{Float64}}(),
        Dict{Any, Dict{Float64, Tuple{Float64, Float64}}}(),
        Dict{Any, FEMLoop}(),
        Float64(mesh_size)
    )
end

_coordinate_key(value) = round(Float64(value); sigdigits = 15)

function _matching_point_key(points, point)
    key = (_coordinate_key(point[1]), _coordinate_key(point[2]))
    haskey(points, key) && return key
    scale = max(abs(Float64(point[1])), abs(Float64(point[2])), 1.0)
    tolerance = 64eps(scale)
    for candidate in keys(points)
        abs(candidate[1] - key[1]) <= tolerance &&
            abs(candidate[2] - key[2]) <= tolerance && return candidate
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
    shape isa DataModel._DifferencePrimitive || return nothing
    outer = shape.outer
    outer isa DataModel.Annulus || return nothing
    isempty(shape.holes) && return nothing
    iszero(outer.ri) && return nothing

    centre = (Float64(outer.at.x), Float64(outer.at.y))
    disks = filter(hole -> hole isa DataModel.Disk, shape.holes)
    cutouts = filter(hole -> !(hole isa DataModel.Disk), shape.holes)
    isempty(disks) && return nothing
    all(hole -> hole isa DataModel.Sector, cutouts) || return nothing
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
    all(radius -> isapprox(
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
        separations = [
            mod(angles[mod1(index + 1, end)] - angles[index], 2π)
            for index in eachindex(angles)
        ]
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
        start = _angle_key(Float64(cutout.at.φ + cutout.φ0))
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
    push!(curves, _line!(registry, outer_stop, inner_stop))
    append!(curves, _circle_arc_path!(
        registry,
        centre,
        inner_radius,
        start_angle + span,
        -span;
        mesh_size,
        first_point = inner_stop_point,
        last_point = inner_start_point
    ))
    push!(curves, _line!(registry, inner_start, outer_start))
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

function _register_shape_breaks!(registry::FEMLoopRegistry, shape)
    if shape isa DataModel.Disk
        _register_full_circle_breaks!(
            registry, (shape.at.x, shape.at.y), shape.r
        )
        return nothing
    end
    if shape isa DataModel.Annulus
        centre = (shape.at.x, shape.at.y)
        _register_full_circle_breaks!(registry, centre, shape.ri)
        _register_full_circle_breaks!(registry, centre, shape.ro)
        return nothing
    end
    if shape isa DataModel.Sector
        centre = (shape.at.x, shape.at.y)
        start_angle = shape.at.φ + shape.φ0
        stop_angle = start_angle + shape.span
        _register_circle_break!(registry, centre, shape.ro, start_angle)
        _register_circle_break!(registry, centre, shape.ro, stop_angle)
        if !iszero(shape.ri)
            _register_circle_break!(registry, centre, shape.ri, start_angle)
            _register_circle_break!(registry, centre, shape.ri, stop_angle)
        end
        return nothing
    end
    if shape isa DataModel.RoundedSectorShape
        for arc in (
                shape.contacts.arcs.base,
                shape.contacts.arcs.lower,
                shape.contacts.arcs.back,
                shape.contacts.arcs.upper
        )
            iszero(arc.radius) && continue
            centre = _transform_point(arc.center, shape.at)
            _register_circle_break!(
                registry, centre, arc.radius, shape.at.φ + arc.start
            )
            _register_circle_break!(
                registry, centre, arc.radius, shape.at.φ + arc.stop
            )
        end
        return nothing
    end
    if hasproperty(shape, :inner) && hasproperty(shape, :outer)
        _register_shape_breaks!(registry, getproperty(shape, :inner))
        _register_shape_breaks!(registry, getproperty(shape, :outer))
        return nothing
    end
    if hasproperty(shape, :holes) && hasproperty(shape, :outer)
        _register_shape_breaks!(registry, getproperty(shape, :outer))
        foreach(
            hole -> _register_shape_breaks!(registry, hole),
            getproperty(shape, :holes)
        )
        _register_tangent_fill_contacts!(registry, shape)
        return nothing
    end
    if hasproperty(shape, :members)
        foreach(
            member -> _register_shape_breaks!(registry, member),
            getproperty(shape, :members)
        )
    end
    return nothing
end

function _register_circle_contacts!(registry::FEMLoopRegistry)
    circles = collect(keys(registry.circle_breaks))
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
    key = _matching_point_key(registry.points, point)
    tag = get!(registry.points, key) do
        gmsh.model.geo.add_point(key[1], key[2], 0.0, Float64(mesh_size))
    end
    registry.point_sizes[tag] = min(
        get(registry.point_sizes, tag, Inf), Float64(mesh_size)
    )
    return tag
end

function _line!(registry::FEMLoopRegistry, first_point::Int, last_point::Int)
    key = minmax(first_point, last_point)
    tag = get!(registry.lines, key) do
        gmsh.model.geo.add_line(key[1], key[2])
    end
    points = get!(registry.curve_points, abs(tag)) do
        Set{Int}()
    end
    push!(points, first_point, last_point)
    return first_point == key[1] ? tag : -tag
end

function _refine_loop_mesh_size!(
        registry::FEMLoopRegistry, loop::FEMLoop, mesh_size
)
    value = Float64(mesh_size)
    for curve in loop.curves
        for point in get(registry.curve_points, abs(curve), Set{Int}())
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
    return sum(eachindex(points)) do index
        next = mod1(index + 1, length(points))
        points[index][1] * points[next][2] -
        points[next][1] * points[index][2]
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
        point_tags = [_point!(registry, point; mesh_size) for point in values]
        curves = [_line!(registry, point_tags[index], point_tags[mod1(index + 1, end)])
                  for index in eachindex(point_tags)]
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        FEMLoop(ccw, cw, abs.(curves))
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
        FEMLoop(ccw, cw, abs.(curves))
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
        point_tags = [_point!(registry, _transform_point(point, shape.at); mesh_size)
                      for point in local_points]
        major_point = shape.a >= shape.b ? local_points[1] : local_points[2]
        major_tag = _point!(
            registry, _transform_point(major_point, shape.at); mesh_size
        )
        curves = Int[]
        for index in 1:4
            next_index = mod1(index + 1, 4)
            curve = gmsh.model.geo.add_ellipse_arc(
                point_tags[index], centre_tag, major_tag, point_tags[next_index]
            )
            registry.curve_points[curve] = Set((
                point_tags[index], point_tags[next_index]
            ))
            push!(curves, curve)
        end
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        FEMLoop(ccw, cw, curves)
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
                tag, stored_first, stored_last = get!(registry.circle_arcs, key) do
                    (
                        gmsh.model.geo.add_circle_arc(
                            first_tag, centre_tag, last_tag
                        ),
                        first_tag,
                        last_tag
                    )
                end
                points_for_curve = get!(registry.curve_points, abs(tag)) do
                    Set{Int}()
                end
                push!(points_for_curve, first_tag, last_tag)
                first_tag == stored_first && last_tag == stored_last ? tag : -tag
            end
            for index in 1:(length(point_tags) - 1)]
end

function _sector_loop!(
        registry::FEMLoopRegistry,
        shape::DataModel.Sector;
        mesh_size = registry.mesh_size
)
    full = isapprox(shape.span, 2π; rtol = 0, atol = 128eps(Float64))
    centre = (shape.at.x, shape.at.y)
    full && return _circle_loop!(registry, centre, shape.ro; mesh_size)

    start_angle = shape.at.φ + shape.φ0
    key = (
        :sector,
        _coordinate_key(centre[1]),
        _coordinate_key(centre[2]),
        _coordinate_key(start_angle),
        _coordinate_key(shape.span),
        _coordinate_key(shape.ri),
        _coordinate_key(shape.ro)
    )
    loop = get!(registry.loops, key) do
        outer = _circle_arc_path!(
            registry,
            centre,
            shape.ro,
            start_angle,
            shape.span;
            mesh_size
        )
        outer_start = _point!(
            registry,
            _circle_point(centre, shape.ro, start_angle);
            mesh_size
        )
        outer_stop = _point!(
            registry,
            _circle_point(centre, shape.ro, start_angle + shape.span);
            mesh_size
        )

        curves = copy(outer)
        if iszero(shape.ri)
            centre_tag = _point!(registry, centre; mesh_size)
            push!(curves, _line!(registry, outer_stop, centre_tag))
            push!(curves, _line!(registry, centre_tag, outer_start))
        else
            inner_stop_point = _circle_point(
                centre, shape.ri, start_angle + shape.span
            )
            inner_start_point = _circle_point(centre, shape.ri, start_angle)
            inner_stop = _point!(registry, inner_stop_point; mesh_size)
            inner_start = _point!(registry, inner_start_point; mesh_size)
            push!(curves, _line!(registry, outer_stop, inner_stop))
            append!(curves,
                _circle_arc_path!(
                    registry,
                    centre,
                    shape.ri,
                    start_angle + shape.span,
                    -shape.span;
                    mesh_size,
                    first_point = inner_stop_point,
                    last_point = inner_start_point
                ))
            push!(curves, _line!(registry, inner_start, outer_start))
        end
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        FEMLoop(ccw, cw, abs.(curves))
    end
    return _refine_loop_mesh_size!(registry, loop, mesh_size)
end

function _rounded_sector_loop!(
        registry::FEMLoopRegistry,
        shape::DataModel.RoundedSectorShape;
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
    key = (:rounded_sector,
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
                                   Int[_line!(
                    registry, point_tags[first_index], point_tags[last_index]
                )]
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
        push!(curves, _line!(registry, point_tags[2], point_tags[3]))
        append!(curves, arc_curves[:lower])
        append!(curves, arc_curves[:back])
        append!(curves, arc_curves[:upper])
        push!(curves, _line!(registry, point_tags[6], point_tags[1]))
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        FEMLoop(ccw, cw, abs.(curves))
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

function _shape_points(shape)
    if hasproperty(shape, :contacts) && hasproperty(shape, :primitive)
        return DataModel.tessellate(shape; points_per_arc = 24)
    end
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
    shape isa DataModel.Sector && return _sector_loop!(
        registry, shape; mesh_size
    )
    shape isa DataModel.RoundedSectorShape && return _rounded_sector_loop!(
        registry, shape; mesh_size
    )
    return _polygon_loop!(registry, _shape_points(shape); mesh_size)
end

function _surface_boundaries(shape)
    if shape isa DataModel.Annulus
        outer = DataModel.Disk(shape.ro, shape.at)
        inner = DataModel.Disk(shape.ri, shape.at)
        return outer, Any[inner]
    end
    if shape isa DataModel.Sector &&
       isapprox(
           shape.span, 2π; rtol = 0, atol = 128eps(Float64)
       ) && !iszero(shape.ri)
        return DataModel.Disk(shape.ro, shape.at), Any[DataModel.Disk(shape.ri, shape.at)]
    end
    if hasproperty(shape, :inner) && hasproperty(shape, :outer)
        return getproperty(shape, :outer), Any[getproperty(shape, :inner)]
    end
    if hasproperty(shape, :holes) && hasproperty(shape, :outer)
        outer, outer_holes = _surface_boundaries(getproperty(shape, :outer))
        return outer, Any[outer_holes..., getproperty(shape, :holes)...]
    end
    return shape, Any[]
end

function _surface!(registry::FEMLoopRegistry, shape, mesh_size)
    outer, holes = _surface_boundaries(shape)
    outer_loop = _boundary_loop!(registry, outer; mesh_size)
    hole_loops = [_boundary_loop!(registry, hole; mesh_size).cw
                  for hole in holes]
    return gmsh.model.geo.add_plane_surface([outer_loop.ccw; hole_loops])
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
        append!(curves, _circle_arc_path!(
            registry,
            plan.centre,
            plan.outer_radius,
            angle,
            span;
            mesh_size,
            first_point = outer_point,
            last_point = next_outer_point
        ))
        append!(curves, _circle_arc_path!(
            registry,
            next_hole_centre,
            Float64(next_hole.r),
            next_angle,
            -π;
            mesh_size,
            first_point = next_outer_point,
            last_point = next_inner_point
        ))
        append!(curves, _circle_arc_path!(
            registry,
            plan.centre,
            plan.inner_radius,
            next_angle,
            -span;
            mesh_size,
            first_point = next_inner_point,
            last_point = inner_point
        ))
        append!(curves, _circle_arc_path!(
            registry,
            hole_centre,
            Float64(hole.r),
            angle + π,
            -π;
            mesh_size,
            first_point = inner_point,
            last_point = outer_point
        ))
        loop = gmsh.model.geo.add_curve_loop(curves)
        push!(surfaces, gmsh.model.geo.add_plane_surface([loop]))
    end
    for (start_angle, span) in plan.outer_spans
        push!(surfaces, _annular_sector_surface!(
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
    return Int[_surface!(registry, shape, mesh_size)]
end

function _boundary_components(shape)
    if hasproperty(shape, :members)
        return collect(getproperty(shape, :members))
    end
    return Any[shape]
end

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
    curves = Int[]
    for surfaces in material_surfaces
        append!(curves, _entity_boundary(surfaces))
    end
    offending = Int[]
    for curve in sort!(unique(curves))
        adjacent, _ = gmsh.model.get_adjacencies(1, curve)
        length(unique(adjacent)) >= 2 || push!(offending, curve)
    end
    isempty(offending) || _fem_error(
        :geometry,
        model.problem.system.system_id,
        :material_partition,
        "material interfaces are not conformal; curves with no adjacent " *
        "field surface: $(join(offending, ", "))"
    )
    return nothing
end

function _build_geometry!(
        model::FEMResolvedModel,
        model_name::String,
        mesh_plan::FEMMeshPlan = last(model.mesh_plans)
)
    gmsh.model.add(model_name)
    registry = FEMLoopRegistry(model.fine_mesh_size)

    for region in model.region_plans
        _register_shape_breaks!(registry, region.shape)
    end
    foreach(
        boundary -> _register_shape_breaks!(registry, boundary),
        model.cable_boundaries
    )
    _register_circle_contacts!(registry)

    material_surfaces = [Int[] for _ in model.material_plans]
    terminal_surfaces = [Int[] for _ in model.terminal_ids]
    for region in model.region_plans
        surfaces = _surfaces!(registry, region.shape, region.mesh_size)
        append!(material_surfaces[region.material_index], surfaces)
        region.terminal_index > 0 && append!(
            terminal_surfaces[region.terminal_index], surfaces
        )
    end

    cable_curves = [Int[] for _ in model.cable_boundaries]
    cable_loops = [FEMLoop[] for _ in model.cable_boundaries]
    for cable_index in eachindex(model.cable_boundaries)
        for component in _boundary_components(model.cable_boundaries[cable_index])
            loop = _boundary_loop!(
                registry,
                component;
                mesh_size = model.cable_outer_mesh_sizes[cable_index]
            )
            push!(cable_loops[cable_index], loop)
            append!(cable_curves[cable_index], loop.curves)
        end
    end

    centre_x, _ = model.centre
    radius = mesh_plan.domain_radius
    shell_radius = mesh_plan.shell_outer_radius
    inner = _circle_loop!(
        registry,
        (centre_x, 0.0),
        radius;
        mesh_size = mesh_plan.domain_mesh_size
    )
    outer = _circle_loop!(
        registry,
        (centre_x, 0.0),
        shell_radius;
        mesh_size = mesh_plan.infinite_mesh_size
    )
    inner_left = _point!(
        registry,
        (centre_x - radius, 0.0);
        mesh_size = mesh_plan.domain_mesh_size
    )
    inner_right = _point!(
        registry,
        (centre_x + radius, 0.0);
        mesh_size = mesh_plan.domain_mesh_size
    )
    outer_left = _point!(
        registry,
        (centre_x - shell_radius, 0.0);
        mesh_size = mesh_plan.infinite_mesh_size
    )
    outer_right = _point!(
        registry,
        (centre_x + shell_radius, 0.0);
        mesh_size = mesh_plan.infinite_mesh_size
    )
    interface_sizes = Dict{Float64, Float64}()
    function register_interface_size(x, mesh_size)
        key = _coordinate_key(x)
        interface_sizes[key] = min(
            get(interface_sizes, key, Inf), Float64(mesh_size)
        )
        return nothing
    end
    register_interface_size(centre_x, mesh_plan.interface_mesh_size)
    for offset in (-2.0, 2.0)
        x = centre_x + offset
        centre_x - radius < x < centre_x + radius &&
            register_interface_size(x, mesh_plan.domain_mesh_size)
    end
    for (cable_index, position) in enumerate(model.problem.system.positions)
        centre_x - radius < position.x < centre_x + radius &&
            register_interface_size(
                position.x,
                mesh_plan.cable_interface_mesh_sizes[cable_index]
            )
    end
    interface_points = Int[inner_left]
    for (x, mesh_size) in sort!(collect(interface_sizes); by = first)
        centre_x - radius < x < centre_x + radius || continue
        push!(interface_points, _point!(registry, (x, 0.0); mesh_size))
    end
    push!(interface_points, inner_right)
    finite_interfaces = [
        _line!(registry, interface_points[index], interface_points[index + 1])
        for index in 1:(length(interface_points) - 1)
    ]
    left_connector = _line!(registry, outer_left, inner_left)
    right_connector = _line!(registry, inner_right, outer_right)
    air_holes = Int[]
    earth_holes = Int[]
    for cable_index in eachindex(cable_loops)
        target = model.cable_hosts[cable_index] === :air ? air_holes : earth_holes
        append!(target, getproperty.(cable_loops[cable_index], :cw))
    end
    air_loop = gmsh.model.geo.add_curve_loop([
        finite_interfaces; inner.curves[1]; inner.curves[2]
    ])
    earth_loop = gmsh.model.geo.add_curve_loop([
        -reverse(finite_interfaces); inner.curves[3]; inner.curves[4]
    ])
    air_surface = gmsh.model.geo.add_plane_surface([air_loop; air_holes])
    earth_surface = gmsh.model.geo.add_plane_surface([earth_loop; earth_holes])

    air_infinite_loop = gmsh.model.geo.add_curve_loop([
        right_connector,
        outer.curves[1],
        outer.curves[2],
        left_connector,
        -inner.curves[2],
        -inner.curves[1]
    ])
    earth_infinite_loop = gmsh.model.geo.add_curve_loop([
        outer.curves[3],
        outer.curves[4],
        -right_connector,
        -inner.curves[4],
        -inner.curves[3],
        -left_connector
    ])
    air_infinite_surface = gmsh.model.geo.add_plane_surface([air_infinite_loop])
    earth_infinite_surface = gmsh.model.geo.add_plane_surface([
        earth_infinite_loop
    ])

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
    _physical_group(
        2,
        [air_infinite_surface],
        model.tags.air_infinite,
        "LCM/domain/air_infinite"
    )
    _physical_group(
        2,
        [earth_infinite_surface],
        model.tags.earth_infinite,
        "LCM/domain/earth_infinite"
    )
    _physical_group(
        2,
        [air_infinite_surface, earth_infinite_surface],
        model.tags.infinite_domain,
        "LCM/domain/infinite_shell"
    )

    interface_curves = sort!(unique(abs.([
        finite_interfaces; left_connector; right_connector
    ])))
    outer_curves = sort!(unique(outer.curves))
    outer_air_curves = intersect(
        outer_curves,
        _entity_boundary([air_infinite_surface])
    )
    outer_earth_curves = intersect(
        outer_curves,
        _entity_boundary([earth_infinite_surface])
    )
    isempty(intersect(outer_air_curves, outer_earth_curves)) || _fem_error(
        :geometry,
        model.problem.system.system_id,
        :outer_electric_boundary,
        "air and earth outer electric boundaries overlap"
    )
    sort!(unique([outer_air_curves; outer_earth_curves])) == outer_curves ||
        _fem_error(
            :geometry,
            model.problem.system.system_id,
            :outer_electric_boundary,
            "air and earth outer electric boundaries do not partition the outer shell"
        )
    inner_shell_curves = sort!(unique(inner.curves))
    _physical_group(
        1,
        outer_curves,
        model.tags.outer_boundary,
        "LCM/boundary/magnetic_dirichlet"
    )
    _physical_group(
        1,
        outer_air_curves,
        model.tags.outer_air_boundary,
        "LCM/boundary/electric_insulation_air"
    )
    _physical_group(
        1,
        outer_earth_curves,
        model.tags.outer_earth_boundary,
        "LCM/boundary/electric_reference_earth"
    )
    _physical_group(
        1,
        inner_shell_curves,
        model.tags.inner_shell_boundary,
        "LCM/boundary/inner_infinite_shell"
    )
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
         air_infinite_surface;
         earth_infinite_surface],
        6_003,
        "LCM/domain/field_maps"
    )

    return FEMGeometry(
        model_name,
        terminal_surfaces,
        terminal_curves,
        material_surfaces,
        cable_curves,
        [air_surface, air_infinite_surface],
        [earth_surface, earth_infinite_surface],
        [air_infinite_surface, earth_infinite_surface],
        outer_curves,
        outer_air_curves,
        outer_earth_curves,
        inner_shell_curves,
        interface_curves
    )
end
