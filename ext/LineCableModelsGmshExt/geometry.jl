struct FEMLoop
    ccw::Int
    cw::Int
    curves::Vector{Int}
end

mutable struct FEMLoopRegistry
    points::Dict{Tuple{Float64, Float64}, Int}
    lines::Dict{Tuple{Int, Int}, Int}
    loops::Dict{Any, FEMLoop}
    mesh_size::Float64
end

function FEMLoopRegistry(mesh_size)
    FEMLoopRegistry(
        Dict{Tuple{Float64, Float64}, Int}(),
        Dict{Tuple{Int, Int}, Int}(),
        Dict{Any, FEMLoop}(),
        Float64(mesh_size)
    )
end

_coordinate_key(value) = round(Float64(value); sigdigits = 15)

function _transform_point(point, at)
    cosine = cos(at.φ)
    sine = sin(at.φ)
    return (
        at.x + cosine * point[1] - sine * point[2],
        at.y + sine * point[1] + cosine * point[2]
    )
end

function _point!(registry::FEMLoopRegistry, point; mesh_size = registry.mesh_size)
    key = (_coordinate_key(point[1]), _coordinate_key(point[2]))
    return get!(registry.points, key) do
        gmsh.model.geo.add_point(key[1], key[2], 0.0, Float64(mesh_size))
    end
end

function _line!(registry::FEMLoopRegistry, first_point::Int, last_point::Int)
    key = minmax(first_point, last_point)
    tag = get!(registry.lines, key) do
        gmsh.model.geo.add_line(key[1], key[2])
    end
    return first_point == key[1] ? tag : -tag
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
    return get!(registry.loops, key) do
        point_tags = [_point!(registry, point; mesh_size) for point in values]
        curves = [_line!(registry, point_tags[index], point_tags[mod1(index + 1, end)])
                  for index in eachindex(point_tags)]
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        FEMLoop(ccw, cw, abs.(curves))
    end
end

function _circle_loop!(
        registry::FEMLoopRegistry,
        centre,
        radius;
        mesh_size = registry.mesh_size
)
    key = (
        :circle,
        _coordinate_key(centre[1]),
        _coordinate_key(centre[2]),
        _coordinate_key(radius)
    )
    return get!(registry.loops, key) do
        centre_tag = _point!(registry, centre; mesh_size)
        points = [
            (centre[1] + radius, centre[2]),
            (centre[1], centre[2] + radius),
            (centre[1] - radius, centre[2]),
            (centre[1], centre[2] - radius)
        ]
        point_tags = [_point!(registry, point; mesh_size) for point in points]
        curves = [gmsh.model.geo.add_circle_arc(
                      point_tags[index], centre_tag, point_tags[mod1(index + 1, 4)]
                  )
                  for index in 1:4]
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        FEMLoop(ccw, cw, curves)
    end
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
    return get!(registry.loops, key) do
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
        curves = [gmsh.model.geo.add_ellipse_arc(
                      point_tags[index],
                      centre_tag,
                      major_tag,
                      point_tags[mod1(index + 1, 4)]
                  )
                  for index in 1:4]
        ccw = gmsh.model.geo.add_curve_loop(curves)
        cw = gmsh.model.geo.add_curve_loop(-reverse(curves))
        FEMLoop(ccw, cw, curves)
    end
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
    count = max(1, ceil(Int, abs(Float64(span)) / (π / 2)))
    angles = range(start_angle, start_angle + span; length = count + 1)
    points = [(
                  centre[1] + radius * cos(angle),
                  centre[2] + radius * sin(angle)
              )
              for angle in angles]
    first_point === nothing || (points[1] = first_point)
    last_point === nothing || (points[end] = last_point)
    centre_tag = _point!(registry, centre; mesh_size)
    point_tags = [_point!(registry, point; mesh_size) for point in points]
    return [gmsh.model.geo.add_circle_arc(
                point_tags[index], centre_tag, point_tags[index + 1]
            )
            for index in 1:count]
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
    return get!(registry.loops, key) do
        outer = _circle_arc_path!(
            registry,
            centre,
            shape.ro,
            start_angle,
            shape.span;
            mesh_size
        )
        outer_start = _point!(registry,
            (
                centre[1] + shape.ro * cos(start_angle),
                centre[2] + shape.ro * sin(start_angle)
            );
            mesh_size)
        outer_stop = _point!(registry,
            (
                centre[1] + shape.ro * cos(start_angle + shape.span),
                centre[2] + shape.ro * sin(start_angle + shape.span)
            );
            mesh_size)

        curves = copy(outer)
        if iszero(shape.ri)
            centre_tag = _point!(registry, centre; mesh_size)
            push!(curves, _line!(registry, outer_stop, centre_tag))
            push!(curves, _line!(registry, centre_tag, outer_start))
        else
            inner_stop_point = (
                centre[1] + shape.ri * cos(start_angle + shape.span),
                centre[2] + shape.ri * sin(start_angle + shape.span)
            )
            inner_start_point = (
                centre[1] + shape.ri * cos(start_angle),
                centre[2] + shape.ri * sin(start_angle)
            )
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
    return get!(registry.loops, key) do
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
        return getproperty(shape, :outer), Any[getproperty(shape, :holes)...]
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

function _build_geometry!(model::FEMResolvedModel, model_name::String)
    gmsh.model.add(model_name)
    registry = FEMLoopRegistry(model.fine_mesh_size)

    material_surfaces = [Int[] for _ in model.material_plans]
    terminal_surfaces = [Int[] for _ in model.terminal_ids]
    for region in model.region_plans
        surface = _surface!(registry, region.shape, model.fine_mesh_size)
        push!(material_surfaces[region.material_index], surface)
        region.terminal_index > 0 && push!(
            terminal_surfaces[region.terminal_index], surface
        )
    end

    cable_curves = [Int[] for _ in model.cable_boundaries]
    cable_loops = [FEMLoop[] for _ in model.cable_boundaries]
    for cable_index in eachindex(model.cable_boundaries)
        for component in _boundary_components(model.cable_boundaries[cable_index])
            loop = _boundary_loop!(registry, component; mesh_size = model.fine_mesh_size)
            push!(cable_loops[cable_index], loop)
            append!(cable_curves[cable_index], loop.curves)
        end
    end

    centre_x, _ = model.centre
    radius = model.domain_radius
    shell_radius = model.shell_outer_radius
    inner = _circle_loop!(
        registry, (centre_x, 0.0), radius; mesh_size = model.coarse_mesh_size
    )
    outer = _circle_loop!(
        registry,
        (centre_x, 0.0),
        shell_radius;
        mesh_size = model.coarse_mesh_size
    )
    inner_left = _point!(
        registry, (centre_x - radius, 0.0); mesh_size = model.coarse_mesh_size
    )
    inner_right = _point!(
        registry, (centre_x + radius, 0.0); mesh_size = model.coarse_mesh_size
    )
    outer_left = _point!(
        registry,
        (centre_x - shell_radius, 0.0);
        mesh_size = model.coarse_mesh_size
    )
    outer_right = _point!(
        registry,
        (centre_x + shell_radius, 0.0);
        mesh_size = model.coarse_mesh_size
    )
    finite_interface = _line!(registry, inner_left, inner_right)
    left_connector = _line!(registry, outer_left, inner_left)
    right_connector = _line!(registry, inner_right, outer_right)
    air_holes = Int[]
    earth_holes = Int[]
    for cable_index in eachindex(cable_loops)
        target = model.cable_hosts[cable_index] === :air ? air_holes : earth_holes
        append!(target, getproperty.(cable_loops[cable_index], :cw))
    end
    air_loop = gmsh.model.geo.add_curve_loop([
        finite_interface, inner.curves[1], inner.curves[2]
    ])
    earth_loop = gmsh.model.geo.add_curve_loop([
        -finite_interface, inner.curves[3], inner.curves[4]
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

    gmsh.model.geo.synchronize()

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

    interface_curves = sort!(unique([
        abs(finite_interface), abs(left_connector), abs(right_connector)
    ]))
    outer_curves = sort!(unique(outer.curves))
    inner_shell_curves = sort!(unique(inner.curves))
    _physical_group(
        1, outer_curves, model.tags.outer_boundary, "LCM/boundary/dirichlet"
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
        inner_shell_curves,
        interface_curves
    )
end
