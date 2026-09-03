const _RHO_MIN = 1.0e-9
const _RHO_MAX = 1.0e10

_finite_point(point) = isfinite(point[1]) && isfinite(point[2])

function _circle_points(radius, xcenter, ycenter; count::Int = 128)
    angles = range(0, 2π; length = count)
    return filter(
        _finite_point,
        Point2f.(xcenter .+ radius .* cos.(angles), ycenter .+ radius .* sin.(angles))
    )
end

function _annulus_polygon(inner_radius, outer_radius, xcenter, ycenter; count::Int = 256)
    outer = _circle_points(outer_radius, xcenter, ycenter; count)
    iszero(inner_radius) && return GeometryBasics.Polygon(outer)
    inner = reverse(_circle_points(inner_radius, xcenter, ycenter; count))
    return GeometryBasics.Polygon(outer, [inner])
end

function _transform_preview_points(points, pose::Pose2)
    cosine = cos(nominal(pose.φ))
    sine = sin(nominal(pose.φ))
    x = nominal(pose.x)
    y = nominal(pose.y)
    return Point2f[(
                       x + cosine * nominal(point[1]) - sine * nominal(point[2]),
                       y + sine * nominal(point[1]) + cosine * nominal(point[2])
                   ) for point in points]
end

function _primitive_geometry(primitive::Disk)
    return GeometryBasics.Polygon(_circle_points(nominal(primitive.r), 0.0, 0.0))
end

function _primitive_geometry(primitive::Annulus)
    return _annulus_polygon(
        nominal(primitive.ri), nominal(primitive.ro), 0.0, 0.0
    )
end

function _primitive_geometry(primitive::Sector)
    outer_angles = range(primitive.φ0, primitive.φ0 + primitive.span; length = 96)
    inner_angles = reverse(outer_angles)
    points = Point2f[(primitive.ro * cos(angle), primitive.ro * sin(angle))
                     for angle in outer_angles]
    append!(points,
        Point2f[(primitive.ri * cos(angle), primitive.ri * sin(angle))
                for angle in inner_angles])
    return GeometryBasics.Polygon(points)
end

function _primitive_geometry(primitive::Rectangle)
    x = nominal(primitive.w) / 2
    y = nominal(primitive.h) / 2
    return GeometryBasics.Polygon(Point2f[(-x, -y), (x, -y), (x, y), (-x, y)])
end

function _primitive_geometry(primitive::Ellipse)
    angles = range(0, 2pi; length = 128)
    return GeometryBasics.Polygon(Point2f[(nominal(primitive.a) * cos(angle),
                                              nominal(primitive.b) * sin(angle))
                                          for angle in angles])
end

function _primitive_geometry(primitive::Polygon)
    return GeometryBasics.Polygon(Point2f[point for point in primitive.points])
end

function _shape_geometry(primitive::AbstractPrimitive)
    geometry = _primitive_geometry(primitive)
    exterior = _transform_preview_points(geometry.exterior, primitive.at)
    interiors = [_transform_preview_points(interior, primitive.at)
                 for interior in geometry.interiors]
    return GeometryBasics.Polygon(exterior, interiors)
end

function _shape_geometry(primitive::DifferenceShape)
    outer = _shape_geometry(primitive.outer)
    interiors = copy(outer.interiors)
    for hole in primitive.holes
        push!(interiors, reverse(_shape_geometry(hole).exterior))
    end
    return GeometryBasics.Polygon(outer.exterior, interiors)
end

function _shape_geometry(shape::RoundedSectorShape)
    return GeometryBasics.Polygon(Point2f.(tessellate(shape)))
end

function _shape_geometry(shape::ShellShape)
    coordinates = tessellate(shape)
    outer = Point2f.(coordinates.outer)
    inner = reverse(Point2f.(coordinates.inner))
    return GeometryBasics.Polygon(outer, [inner])
end

preview_materials(region::PlacedRegion) = (region.source.material,)

function preview_shapes(region::PlacedRegion)
    return PreviewShape[
        PreviewShape(
        _shape_geometry(region.primitive),
        region.source.material,
        region.source.tag
    ),
    ]
end

function _each_material(callback, design::CableDesign)
    for source in design.geometry.regions
        foreach(callback, preview_materials(source))
    end
    return nothing
end

function _each_material(callback, designs::AbstractVector{<:CableDesign})
    for design in designs
        _each_material(callback, design)
    end
    return nothing
end

function _property_ranges(design)
    resistivities = Float64[]
    permeabilities = Float64[]
    permittivities = Float64[]
    _each_material(design) do material
        resistivity = nominal(material.rho)
        permeability = nominal(material.mu_r)
        permittivity = nominal(material.eps_r)
        isfinite(resistivity) && push!(resistivities, resistivity)
        isfinite(permeability) && push!(permeabilities, permeability)
        isfinite(permittivity) && push!(permittivities, permittivity)
    end
    return (;
        rho = isempty(resistivities) ?
              (_RHO_MIN, _RHO_MAX) : extrema(resistivities),
        mu_r = isempty(permeabilities) ?
               (1.0, 300.0) : extrema(permeabilities),
        eps_r = isempty(permittivities) ?
                (1.0, 1000.0) : extrema(permittivities)
    )
end

"""
    material_property_ranges()
    material_property_ranges(design::CableDesign)
    material_property_ranges(designs::AbstractVector{<:CableDesign})

Return physical resistivity, relative permeability, and relative permittivity
ranges without assigning any presentation colour scheme.
"""
function material_property_ranges()
    return (; rho = (_RHO_MIN, _RHO_MAX), mu_r = (1.0, 300.0), eps_r = (1.0, 1000.0))
end

material_property_ranges(design::CableDesign) = _property_ranges(design)
material_property_ranges(designs::AbstractVector{<:CableDesign}) = _property_ranges(designs)
