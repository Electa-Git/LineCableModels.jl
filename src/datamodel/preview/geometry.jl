"""
$(TYPEDEF)

Detached geometry and material identity for one physical cable region.

This is domain data, not a plotting instruction. Renderers choose colours,
strokes, labels, and legend grouping themselves.

$(TYPEDFIELDS)
"""
struct PreviewShape{G, M}
    "Closed two-dimensional region geometry in metres."
    geometry::G
    "Physical material represented by the region."
    material::M
    "Stable construction tag carried by the source region."
    tag::Symbol
    function PreviewShape(geometry::G, material::M, tag::Symbol) where {G, M}
        return validate(new{G, M}(geometry, material, tag))
    end
end

function validate(shape::PreviewShape)
    points = shape.geometry isa GeometryBasics.Polygon ?
             Iterators.flatten((shape.geometry.exterior, shape.geometry.interiors...)) :
             shape.geometry
    applicable(iterate, points) || throw(ArgumentError(
        "PreviewShape.geometry must be a collection of points",
    ))
    found = false
    for point in points
        found = true
        coordinates = collect(point)
        length(coordinates) >= 2 || throw(ArgumentError(
            "PreviewShape.geometry points require at least two coordinates",
        ))
        all(value -> value isa Real && isfinite(value), coordinates) ||
            throw(ArgumentError(
                "PreviewShape.geometry coordinates must be finite real numbers",
            ))
    end
    found || throw(ArgumentError("PreviewShape.geometry cannot be empty"))
    isempty(String(shape.tag)) && throw(ArgumentError("PreviewShape.tag cannot be empty"))
    return shape
end
