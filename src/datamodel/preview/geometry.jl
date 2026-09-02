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
        points = geometry isa GeometryBasics.Polygon ?
                 Iterators.flatten((geometry.exterior, geometry.interiors...)) : geometry
        applicable(iterate, points) || throw(ArgumentError(
            "preview geometry must be a collection of points",
        ))
        found = false
        for point in points
            found = true
            coordinates = collect(point)
            length(coordinates) >= 2 || throw(ArgumentError(
                "preview geometry points require at least two coordinates",
            ))
            all(value -> value isa Real && isfinite(nominal(value)), coordinates) ||
                throw(ArgumentError(
                    "preview geometry coordinates must be finite real numbers",
                ))
        end
        found || throw(ArgumentError("preview geometry cannot be empty"))
        isempty(String(tag)) && throw(ArgumentError("preview shape tags cannot be empty"))
        return new{G, M}(geometry, material, tag)
    end
end
