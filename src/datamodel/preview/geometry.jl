"""
$(TYPEDEF)

Detached geometry and material identity for one physical cable region.

Renderers choose colors, strokes, labels, and legend grouping for each region.

$(TYPEDFIELDS)
"""
struct PreviewShape{G, M}
    "Closed two-dimensional region geometry in meters."
    geometry::G
    "Physical material represented by the region."
    material::M
    "Stable construction tag carried by the source region."
    tag::Symbol
end
