struct CablePreviewPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

"""
$(TYPEDEF)

Prepare one canvas containing a titled preview panel for each cable design in a
collection.
"""
struct CableCollectionPreviewPlotDefinition <: PlotBuilder.AbstractPlotDefinition end

struct SystemPreviewPlotDefinition <: PlotBuilder.AbstractPlotDefinition end
struct MaterialScalePlotDefinition <: PlotBuilder.AbstractPlotDefinition end

"""
$(TYPEDEF)

Store one filled polygon prepared for a cable or system preview.

$(TYPEDFIELDS)
"""
struct PreviewPolygon{G, C, S}
    "Renderer-neutral polygon geometry."
    geometry::G
    "Legend text, or `nothing` when the polygon introduces no entry."
    label::Union{Nothing, String}
    "Identity shared by polygons controlled by one legend entry."
    group::Symbol
    "Final fill color after material-property mapping."
    color::C
    "Polygon outline color."
    stroke::S
    "Polygon outline width in pixels."
    width::Float64
    function PreviewPolygon(
            geometry::G,
            label::Union{Nothing, String},
            group::Symbol,
            fill_color::C,
            stroke::S,
            width::Float64
    ) where {G, C, S}
        _validate_preview_geometry(geometry)
        _validate_preview_group(group)
        isfinite(width) && width >= 0 || throw(ArgumentError(
            "preview polygon widths must be finite and nonnegative",
        ))
        return new{G, C, S}(geometry, label, group, fill_color, stroke, width)
    end
end

function PreviewPolygon(geometry, label, group, color, stroke, width::Real)
    resolved_label = label === nothing ? nothing : String(label)
    return PreviewPolygon(
        geometry,
        resolved_label,
        Symbol(group),
        color,
        stroke,
        Float64(width)
    )
end

"""
$(TYPEDEF)

Store one horizontal reference prepared for a system preview.

$(TYPEDFIELDS)
"""
struct PreviewReferenceLine{V, C}
    "Vertical coordinates in meters."
    values::V
    "Legend-group identity."
    group::Symbol
    "Line color."
    color::C
    "Line width in pixels."
    width::Float64
    function PreviewReferenceLine(
            values::V,
            group::Symbol,
            color::C,
            width::Float64
    ) where {V, C}
        applicable(iterate, values) && !isempty(values) || throw(ArgumentError(
            "preview reference values must be a nonempty collection",
        ))
        all(value -> value isa Real && isfinite(nominal(value)), values) || throw(
            ArgumentError("preview reference values must be finite real numbers"),
        )
        _validate_preview_group(group)
        isfinite(width) && width >= 0 || throw(ArgumentError(
            "preview reference widths must be finite and nonnegative",
        ))
        return new{V, C}(values, group, color, width)
    end
end

function PreviewReferenceLine(values, group, color, width::Real)
    return PreviewReferenceLine(values, Symbol(group), color, Float64(width))
end

"""
$(TYPEDEF)

Store the detached geometry needed to redraw one preview.

$(TYPEDFIELDS)
"""
struct PreviewPayload{P, R, L, S}
    "Prepared filled polygons in draw order."
    polygons::P
    "Prepared horizontal references in draw order."
    references::R
    "Explicit `(x, y)` limits, or `nothing` for fitted limits."
    limits::L
    "Captured runtime state used for current-state SVG replay."
    runtime::S
    function PreviewPayload(
            polygons::P,
            references::R,
            limits::L,
            runtime::S
    ) where {P, R, L, S}
        _validate_preview_entries(polygons, PreviewPolygon, "polygons")
        _validate_preview_entries(references, PreviewReferenceLine, "references")
        _validate_preview_limits(limits)
        _validate_preview_runtime(runtime, polygons, references)
        return new{P, R, L, S}(polygons, references, limits, runtime)
    end
end

function _validate_preview_entries(entries, type, name::AbstractString)
    applicable(iterate, entries) || throw(ArgumentError(
        "preview $name must be a collection",
    ))
    all(entry -> entry isa type, entries) || throw(ArgumentError(
        "preview $name must contain $(nameof(type)) values",
    ))
    return entries
end

function _validate_preview_group(group::Symbol)
    isempty(String(group)) && throw(ArgumentError(
        "preview group identities cannot be empty",
    ))
    return group
end

_preview_geometry_points(geometry::GeometryBasics.Polygon) =
    Iterators.flatten((geometry.exterior, geometry.interiors...))
_preview_geometry_points(geometry) = geometry

function _validate_preview_geometry(geometry)
    points = _preview_geometry_points(geometry)
    applicable(iterate, points) || throw(ArgumentError(
        "preview geometry must be a collection of points",
    ))
    found = false
    for point in points
        found = true
        applicable(iterate, point) || throw(ArgumentError(
            "preview geometry entries must be coordinate points",
        ))
        coordinates = collect(point)
        length(coordinates) >= 2 || throw(ArgumentError(
            "preview geometry points require at least two coordinates",
        ))
        all(value -> value isa Real && isfinite(nominal(value)), coordinates) || throw(
            ArgumentError("preview geometry coordinates must be finite real numbers"),
        )
    end
    found || throw(ArgumentError("preview geometry cannot be empty"))
    return geometry
end

function _validate_preview_limits(limits)
    limits === nothing && return limits
    limits isa Tuple && length(limits) == 2 || throw(ArgumentError(
        "preview limits must contain x and y intervals or be nothing",
    ))
    for interval in limits
        interval isa Tuple && length(interval) == 2 || throw(ArgumentError(
            "each preview limit must be a two-value interval",
        ))
        lower, upper = interval
        lower isa Real && upper isa Real &&
            isfinite(nominal(lower)) && isfinite(nominal(upper)) || throw(
            ArgumentError("preview limits must be finite real numbers"),
        )
        lower < upper || throw(ArgumentError(
            "preview limit intervals must be strictly increasing",
        ))
    end
    return limits
end

function _validate_preview_runtime(runtime, polygons, references)
    runtime === nothing && return runtime
    runtime isa NamedTuple && keys(runtime) == (:hidden_groups, :current_limits) || throw(
        ArgumentError(
            "preview runtime must contain hidden_groups and current_limits",
        ),
    )
    runtime.hidden_groups isa Tuple &&
        all(group -> group isa Symbol, runtime.hidden_groups) || throw(
        ArgumentError("preview hidden groups must be a tuple of symbols"),
    )
    groups = Set((entry.group for entry in (polygons..., references...)))
    all(group -> group in groups, runtime.hidden_groups) || throw(ArgumentError(
        "preview runtime contains a group absent from its geometry",
    ))
    _validate_preview_limits(runtime.current_limits)
    return runtime
end
