"""
    AbstractPlotDefinition

Supertype for renderer-independent PlotBuilder recipe identifiers.
"""
abstract type AbstractPlotDefinition end

const LEGEND_OVERFLOW_MODES = (:ellipsis, :show_all)

"""
$(TYPEDEF)

Control the standard legend supplied by the plotting shell.

$(TYPEDFIELDS)
"""
struct LegendDefinition
    "Whether the shell renders a legend."
    enabled::Bool
    "Whether legend entries control the visibility of registered plots."
    interactive::Bool
    "How entries exceeding the available height are handled."
    overflow::Symbol
    function LegendDefinition(
            enabled::Bool,
            interactive::Bool,
            overflow::Symbol
    )
        overflow in LEGEND_OVERFLOW_MODES || throw(
            ArgumentError("legend overflow must be :ellipsis or :show_all"),
        )
        return new(enabled, interactive, overflow)
    end
end

function LegendDefinition(;
        enabled::Bool = true,
        interactive::Bool = true,
        overflow::Symbol = :ellipsis
)
    return LegendDefinition(enabled, interactive, overflow)
end

"""
$(TYPEDEF)

Define one color scale displayed in the standard side dock.

$(TYPEDFIELDS)
"""
struct ColorbarDefinition{C, T}
    "Displayed scale label."
    label::String
    "Ordered colors used by the scale."
    colormap::C
    "Finite, strictly increasing color limits."
    limits::Tuple{Float64, Float64}
    "Tick positions and their displayed labels."
    ticks::T
end

function ColorbarDefinition(
        label::AbstractString,
        colormap,
        limits::Tuple{<:Real, <:Real},
        ticks
)
    lower, upper = Float64.(limits)
    isfinite(lower) && isfinite(upper) && lower < upper || throw(
        ArgumentError("colorbar limits must be finite and strictly increasing"),
    )
    ticks isa Tuple && length(ticks) == 2 || throw(
        ArgumentError("colorbar ticks must be a tuple of positions and labels"),
    )
    positions, labels = ticks
    applicable(length, positions) && applicable(length, labels) || throw(
        ArgumentError("colorbar tick positions and labels must be collections"),
    )
    length(positions) == length(labels) || throw(
        DimensionMismatch("colorbar tick positions and labels must have equal lengths"),
    )
    all(position -> position isa Real && isfinite(position), positions) || throw(
        ArgumentError("colorbar tick positions must be finite real values"),
    )
    all(position -> lower <= position <= upper, positions) || throw(
        ArgumentError("colorbar tick positions must lie within the color limits"),
    )
    return ColorbarDefinition(String(label), colormap, (lower, upper), ticks)
end

"""
$(TYPEDEF)

Define the filename and theme used when the standard shell exports SVG.

$(TYPEDFIELDS)
"""
struct ExportDefinition
    "Export theme, either `:default` or `:publication`."
    theme::Symbol
    "Unsanitized base filename."
    name::String
    "Whether to open the completed file with the system application."
    open_file::Bool
end

function ExportDefinition(;
        theme::Symbol = :default,
        name::AbstractString = "linecablemodels_plot",
        open_file::Bool = true
)
    _validate_export_theme(theme)
    isempty(strip(name)) && throw(ArgumentError("export name cannot be empty"))
    return ExportDefinition(theme, String(name), open_file)
end

"Supertype for definition-owned controls added to the standard toolbar."
abstract type AbstractWidgetDefinition end

"Abstract supertype for renderer-independent grid track sizes."
abstract type AbstractTrackSize end

"""
    FixedTrack(value)

Define a nonnegative fixed-size grid track in pixels.
"""
struct FixedTrack <: AbstractTrackSize
    "Track size in pixels."
    value::Float64
    function FixedTrack(value::Real)
        isfinite(value) && value >= 0 || throw(
            ArgumentError("fixed track size must be finite and nonnegative"),
        )
        return new(Float64(value))
    end
end

"""
    RelativeTrack([weight])

Define a grid track that receives a positive share of available space.
"""
struct RelativeTrack <: AbstractTrackSize
    "Relative share of available space."
    weight::Float64
    function RelativeTrack(weight::Real = 1)
        isfinite(weight) && weight > 0 || throw(
            ArgumentError("relative track weight must be finite and positive"),
        )
        return new(Float64(weight))
    end
end

"Define a grid track sized from its rendered content."
struct ContentTrack <: AbstractTrackSize end

"""
    GridArea(rows, columns)

Select positive, one-based row and column spans in a named grid.
"""
struct GridArea
    "Inclusive row span."
    rows::UnitRange{Int}
    "Inclusive column span."
    columns::UnitRange{Int}
    function GridArea(rows::UnitRange{<:Integer}, columns::UnitRange{<:Integer})
        !isempty(rows) && !isempty(columns) || throw(
            ArgumentError("grid area spans cannot be empty"),
        )
        first(rows) > 0 && first(columns) > 0 || throw(
            ArgumentError("grid areas use positive one-based indices"),
        )
        return new(Int(first(rows)):Int(last(rows)), Int(first(columns)):Int(last(columns)))
    end
end

function GridArea(row::Integer, column::Integer)
    GridArea(Int(row):Int(row), Int(column):Int(column))
end
function GridArea(rows::UnitRange{<:Integer}, column::Integer)
    GridArea(rows, Int(column):Int(column))
end
GridArea(row::Integer, columns::UnitRange{<:Integer}) = GridArea(Int(row):Int(row), columns)

"""
    GridSpec(name; parent, area, rows, columns, rowgap, columngap, padding)

Declare one grid in a renderer-independent named layout tree.
"""
struct GridSpec
    "Unique grid name."
    name::Symbol
    "Parent grid name, or `nothing` for the root."
    parent::Union{Nothing, Symbol}
    "Area occupied in the parent grid."
    area::Union{Nothing, GridArea}
    "Row track declarations."
    rows::Vector{AbstractTrackSize}
    "Column track declarations."
    columns::Vector{AbstractTrackSize}
    "Gap between rows in pixels."
    rowgap::Float64
    "Gap between columns in pixels."
    columngap::Float64
    "Left, right, bottom, and top padding in pixels."
    padding::NTuple{4, Float64}
end

function GridSpec(
        name::Symbol;
        parent::Union{Nothing, Symbol} = nothing,
        area::Union{Nothing, GridArea} = nothing,
        rows::AbstractVector{<:AbstractTrackSize} = AbstractTrackSize[RelativeTrack()],
        columns::AbstractVector{<:AbstractTrackSize} = AbstractTrackSize[RelativeTrack()],
        rowgap::Real = 0,
        columngap::Real = 0,
        padding::NTuple{4, <:Real} = (0, 0, 0, 0)
)
    isempty(rows) && throw(ArgumentError("grid rows cannot be empty"))
    isempty(columns) && throw(ArgumentError("grid columns cannot be empty"))
    isfinite(rowgap) && rowgap >= 0 || throw(
        ArgumentError("grid row gap must be finite and nonnegative"),
    )
    isfinite(columngap) && columngap >= 0 || throw(
        ArgumentError("grid column gap must be finite and nonnegative"),
    )
    all(value -> isfinite(value) && value >= 0, padding) || throw(
        ArgumentError("grid padding must be finite and nonnegative"),
    )
    return GridSpec(
        name,
        parent,
        area,
        AbstractTrackSize[rows...],
        AbstractTrackSize[columns...],
        Float64(rowgap),
        Float64(columngap),
        Tuple(Float64.(padding))
    )
end

"""
    SlotSpec(name, parent, area; halign, valign)

Declare a named content destination inside a grid.
"""
struct SlotSpec
    "Unique slot name."
    name::Symbol
    "Parent grid name."
    parent::Symbol
    "Area occupied in the parent grid."
    area::GridArea
    "Horizontal content alignment."
    halign::Symbol
    "Vertical content alignment."
    valign::Symbol
end

function SlotSpec(
        name::Symbol,
        parent::Symbol,
        area::GridArea;
        halign::Symbol = :center,
        valign::Symbol = :center
)
    halign in (:left, :center, :right, :stretch) || throw(
        ArgumentError("slot horizontal alignment must be :left, :center, :right, or :stretch"),
    )
    valign in (:top, :center, :bottom, :stretch) || throw(
        ArgumentError("slot vertical alignment must be :top, :center, :bottom, or :stretch"),
    )
    return SlotSpec(name, parent, area, halign, valign)
end

"""
    LayoutSpec(name, grids, slots)

Define and validate a renderer-independent named grid tree.
"""
struct LayoutSpec
    "Layout identity."
    name::Symbol
    "Named root and nested grids."
    grids::Vector{GridSpec}
    "Named content destinations."
    slots::Vector{SlotSpec}
    function LayoutSpec(name::Symbol, grids::AbstractVector, slots::AbstractVector)
        layout = new(name, GridSpec[grids...], SlotSpec[slots...])
        validate(layout)
        return layout
    end
end

"""
    PlacementSpec([slot], [area])

Assign a view to a named slot with automatic or explicit grid placement.
"""
struct PlacementSpec
    "Destination slot name."
    slot::Symbol
    "Explicit area inside the slot, or `nothing` for automatic placement."
    area::Union{Nothing, GridArea}
end

PlacementSpec(slot::Symbol = :canvas) = PlacementSpec(slot, nothing)

"Declare page-level reset and SVG controls and their destination slot."
struct ControlSpec
    "Whether to show the reset control."
    reset::Bool
    "Whether to show the SVG export control."
    export_svg::Bool
    "Destination toolbar slot."
    slot::Symbol
end

function ControlSpec(; reset::Bool = true, export_svg::Bool = true, slot::Symbol = :toolbar)
    ControlSpec(reset, export_svg, slot)
end

"""
$(TYPEDEF)

Declare legend visibility, interactivity, destination slot, and overflow mode.

$(TYPEDFIELDS)
"""
struct LegendSpec
    "Whether to render a legend."
    enabled::Bool
    "Whether legend entries control series visibility."
    interactive::Bool
    "Destination legend slot."
    slot::Symbol
    "Overflow mode, either `:ellipsis` or `:show_all`."
    overflow::Symbol
end

"""
$(TYPEDSIGNATURES)

Construct a [`LegendSpec`](@ref).

# Keywords

- `enabled`: Render the legend. Default: `true`.
- `interactive`: Let legend entries control series visibility. Default: `true`.
- `slot`: Destination layout slot. Default: `:legend`.
- `overflow`: Use `:ellipsis` to show the largest fitting entry prefix followed
  by `(...)`, or `:show_all` to render every entry. Default: `:ellipsis`.

# Returns

- A validated [`LegendSpec`](@ref).

# Errors

- Throws `ArgumentError` when `overflow` is unsupported.
"""
function LegendSpec(;
        enabled::Bool = true,
        interactive::Bool = true,
        slot::Symbol = :legend,
        overflow::Symbol = :ellipsis
)
    overflow in LEGEND_OVERFLOW_MODES || throw(
        ArgumentError("legend overflow must be :ellipsis or :show_all"),
    )
    return LegendSpec(enabled, interactive, slot, overflow)
end

"""
    ColorbarSpec(label, colormap, limits, ticks; slot=:colorbars)

Declare one renderer-independent colour bar and its destination slot.
"""
struct ColorbarSpec{C, T}
    "Displayed colorbar label."
    label::String
    "Renderer-independent colour-map value."
    colormap::C
    "Finite, strictly increasing colour limits."
    limits::Tuple{Float64, Float64}
    "Tick positions and labels."
    ticks::T
    "Destination colorbar slot."
    slot::Symbol
end

function ColorbarSpec(
        label::AbstractString,
        colormap,
        limits::Tuple{<:Real, <:Real},
        ticks;
        slot::Symbol = :colorbars
)
    lower, upper = Float64.(limits)
    isfinite(lower) && isfinite(upper) && lower < upper || throw(
        ArgumentError("colorbar limits must be finite and strictly increasing"),
    )
    ticks isa Tuple && length(ticks) == 2 || throw(
        ArgumentError("colorbar ticks must be a tuple of positions and labels"),
    )
    positions, labels = ticks
    applicable(length, positions) && applicable(length, labels) || throw(
        ArgumentError("colorbar tick positions and labels must be collections"),
    )
    length(positions) == length(labels) || throw(
        DimensionMismatch("colorbar tick positions and labels must have equal lengths"),
    )
    all(position -> position isa Real && isfinite(position), positions) || throw(
        ArgumentError("colorbar tick positions must be finite real values"),
    )
    all(position -> lower <= position <= upper, positions) || throw(
        ArgumentError("colour-bar tick positions must lie within the colour limits"),
    )
    return ColorbarSpec(String(label), colormap, (lower, upper), ticks, slot)
end

"Declare status-line visibility, initial text, and destination slot."
struct StatusSpec
    "Whether to render a status line."
    enabled::Bool
    "Initial status message."
    initial::String
    "Destination status slot."
    slot::Symbol
end

function StatusSpec(; enabled::Bool = true, initial::AbstractString = "Ready.", slot::Symbol = :status)
    StatusSpec(enabled, String(initial), slot)
end

"Declare the SVG theme, base filename, and automatic-open behaviour."
struct ExportSpec
    "Export theme, either `:default` or `:publication`."
    theme::Symbol
    "Unsanitised base filename."
    name::String
    "Whether to ask the operating system to open the exported file."
    open_file::Bool
end

function ExportSpec(;
        theme::Symbol = :default,
        name::AbstractString = "linecablemodels_plot",
        open_file::Bool = true
)
    _validate_export_theme(theme)
    isempty(strip(name)) && throw(ArgumentError("export name cannot be empty"))
    return ExportSpec(theme, String(name), open_file)
end

const AXIS_SCALES = (:linear, :log10)
const AXIS_RESERVED_ATTRIBUTES = (
    :scale,
    :allowed_scales,
    :exponent,
    :label,
    :xlabel,
    :ylabel,
    :zlabel,
    :xscale,
    :yscale,
    :zscale,
    :xticks,
    :yticks,
    :zticks,
    :xtickformat,
    :ytickformat,
    :ztickformat,
    :title,
    :aspect
)
const SERIES_RESERVED_ATTRIBUTES = (:group, :visible, :label)
const VIEW_RESERVED_ATTRIBUTES = (
    :placement,
    :aspect,
    :limits,
    :xlabel,
    :ylabel,
    :zlabel,
    :xscale,
    :yscale,
    :zscale,
    :xticks,
    :yticks,
    :zticks,
    :xtickformat,
    :ytickformat,
    :ztickformat,
    :title
)

function _reject_reserved(attributes::NamedTuple, reserved::Tuple, owner::AbstractString)
    collisions = Tuple(key for key in keys(attributes) if key in reserved)
    isempty(collisions) || throw(
        ArgumentError("$owner attributes contain reserved semantic keys: $(join(collisions, ", "))"),
    )
    return attributes
end

"""
    AxisSpec(dim, quantity, units, label[, scale]; allowed_scales, exponent, attributes)

Describe one renderer-independent plot axis.
"""
struct AxisSpec{A <: NamedTuple}
    "Axis dimension, one of `:x`, `:y`, or `:z`."
    dim::Symbol
    "Typed quantity identity."
    quantity::Quantity
    "Display units."
    units::UnitExpr
    "Displayed axis label."
    label::String
    "Current scale."
    scale::Symbol
    "Scales that interactive controls may select."
    allowed_scales::Tuple{Vararg{Symbol}}
    "Base-ten display exponent for linear tick labels."
    exponent::Int
    "Validated visual axis attributes."
    attributes::A
end

function AxisSpec(
        dim::Symbol,
        quantity::Quantity,
        units::UnitExpr,
        label::AbstractString,
        scale::Symbol = :linear;
        allowed_scales = (scale,),
        exponent::Integer = 0,
        attributes::NamedTuple = (;)
)
    dim in (:x, :y, :z) || throw(ArgumentError("axis dimension must be :x, :y, or :z"))
    allowed_scales isa Tuple && !isempty(allowed_scales) || throw(
        ArgumentError("allowed axis scales must be a nonempty tuple"),
    )
    all(item -> item isa Symbol && item in AXIS_SCALES, allowed_scales) || throw(
        ArgumentError("axis scales must be :linear or :log10"),
    )
    length(unique(allowed_scales)) == length(allowed_scales) || throw(
        ArgumentError("allowed axis scales cannot contain duplicates"),
    )
    scale in allowed_scales || throw(
        ArgumentError("current axis scale :$scale is not allowed"),
    )
    _reject_reserved(attributes, AXIS_RESERVED_ATTRIBUTES, "AxisSpec")
    return AxisSpec(
        dim,
        quantity,
        units,
        String(label),
        scale,
        allowed_scales,
        Int(exponent),
        attributes
    )
end

"""
    SeriesSpec(kind, xdata, ydata, zdata, label; group, visible, attributes)

Describe one renderer-independent plotting primitive and its data.
"""
struct SeriesSpec{X, Y, Z, A <: NamedTuple}
    "Primitive symbol rendered through `Val` dispatch."
    kind::Symbol
    "Data associated with the x dimension."
    xdata::X
    "Data associated with the y dimension."
    ydata::Y
    "Data associated with the z dimension or geometry."
    zdata::Z
    "Legend label, or `nothing`."
    label::Union{Nothing, String}
    "Visibility-group identity, or `nothing`."
    group::Union{Nothing, Symbol}
    "Initial series visibility."
    visible::Bool
    "Validated visual primitive attributes."
    attributes::A
end

function SeriesSpec(
        kind::Symbol,
        xdata,
        ydata,
        zdata,
        label;
        group::Union{Nothing, Symbol} = nothing,
        visible::Bool = true,
        attributes::NamedTuple = (;)
)
    _reject_reserved(attributes, SERIES_RESERVED_ATTRIBUTES, "SeriesSpec")
    resolved_label = label === nothing ? nothing : String(label)
    return SeriesSpec(kind, xdata, ydata, zdata, resolved_label, group, visible, attributes)
end

"""
    ViewSpec(xaxis, yaxis, zaxis, title, series, key; placement, aspect, limits, attributes)

Describe one plot panel and its placement.
"""
struct ViewSpec{A <: NamedTuple}
    "Specification for the x axis, or `nothing`."
    xaxis::Union{Nothing, AxisSpec}
    "Specification for the y axis, or `nothing`."
    yaxis::Union{Nothing, AxisSpec}
    "Specification for the z axis, or `nothing`."
    zaxis::Union{Nothing, AxisSpec}
    "Displayed panel title."
    title::String
    "Renderer-independent series declarations."
    series::Vector{SeriesSpec}
    "Semantic panel identity."
    key::NamedTuple
    "Named-slot placement."
    placement::PlacementSpec
    "Renderer-independent aspect declaration."
    aspect::Any
    "Explicit axis limits, or `nothing`."
    limits::Any
    "Validated visual panel attributes."
    attributes::A
end

function ViewSpec(
        xaxis,
        yaxis,
        zaxis,
        title,
        series,
        key::NamedTuple;
        placement::PlacementSpec = PlacementSpec(),
        aspect = nothing,
        limits = nothing,
        attributes::NamedTuple = (;)
)
    _reject_reserved(attributes, VIEW_RESERVED_ATTRIBUTES, "ViewSpec")
    return ViewSpec(
        xaxis,
        yaxis,
        zaxis,
        String(title),
        SeriesSpec[series...],
        key,
        placement,
        aspect,
        limits,
        attributes
    )
end

"""
    PageSpec(title, size, key, layout, views; controls, legend, colorbars, status, export_spec)

Describe one complete render page using typed renderer-independent components.
"""
struct PageSpec
    "Displayed page title."
    title::String
    "Figure width and height in pixels."
    size::Tuple{Int, Int}
    "Semantic page identity."
    key::NamedTuple
    "Named layout tree."
    layout::LayoutSpec
    "Plot panels on the page."
    views::Vector{ViewSpec}
    "Interactive control declaration."
    controls::ControlSpec
    "Legend declaration."
    legend::LegendSpec
    "Colorbar declarations."
    colorbars::Vector{ColorbarSpec}
    "Status-line declaration."
    status::StatusSpec
    "SVG export declaration."
    export_spec::ExportSpec
end

function PageSpec(
        title::AbstractString,
        size::Tuple{<:Integer, <:Integer},
        key::NamedTuple,
        layout::LayoutSpec,
        views::AbstractVector;
        controls::ControlSpec = ControlSpec(),
        legend::LegendSpec = LegendSpec(),
        colorbars::AbstractVector = ColorbarSpec[],
        status::StatusSpec = StatusSpec(),
        export_spec::ExportSpec = ExportSpec(name = title)
)
    all(>(0), size) || throw(ArgumentError("page dimensions must be positive"))
    page = PageSpec(
        String(title),
        Tuple(Int.(size)),
        key,
        layout,
        ViewSpec[views...],
        controls,
        legend,
        ColorbarSpec[colorbars...],
        status,
        export_spec
    )
    validate(page)
    return page
end

"""
    PlotRecipe(spec, object, input, renderer, figures)

Store a domain declaration, its resolved options, and its completed validated
plot pages.
"""
struct PlotRecipe{S <: AbstractPlotDefinition, O, I <: NamedTuple, R <: NamedTuple}
    "Plot definition type."
    spec::Type{S}
    "Domain object being plotted."
    object::O
    "Validated semantic recipe options."
    input::I
    "Validated renderer options."
    renderer::R
    "Completed plot pages."
    figures::Vector{PageSpec}
end

function PlotRecipe(
        spec::Type{S}, object::O, input::I, renderer::R,
        figures::AbstractVector = PageSpec[]
) where {S <: AbstractPlotDefinition, O, I <: NamedTuple, R <: NamedTuple}
    recipe = PlotRecipe{S, O, I, R}(spec, object, input, renderer, PageSpec[figures...])
    validate(recipe)
    return recipe
end

function PlotRecipe(spec::Type{S}, figures::AbstractVector) where {S <:
                                                                   AbstractPlotDefinition}
    return PlotRecipe(spec, nothing, (;), (;), figures)
end
