const SUPPORTED_PRIMITIVES = (
    :line, :scatter, :histogram, :stairs, :heatmap, :polygon, :hline)

function _overlaps(first::GridArea, second::GridArea)
    !isempty(intersect(first.rows, second.rows)) &&
        !isempty(intersect(first.columns, second.columns))
end

function _validate_area(area::GridArea, rows::Int, columns::Int, owner::AbstractString)
    last(area.rows) <= rows && last(area.columns) <= columns || throw(
        ArgumentError("$owner area exceeds its parent grid tracks"),
    )
    return area
end

function _check_sibling_overlap(children)
    length(children) < 2 && return nothing
    for first_index in 1:(length(children) - 1)
        for second_index in (first_index + 1):length(children)
            first_area = children[first_index].area
            second_area = children[second_index].area
            _overlaps(first_area, second_area) && throw(
                ArgumentError("sibling layout areas $first_index and $second_index overlap"),
            )
        end
    end
    return nothing
end

"""
    validate(value)

Validate a layout, page, or completed plot recipe and return it.

# Errors

- `ArgumentError`, `DimensionMismatch`, or `DomainError` when semantic fields,
  data shapes, layout relationships, placements, or logarithmic data are invalid.
"""
function _check_layout(layout::LayoutSpec)
    isempty(layout.grids) && throw(ArgumentError("a layout requires at least one grid"))
    grid_names = getfield.(layout.grids, :name)
    slot_names = getfield.(layout.slots, :name)
    length(unique(grid_names)) == length(grid_names) || throw(
        ArgumentError("layout grid names must be unique"),
    )
    length(unique(slot_names)) == length(slot_names) || throw(
        ArgumentError("layout slot names must be unique"),
    )
    isempty(intersect(grid_names, slot_names)) || throw(
        ArgumentError("layout grid and slot names must be globally unique"),
    )
    roots = filter(grid -> grid.parent === nothing, layout.grids)
    length(roots) == 1 || throw(ArgumentError("a layout requires exactly one root grid"))
    only(roots).area === nothing ||
        throw(ArgumentError("the root grid cannot have an area"))
    grids = Dict(grid.name => grid for grid in layout.grids)
    for grid in layout.grids
        grid.parent === nothing && continue
        haskey(grids, grid.parent) || throw(
            ArgumentError("grid :$(grid.name) references missing parent :$(grid.parent)"),
        )
        grid.area === nothing && throw(
            ArgumentError("nested grid :$(grid.name) requires an area"),
        )
        parent = grids[grid.parent]
        _validate_area(grid.area, length(parent.rows), length(parent.columns), "grid :$(grid.name)")
        visited = Set{Symbol}((grid.name,))
        ancestor = grid.parent
        while ancestor !== nothing
            ancestor in visited &&
                throw(ArgumentError("layout grid hierarchy contains a cycle"))
            push!(visited, ancestor)
            ancestor = grids[ancestor].parent
        end
    end
    for slot in layout.slots
        haskey(grids, slot.parent) || throw(
            ArgumentError("slot :$(slot.name) references missing grid :$(slot.parent)"),
        )
        parent = grids[slot.parent]
        _validate_area(slot.area, length(parent.rows), length(parent.columns), "slot :$(slot.name)")
    end
    for parent in layout.grids
        children = Any[grid for grid in layout.grids if grid.parent === parent.name]
        append!(children, [slot for slot in layout.slots if slot.parent === parent.name])
        _check_sibling_overlap(children)
    end
    return nothing
end

function Validation.rules(::Type{<:LayoutSpec})
    (Validation.OwnerRule(:plot_layout, _check_layout),)
end

function _validate_series(series::SeriesSpec)
    series.kind in SUPPORTED_PRIMITIVES || throw(
        ArgumentError("unsupported PlotBuilder primitive :$(series.kind)"),
    )
    if series.kind in (:line, :scatter, :stairs)
        series.xdata === nothing && throw(ArgumentError(":$(series.kind) requires x data"))
        series.ydata === nothing && throw(ArgumentError(":$(series.kind) requires y data"))
        length(series.xdata) == length(series.ydata) || throw(
            DimensionMismatch(":$(series.kind) x and y data must have equal lengths"),
        )
    elseif series.kind === :histogram
        series.xdata === nothing && throw(ArgumentError(":histogram requires sample data"))
    elseif series.kind === :heatmap
        any(isnothing, (series.xdata, series.ydata, series.zdata)) && throw(
            ArgumentError(":heatmap requires x, y, and z data"),
        )
        size(series.zdata) == (length(series.xdata), length(series.ydata)) || throw(
            DimensionMismatch(":heatmap z data must match x and y dimensions"),
        )
    elseif series.kind === :polygon
        series.zdata === nothing && throw(ArgumentError(":polygon requires geometry data"))
    elseif series.kind === :hline
        series.ydata === nothing && throw(ArgumentError(":hline requires y data"))
    end
    return series
end

function _validate_log_axis(view::ViewSpec, axis::AxisSpec)
    :log10 in axis.allowed_scales || return axis
    found = false
    for series in view.series
        samples = axis.dim === :x ? series.xdata :
                  axis.dim === :y ? series.ydata : series.zdata
        samples === nothing && continue
        for sample in samples
            nominal_value = nominal(sample)
            nominal_value isa Real || continue
            found = true
            lower = nominal_value - abs(standard_uncertainty(sample))
            isfinite(nominal_value) && isfinite(lower) && lower > 0 || throw(
                DomainError(sample, "logarithmic axes require positive finite data and uncertainty bounds"),
            )
        end
    end
    found || throw(
        ArgumentError("axis :$(axis.dim) declares logarithmic scale without plottable data"),
    )
    return axis
end

function _validate_view_limits(view::ViewSpec)
    view.limits === nothing && return view
    view.limits isa Tuple && length(view.limits) == 2 || throw(
        ArgumentError("view limits must be `(xlimits, ylimits)` or `nothing`"),
    )
    for (axis, limits) in zip((view.xaxis, view.yaxis), view.limits)
        axis === nothing && throw(
            ArgumentError("view limits require both x and y axes"),
        )
        limits isa Tuple && length(limits) == 2 || throw(
            ArgumentError("each view limit must be a two-value tuple"),
        )
        lower, upper = limits
        lower isa Real && upper isa Real && isfinite(lower) && isfinite(upper) || throw(
            ArgumentError("view limits must contain finite real values"),
        )
        lower < upper || throw(
            ArgumentError("view limits must be strictly increasing"),
        )
        axis.scale === :log10 && lower <= 0 &&
            throw(
                DomainError(limits, "logarithmic view limits must be positive"),
            )
    end
    return view
end

function _validate_view_aspect(view::ViewSpec)
    view.aspect === nothing && return view
    view.aspect === :data && return view
    view.aspect isa Real && isfinite(view.aspect) && view.aspect > 0 && return view
    throw(ArgumentError("view aspect must be nothing, :data, or a positive finite number"))
end

function _required_slots(page::PageSpec)
    required = Symbol[]
    append!(required, unique(view.placement.slot for view in page.views))
    scale_controls = any(
        axis -> axis !== nothing && :log10 in axis.allowed_scales,
        (view_axis for view in page.views for view_axis in (view.xaxis, view.yaxis))
    )
    (page.controls.reset || page.controls.export_svg || scale_controls) &&
        push!(required, page.controls.slot)
    page.legend.enabled && push!(required, page.legend.slot)
    page.status.enabled && push!(required, page.status.slot)
    append!(required, unique(colorbar.slot for colorbar in page.colorbars))
    return unique(required)
end

function _validate_legend_slot(page::PageSpec)
    page.legend.enabled || return page
    page.legend.overflow === :ellipsis || return page
    slot = only(filter(item -> item.name === page.legend.slot, page.layout.slots))
    parent = only(filter(grid -> grid.name === slot.parent, page.layout.grids))
    tracks = parent.rows[slot.area.rows]
    any(track -> track isa ContentTrack, tracks) && throw(
        ArgumentError(
        "responsive legend slot :$(slot.name) must use fixed or relative row tracks; " *
        "use `overflow=:show_all` for a content-sized legend row",
    ),
    )
    return page
end

function _check_page(page::PageSpec)
    validate(page.layout)
    slots = Set(getfield.(page.layout.slots, :name))
    missing = setdiff(_required_slots(page), collect(slots))
    isempty(missing) || throw(
        ArgumentError("page content references missing layout slots: $(join(missing, ", "))"),
    )
    _validate_legend_slot(page)
    for view in page.views
        foreach(_validate_series, view.series)
        _validate_view_limits(view)
        _validate_view_aspect(view)
        view.xaxis === nothing || view.xaxis.dim === :x ||
            throw(
                ArgumentError("a view x-axis must declare dimension :x"),
            )
        view.yaxis === nothing || view.yaxis.dim === :y ||
            throw(
                ArgumentError("a view y-axis must declare dimension :y"),
            )
        view.zaxis === nothing || view.zaxis.dim === :z ||
            throw(
                ArgumentError("a view z-axis must declare dimension :z"),
            )
        axes = [axis for axis in (view.xaxis, view.yaxis, view.zaxis) if axis !== nothing]
        foreach(axis -> _validate_log_axis(view, axis), axes)
        length(unique(axis.dim for axis in axes)) == length(axes) || throw(
            ArgumentError("view axes must have unique dimensions"),
        )
    end
    view_keys = [view.key for view in page.views if !isempty(view.key)]
    length(unique(view_keys)) == length(view_keys) || throw(
        ArgumentError("page views must have unique nonempty semantic keys"),
    )
    for slot in unique(view.placement.slot for view in page.views)
        placements = [view.placement for view in page.views if view.placement.slot === slot]
        explicit = map(placement -> placement.area !== nothing, placements)
        all(explicit) || all(!, explicit) ||
            throw(
                ArgumentError("views in slot :$slot cannot mix automatic and explicit placement"),
            )
        all(explicit) && _check_sibling_overlap(placements)
    end
    return nothing
end

function Validation.rules(::Type{<:PageSpec})
    (Validation.OwnerRule(:plot_page, _check_page),)
end

function _check_recipe(recipe::PlotRecipe)
    foreach(validate, recipe.figures)
    keys = [page.key for page in recipe.figures if !isempty(page.key)]
    length(unique(keys)) == length(keys) || throw(
        ArgumentError("render pages must have unique nonempty semantic keys"),
    )
    return nothing
end

function Validation.rules(::Type{<:PlotRecipe})
    (Validation.OwnerRule(:plot_recipe, _check_recipe),)
end
