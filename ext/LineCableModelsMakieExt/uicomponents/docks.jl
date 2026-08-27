mutable struct ResponsiveLegend
    legend::Any
    title::Any
    entries::Any
    ellipsis_entry::Any
    capacity::Int
    heights::Dict{Int, Float64}
    fitting::Bool
end

function _visibility_groups(panels)
    groups = Dict{Symbol, Vector{Any}}()
    labels = Dict{Symbol, String}()
    order = Symbol[]
    for panel in panels
        for group in panel.group_order
            haskey(groups, group) || push!(order, group)
            plots = panel.groups[group]
            append!(get!(groups, group, Any[]), plots)
        end
        merge!(labels, panel.group_labels)
    end
    return groups, labels, order
end

function _rich_scientific_label(label::AbstractString)
    matched = match(r"^([+−-]?[0-9]+(?:\.[0-9]+)?)e([+−-]?[0-9]+)$", label)
    matched === nothing && return label
    coefficient, raw_exponent = matched.captures
    exponent = parse(Int, replace(raw_exponent, "−" => "-"))
    formatted_exponent = replace(string(exponent), "-" => "−")
    prefix = coefficient == "1" ? "" : coefficient == "-1" ? "−" : "$coefficient×"
    return Makie.rich(
        prefix,
        "10",
        Makie.superscript(
            formatted_exponent;
            offset = Makie.Vec2f(0.1, 0.0)
        )
    )
end

function _colorbar_ticks(ticks)
    values, labels = ticks
    return values, map(_rich_scientific_label, labels)
end

_makie_alignment(value::Symbol) = value === :stretch ? :center : value

function _colorbars!(
        slot,
        descriptors;
        halign::Symbol = :left,
        valign::Symbol = :top
)
    isempty(descriptors) && return nothing
    grid = GridLayout(
        width = LEGEND_DOCK_WIDTH,
        tellwidth = true,
        halign = _makie_alignment(halign),
        valign = _makie_alignment(valign)
    )
    grid.default_colgap = Fixed(COLORBAR_LABEL_GAP)
    grid.default_rowgap = Fixed(COLORBAR_ROW_GAP)
    slot[] = grid
    colorbars = Any[]
    for (row, descriptor) in enumerate(descriptors)
        Label(
            grid[row, 1],
            descriptor.label;
            halign = :right,
            valign = :center,
            fontsize = COLORBAR_LABEL_SIZE
        )
        push!(colorbars,
            Colorbar(
                grid[row, 2];
                colormap = descriptor.colormap,
                limits = descriptor.limits,
                ticks = _colorbar_ticks(descriptor.ticks),
                label = "",
                labelvisible = false,
                vertical = false,
                height = 14,
                ticklabelsize = COLORBAR_TICK_LABEL_SIZE,
                alignmode = Mixed(left = 0, right = 0),
                tellwidth = false
            ))
    end
    colsize!(grid, 1, Auto(true))
    colsize!(grid, 2, Fixed(COLORBAR_WIDTH))
    return (; grid, colorbars)
end

function _set_legend_capacity!(state::ResponsiveLegend, capacity::Int)
    total = length(state.entries)
    0 <= capacity <= total || throw(BoundsError(state.entries, capacity))
    capacity == state.capacity && return state
    displayed = copy(state.entries[1:capacity])
    capacity < total && push!(displayed, state.ellipsis_entry)
    state.legend.entrygroups[] = [(state.title, displayed)]
    state.capacity = capacity
    return state
end

function _legend_height!(state::ResponsiveLegend, capacity::Int)
    return get!(state.heights, capacity) do
        _set_legend_capacity!(state, capacity)
        height = state.legend.layoutobservables.autosize[][2]
        height === nothing && return 0.0
        return Float64(height)
    end
end

function _fit_legend!(state::ResponsiveLegend, available_height::Real)
    state.fitting && return state
    state.fitting = true
    try
        total = length(state.entries)
        available = max(0.0, Float64(available_height) - LEGEND_HEIGHT_TOLERANCE)
        if _legend_height!(state, total) <= available
            return _set_legend_capacity!(state, total)
        end
        lower = 0
        upper = max(0, total - 1)
        best = 0
        while lower <= upper
            middle = (lower + upper) ÷ 2
            if _legend_height!(state, middle) <= available
                best = middle
                lower = middle + 1
            else
                upper = middle - 1
            end
        end
        return _set_legend_capacity!(state, best)
    finally
        state.fitting = false
    end
end

function _observe_legend!(
        state::ResponsiveLegend,
        slot_grid,
        context::UIContext
)
    observer = on(slot_grid.layoutobservables.computedbbox) do bounding_box
        _fit_legend!(state, bounding_box.widths[2])
        return nothing
    end
    push!(context.observers, observer)
    return state
end

function _legend!(
        slot,
        panels;
        width = nothing,
        overflow::Symbol = :ellipsis
)
    groups, group_labels, group_order = _visibility_groups(panels)
    entries = Any[]
    labels = String[]
    for group in group_order
        haskey(group_labels, group) || continue
        push!(entries, groups[group])
        push!(labels, group_labels[group])
    end
    isempty(entries) && return nothing
    dimensions = width === nothing ? (;) : (; width)
    ellipsis = PolyElement(color = :transparent, strokecolor = :transparent)
    legend_entries = [entry => (; polystrokecolor = :transparent, polystrokewidth = 0)
                      for entry in entries]
    contents = overflow === :ellipsis ? Any[legend_entries..., ellipsis] : legend_entries
    legend_labels = overflow === :ellipsis ? [labels; "(...)"] : labels
    legend = Legend(
        slot,
        contents,
        legend_labels;
        dimensions...,
        tellheight = overflow === :show_all,
        halign = :left,
        valign = :top
    )
    overflow === :show_all && return (; legend, responsive = nothing)
    title, legend_entries = only(legend.entrygroups[])
    complete_entries = copy(legend_entries[1:(end - 1)])
    responsive = ResponsiveLegend(
        legend,
        title,
        complete_entries,
        last(legend_entries),
        -1,
        Dict{Int, Float64}(),
        false
    )
    _set_legend_capacity!(responsive, length(complete_entries))
    return (; legend, responsive)
end

function _window_padding(padding::NTuple{4, <:Real})
    return ntuple(index -> max(OUTER_WINDOW_INSET, padding[index]), 4)
end

function _page_supports_log(panels, dim::Symbol)
    isempty(panels) && return false
    return all(panels) do panel
        specification = dim === :x ? panel.metadata.xaxis : panel.metadata.yaxis
        specification !== nothing && :log10 in specification.allowed_scales
    end
end
