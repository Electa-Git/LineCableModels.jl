# Apply styles to the native semantic handles before constructing legends and
# visibility controls. Every recipe uses this through the common plot shell.
function _addon_series_styles!(groups, order, attributes)
    attributes === nothing && return nothing
    styles = LineCableModels.PlotBuilder._series_attributes(attributes, length(order))
    for (group, style) in zip(order, styles)
        isempty(style) && continue
        handles = groups[group]
        consumed = Set{Symbol}()
        markers = Any[]
        for handle in handles
            # Native wrappers such as PlotList forward attributes as properties.
            direct = filter(key -> hasproperty(handle, key), keys(style))
            for key in direct
                getproperty(handle, key)[] = getproperty(style, key)
                push!(consumed, key)
            end
            if handle isa Union{Makie.Lines, Makie.LineSegments} && haskey(style, :marker)
                push!(consumed, :marker)
                style.marker === nothing && continue
                marker_keys = filter(key -> key in Makie.attribute_names(Makie.Scatter) &&
                    key ∉ (:color, :visible), keys(style))
                marker_attributes = NamedTuple{marker_keys}(map(key -> getproperty(style, key), marker_keys))
                points = scatter!(handle.parent, handle[1];
                    color=handle.color, visible=handle.visible, marker_attributes...)
                push!(markers, points)
                union!(consumed, marker_keys)
            end
        end
        unused = setdiff(keys(style), consumed)
        isempty(unused) || throw(ArgumentError(
            "series_attributes for $group are unsupported by its native plots: $(join(unused, ", "))"))
        append!(handles, markers)
    end
    return nothing
end
