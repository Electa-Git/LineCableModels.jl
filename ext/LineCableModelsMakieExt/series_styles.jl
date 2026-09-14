# Apply styles to the native semantic handles before constructing legends and
# visibility controls. Every recipe uses this through the common plot shell.
function _addon_series_styles!(groups, order, attributes; defaults=nothing, shared=(;))
    dependents = Pair{Makie.Plot,Makie.Plot}[]
    # Nest coincident intervals in insertion order, without moving or sampling
    # their coordinates. Stroke widths stay in screen units on logarithmic axes.
    error_groups = filter(order) do group
        # Only grouped curve/interval series receive defaults. A plotwindow
        # caller's individually styled native primitives remain caller-owned.
        any(handle -> handle isa Makie.Lines, groups[group]) &&
            any(handle -> handle isa Makie.Errorbars, groups[group])
    end
    if length(error_groups) > 1
        for (index, group) in enumerate(error_groups)
            width = (length(error_groups) - index) / (length(error_groups) - 1)
            for handle in groups[group]
                handle isa Makie.Errorbars || continue
                handle.whiskerwidth[] = 4 + 6width
                handle.linewidth[] = 1 + width
            end
        end
    end
    attributes === nothing && defaults === nothing && isempty(shared) && return dependents
    styles = LineCableModels.PlotBuilder._series_attributes(attributes, length(order))
    shared_consumed = Set{Symbol}()
    drawing_order = defaults === nothing ? collect(eachindex(order)) :
        sortperm(collect(eachindex(order)); by=index -> defaults[index].priority)
    for index in drawing_order
        group, overrides = order[index], merge(shared,styles[index])
        automatic = defaults === nothing ? nothing : defaults[index]
        automatic_attributes = automatic === nothing ? (;) : automatic.attributes
        # Error bars supply the glyphs for an uncertainty overlay, regardless of
        # result owner. Explicit native marker overrides remain intentional.
        isempty(error_groups) || (automatic_attributes = merge(automatic_attributes, (marker=nothing,)))
        style = merge(automatic_attributes, overrides)
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
                coordinates = if automatic === nothing || haskey(overrides, :marker)
                    handle[1]
                else
                    _addon_marker_coordinates(handle, automatic.phase;
                        endpoints=automatic.endpoints)
                end
                hollow = automatic !== nothing && automatic.hollow &&
                    !haskey(overrides, :marker)
                marker_attributes = merge(hollow ?
                    (strokecolor=handle.color[], strokewidth=1.5) : (;), marker_attributes)
                points = scatter!(handle.parent, coordinates;
                    color=hollow ? :transparent : handle.color[],
                    visible=handle.visible[], marker_attributes...)
                # Visibility follows the same owner contract as uncertainty
                # bars; the shared shell binds it once for the complete series.
                push!(dependents, points => handle)
                on(handle.parent, handle.color) do value
                    if hollow
                        haskey(overrides, :strokecolor) || (points.strokecolor[] = value)
                    else
                        points.color[] = value
                    end
                end
                push!(markers, points)
                union!(consumed, marker_keys)
            end
        end
        union!(shared_consumed,consumed)
        unused = setdiff(keys(styles[index]), consumed)
        isempty(unused) || throw(ArgumentError(
            "series_attributes for $group are unsupported by its native plots: $(join(unused, ", "))"))
        append!(handles, markers)
    end
    unused = setdiff(keys(shared),shared_consumed)
    isempty(unused) || throw(ArgumentError("native series attributes are unsupported by these plots: $(join(unused, ", "))"))
    return dependents
end

# The mask selects saved points only. Stable catalogue slots, rather than the
# current visible-series count, keep phases unchanged when formulas are hidden.
function _addon_marker_coordinates(handle, phase; endpoints=false)
    slot, count = phase
    axis_scene = handle.parent
    return lift(handle[1], axis_scene.viewport) do points, viewport
        n = length(points)
        n == 0 && return points
        width = max(1.0, Float64(viewport.widths[1]))
        cycles = clamp(floor(Int, width / max(80, 14count)), 1, 8)
        stride = max(count, cld(n, cycles))
        offset = floor(Int, (slot - 1) * stride / count)
        indices = collect((1 + offset):stride:n)
        isempty(indices) && push!(indices, 1 + mod(slot - 1, n))
        if endpoints
            first(indices) == 1 || pushfirst!(indices, 1)
            last(indices) == n || push!(indices, n)
        end
        return points[indices]
    end
end

function _addon_comparison_styles(indices, roles, count)
    shapes = (:rect, :diamond, :dtriangle, :cross, :xcross, :pentagon, :hexagon)
    return Tuple((attributes=(;
            marker=role === :reference ? :circle : role === :default ? :utriangle :
                shapes[mod1(index, length(shapes))],
            markersize=role === :reference ? 11 : 8,
            linestyle=:solid),
        hollow=role === :reference, endpoints=role === :reference, phase=(index, count),
        priority=role === :default ? 2 : role === :reference ? 1 : 0)
        for (index, role) in zip(indices, roles))
end

# A deterministic farthest-point palette in perceptual space. Candidate RGB
# colors are in gamut and avoid very light strokes on the white axis background.
# The prefix is stable: requesting more colors never changes existing identities.
const _ADDON_CURVE_COLORS = RGB{Float64}[RGB(0.18, 0.18, 0.18), RGB(0.0, 0.36, 0.68)]
function _addon_comparison_color(index::Int)
    index > 0 || throw(ArgumentError("series color indices must be positive"))
    if index > length(_ADDON_CURVE_COLORS)
        candidates = [RGB{Float64}(HSV(hue, saturation, value))
            for hue in 0:5:355 for saturation in (0.55, 0.75, 0.95)
            for value in (0.55, 0.7, 0.85)]
        filter!(color -> 0.38 <= convert(Oklab, color).l <= 0.72, candidates)
        coordinates = map(color -> convert(Oklab, color), candidates)
        distance(a, b) = (a.l-b.l)^2 + (a.a-b.a)^2 + (a.b-b.b)^2
        distances = [minimum(distance(color, convert(Oklab, selected))
            for selected in _ADDON_CURVE_COLORS) for color in coordinates]
        while length(_ADDON_CURVE_COLORS) < index
            next = argmax(distances)
            distances[next] > 0 || throw(ArgumentError("too many series for the distinct curve palette"))
            push!(_ADDON_CURVE_COLORS, candidates[next])
            selected = coordinates[next]
            distances .= min.(distances, (distance(color, selected) for color in coordinates))
        end
    end
    return _ADDON_CURVE_COLORS[index]
end
