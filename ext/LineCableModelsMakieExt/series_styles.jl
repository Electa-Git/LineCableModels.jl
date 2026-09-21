# Apply styles to the native semantic handles before constructing legends and
# visibility controls. Every recipe uses this through the common plot shell.
function _addon_series_styles!(groups, order, attributes; defaults=nothing, shared=(;),
        marker_coordinates=nothing)
    dependents = Pair{Makie.Plot,Makie.Plot}[]
    # Nest coincident intervals in insertion order, without moving or sampling
    # their coordinates. Stroke widths stay in screen units on logarithmic axes.
    error_groups = filter(order) do group
        # Only grouped curve/interval series receive defaults. A plotwindow
        # caller's individually styled native primitives remain caller-owned.
        any(handle -> handle isa Makie.Lines, groups[group]) &&
            any(handle -> handle isa Makie.Errorbars, groups[group])
    end
    defaults === nothing || sort!(error_groups;
        by=group -> -defaults[findfirst(==(group),order)].priority)
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
    styles = _series_attributes(attributes, length(order))
    shared_consumed = Set{Symbol}()
    drawing_order = defaults === nothing ? collect(eachindex(order)) :
        sortperm(collect(eachindex(order)); by=index -> defaults[index].priority)
    for index in drawing_order
        group, overrides = order[index], merge(shared,styles[index])
        automatic = defaults === nothing ? nothing : defaults[index]
        automatic_attributes = automatic === nothing ? (;) : automatic.attributes
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
                marker_defaults = Makie.default_theme(handle.parent, Makie.Scatter)
                marker_keys = filter(key -> haskey(marker_defaults, key) &&
                    key ∉ (:color, :visible), keys(style))
                marker_attributes = NamedTuple{marker_keys}(map(key -> getproperty(style, key), marker_keys))
                coordinates = if automatic === nothing || haskey(overrides, :marker)
                    handle[1]
                else
                    marker_coordinates[handle]
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

# Select retained sample indices for both glyph kinds. Stable catalogue slots,
# rather than visible-series count, preserve phases when formulas are filtered.
# Interval positions have priority; short series may have no marker positions.
function _addon_glyph_indices(n, width, phase; endpoints=false, uncertain_indices=Int[],
        errorbar_sampling=:staggered)
    slot, count = phase
    n == 0 && return (markers=Int[], intervals=Int[])
    cycles = clamp(floor(Int, max(1, width) / max(80, 14count)), 1, 8)
    stride = max(isempty(uncertain_indices) ? count : 2count, cld(n, cycles))
    # Candidate markers start inside each slot, leaving deterministic endpoints
    # available to the independent reference role without shifting candidate slots.
    offset = floor(Int, (slot - 0.5) * stride / count)
    markers = collect((1 + offset):stride:n)
    isempty(markers) && push!(markers, 1 + mod(slot - 1, n))
    if endpoints
        push!(markers, 1, n)
        sort!(unique!(markers))
    end
    isempty(uncertain_indices) && return (markers=markers, intervals=Int[])
    errorbar_sampling === :all && return (markers=Int[], intervals=collect(1:n))

    # Interval slots precede marker slots by half a slot in each cycle.
    interval_offset = floor(Int, (slot - (endpoints ? 0.5 : 1)) * stride / count)
    intervals = collect((1 + interval_offset):stride:n)
    filter!(index -> index in uncertain_indices, intervals)
    if isempty(intervals)
        target = 1 + mod(interval_offset, n)
        push!(intervals, uncertain_indices[argmin(abs.(uncertain_indices .- target))])
    end
    filter!(index -> index ∉ intervals, markers)
    if isempty(markers) && n > length(intervals)
        for index in 1:n
            if index ∉ intervals
                push!(markers, index)
                break
            end
        end
    end
    return (; markers, intervals)
end

function _addon_comparison_styles(indices, roles, count)
    shapes = (:rect, :diamond, :dtriangle, :cross, :xcross, :pentagon, :hexagon)
    return Tuple((attributes=(;
            color=role === :reference ? RGB(0.0, 0.0, 0.0) : _addon_comparison_color(index),
            marker=role === :reference ? :circle :
                shapes[mod1(index, length(shapes))],
            markersize=role === :reference ? 11 : 8,
            linestyle=:solid),
        hollow=role === :reference, endpoints=role === :reference,
        phase=role === :reference ? (1,1) : (index,count),
        priority=role === :reference ? 1 : 0)
        for (index, role) in zip(indices, roles))
end

# A deterministic farthest-point palette in perceptual space. Candidate RGB
# colors are in gamut and avoid very light strokes on the white axis background.
# The prefix is stable: requesting more colors never changes existing identities.
# Oklab display criteria: L in [0.38,0.72], chroma >= 0.10, distance
# from white >= 0.30 and from black >= 0.40. These are finite-palette
# engineering choices, not a guarantee for arbitrarily many simultaneous curves.
_addon_color_distance(a,b) = (a.l-b.l)^2 + (a.a-b.a)^2 + (a.b-b.b)^2
function _addon_candidate_color(color)
    lab=convert(Oklab,color)
    return 0.38 <= lab.l <= 0.72 && hypot(lab.a,lab.b)>=0.10 &&
        _addon_color_distance(lab,convert(Oklab,RGB(1.,1.,1.)))>=0.30^2 &&
        _addon_color_distance(lab,convert(Oklab,RGB(0.,0.,0.)))>=0.40^2
end
const _ADDON_CURVE_COLORS = filter(_addon_candidate_color,RGB{Float64}[RGB(0.0,0.36,0.68)])
function _addon_comparison_color(index::Int)
    index > 0 || throw(ArgumentError("series color indices must be positive"))
    if index > length(_ADDON_CURVE_COLORS)
        candidates = [RGB{Float64}(HSV(hue, saturation, value))
            for hue in 0:5:355 for saturation in (0.55, 0.75, 0.95)
            for value in (0.55, 0.7, 0.85)]
        filter!(_addon_candidate_color,candidates)
        coordinates = map(color -> convert(Oklab, color), candidates)
        distance(a,b) = _addon_color_distance(a,b)
        distances = [minimum(distance(color, convert(Oklab, selected))
            for selected in (RGB(0.,0.,0.),_ADDON_CURVE_COLORS...)) for color in coordinates]
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
