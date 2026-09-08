function _native_cable_title(display::Bool, design)
    display ?
    "Cable design preview: $(design.cable_id)" : "Cable design preview"
end

function _native_cable_colorbars(display::Bool, design)
    display ?
    _material_schemes(LineCableModels.DataModel.material_property_ranges(design)) : ()
end

function _default_preview_legend_group(region, milliken::Bool)
    tag = region.source.tag
    tag === :wire || return tag
    patterns = region.placement.patterns
    bounded = any(patterns) do entry
        entry.pattern isa LineCableModels.DataModel.BoundedPlacement
    end
    milliken && (bounded || isempty(patterns)) && return :stranded_sector
    bounded && return :stranded_core
    return tag
end

function _preview_legend_group(region, legend_group, milliken::Bool)
    tag = region.source.tag
    resolved = if legend_group === nothing
        _default_preview_legend_group(region, milliken)
    elseif legend_group isa Function
        legend_group(region)
    elseif legend_group isa AbstractDict
        get(legend_group, tag, tag)
    else
        throw(ArgumentError(
            "legend_group must be a function, a tag-to-group dictionary, or nothing",
        ))
    end
    resolved isa Union{Symbol, AbstractString} || throw(ArgumentError(
        "preview presentation groups must be symbols or strings",
    ))
    isempty(strip(String(resolved))) && throw(ArgumentError(
        "preview presentation groups cannot be empty",
    ))
    return Symbol(resolved)
end

function _preview_legend_label(group::Symbol, legend_labels)
    fallback = uppercasefirst(replace(String(group), '_' => ' '))
    resolved = if legend_labels === nothing
        fallback
    elseif legend_labels isa Function
        legend_labels(group)
    elseif legend_labels isa AbstractDict
        get(legend_labels, group, fallback)
    else
        throw(ArgumentError(
            "legend_labels must be a function, a group-to-label dictionary, or nothing",
        ))
    end
    resolved isa AbstractString && !isempty(strip(resolved)) || throw(ArgumentError(
        "preview legend labels must be nonempty strings",
    ))
    return String(resolved)
end

function _native_design_shapes(
        design,
        xcenter,
        ycenter;
        display_legend::Bool,
        display_dielectric_pattern::Bool = true,
        legend_group = nothing,
        legend_labels = nothing
)
    offset = LineCableModels.Pose2(xcenter, ycenter, 0)
    milliken = any(region -> region.source.tag === :milliken_fill,
        design.geometry.regions)
    identities = Dict{Symbol, NamedTuple}()
    labelled_groups = Set{Symbol}()
    polygons = PreviewPolygon[]
    for region in design.geometry.regions
        bounded = any(region.placement.patterns) do entry
            entry.pattern isa LineCableModels.DataModel.BoundedPlacement
        end
        enclosure_boundary = any(region.placement.patterns) do entry
            entry.pattern isa LineCableModels.DataModel.EnclosureBoundary
        end
        stroke = enclosure_boundary ? :black :
                 RGB(102 / 255, 109 / 255, 118 / 255)
        width = enclosure_boundary ? 0.8 : (bounded ? 0.5 : 0.45)
        presentation_group = _preview_legend_group(region, legend_group, milliken)
        identity = get!(identities, presentation_group) do
            (;
                label = _preview_legend_label(presentation_group, legend_labels),
                group = Symbol("preview_", presentation_group),
            )
        end
        primitive = LineCableModels.resolve(offset, region.primitive)
        placed = LineCableModels.PlacedRegion(
            region.source,
            primitive,
            region.terminal,
            (patterns = region.placement.patterns,),
            region.paths,
        )
        for shape in LineCableModels.DataModel.preview_shapes(placed)
            label = display_legend && identity.group ∉ labelled_groups ?
                    identity.label : nothing
            push!(polygons, PreviewPolygon(
                shape.geometry,
                label,
                identity.group,
                _material_color(shape.material; pattern = display_dielectric_pattern),
                stroke,
                width,
            ))
        end
        push!(labelled_groups, identity.group)
    end
    return polygons
end

function _native_preview_layout(count::Int, layout)
    count > 0 || throw(ArgumentError(
        "a cable collection preview requires at least one design",
    ))
    if layout === nothing
        columns = ceil(Int, sqrt(count))
        return cld(count, columns), columns
    end
    layout isa Tuple && length(layout) == 2 &&
    all(value -> value isa Integer && !(value isa Bool), layout) || throw(
        ArgumentError("layout must be a tuple of two positive integers or nothing"),
    )
    rows, columns = Int.(layout)
    rows > 0 && columns > 0 || throw(ArgumentError(
        "layout dimensions must be positive",
    ))
    rows * columns >= count || throw(DimensionMismatch(
        "layout provides $(rows * columns) slots for $count cable designs",
    ))
    return rows, columns
end

function _native_system_limits(system, zoom_factor)
    if zoom_factor !== nothing
        zoom_factor isa Real || throw(ArgumentError(
            "zoom_factor must be a positive real value",
        ))
        isfinite(zoom_factor) && zoom_factor > 0 || throw(ArgumentError(
            "zoom_factor must be finite and greater than zero",
        ))
    end
    horizontal = Float64[nominal(position.x) for position in system.positions]
    vertical = Float64[nominal(position.y) for position in system.positions]
    radii = Float64[nominal(outer_radius(design)) for design in system.designs]
    center_x = isempty(horizontal) ? 0.0 : mean(horizontal)
    center_y = isempty(vertical) ? -1.0 : mean(vertical)
    half_x = isempty(horizontal) ? 1.0 :
             max(maximum(horizontal .+ radii) - center_x,
        center_x - minimum(horizontal .- radii))
    half_y = isempty(vertical) ? 1.0 :
             max(maximum(vertical .+ radii) - center_y,
        center_y - minimum(vertical .- radii))
    base_halfspan = max(half_x, half_y, eps(Float64))
    halfspan = base_halfspan * 1.05 * (zoom_factor === nothing ? 1.5 : Float64(zoom_factor))
    return (
        (center_x - halfspan, center_x + halfspan),
        (center_y - halfspan, center_y + halfspan)
    )
end

function _native_earth_colorbars(earth_model)
    earth_model === nothing && return ()
    resistivities = Float64[]
    permeabilities = Float64[]
    permittivities = Float64[]
    for layer in earth_model.layers[2:end]
        push!(resistivities, nominal(layer.rho))
        push!(permeabilities, nominal(layer.mu_r))
        push!(permittivities, nominal(layer.eps_r))
    end
    isempty(resistivities) && return ()
    ranges = (;
        rho = extrema(resistivities),
        mu_r = extrema(permeabilities),
        eps_r = extrema(permittivities)
    )
    return _earth_material_schemes(ranges; alpha = 0.25)
end

function _native_system_shapes(
        system,
        display_legend;
        display_dielectric_pattern::Bool = true,
        legend_group = nothing,
        legend_labels = nothing
)
    polygons = PreviewPolygon[]
    references = (PreviewReference([0.0], :air_earth, :black, 1.5),)
    for (design, position) in zip(system.designs, system.positions)
        append!(polygons,
            _native_design_shapes(
                design,
                nominal(position.x),
                nominal(position.y);
                display_legend = display_legend &&
                                 (legend_group !== nothing || legend_labels !== nothing),
                display_dielectric_pattern,
                legend_group,
                legend_labels
            ))
    end
    return polygons, references
end

function _native_system_title(display::Bool, system)
    display ?
    "Cable system cross-section: $(system.system_id)" : "Cable system cross-section"
end

function _native_system_colorbars(display::Bool, earth_model)
    display ?
    _native_earth_colorbars(earth_model) : ()
end
