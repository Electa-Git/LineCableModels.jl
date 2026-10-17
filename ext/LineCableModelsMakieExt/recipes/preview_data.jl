function _native_cable_title(display::Bool, design)
    display ?
    "Cable design preview: $(design.cable_id)" : "Cable design preview"
end

function _native_cable_colorbars(display::Bool, design)
    display ?
    _material_schemes(DataModel.material_property_ranges(design)) : ()
end

function _default_preview_legend_group(placed, milliken::Bool)
    tag = placed.source.tag
    tag === :wire || return tag
    patterns = placed.placement.patterns
    bounded = any(patterns) do entry
        entry.pattern isa DataModel.BoundedPlacement
    end
    milliken && (bounded || isempty(patterns)) && return :stranded_sector
    bounded && return :stranded_core
    return tag
end

_preview_legend_group(placed, ::Nothing, milliken::Bool) =
    _default_preview_legend_group(placed, milliken)
function _preview_legend_group(placed, legend_group, milliken::Bool)
    resolved = validate(legend_group, DataModel.PlacedRegion)(placed)
    resolved isa Union{Symbol, AbstractString} || throw(ArgumentError(
        "preview presentation groups must be symbols or strings",
    ))
    isempty(strip(String(resolved))) && throw(ArgumentError(
        "preview presentation groups cannot be empty",
    ))
    return Symbol(resolved)
end

# A preview groups its regions by a function of the region.
validate(legend_group::Function, ::Type{DataModel.PlacedRegion}) = legend_group
function validate(legend_group, ::Type{DataModel.PlacedRegion})
    throw(ArgumentError("legend_group must be a function of the region, or nothing; " *
        "got $(typeof(legend_group))"))
end

function _preview_legend_label(group::Symbol, legend_labels)
    fallback = uppercasefirst(replace(String(group), '_' => ' '))
    resolved = legend_labels === nothing ? fallback :
               get(validate(legend_labels, Makie.LegendEntry), group, fallback)
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
    milliken = any(placed -> placed.source.tag === :milliken_fill,
        design.geometry.regions)
    identities = Dict{Symbol, NamedTuple}()
    labelled_groups = Set{Symbol}()
    polygons = PreviewPolygon[]
    for placed in design.geometry.regions
        bounded = any(placed.placement.patterns) do entry
            entry.pattern isa DataModel.BoundedPlacement
        end
        enclosure_boundary = any(placed.placement.patterns) do entry
            entry.pattern isa DataModel.EnclosureBoundary
        end
        stroke = enclosure_boundary ? :black :
                 RGB(102 / 255, 109 / 255, 118 / 255)
        width = enclosure_boundary ? 0.8 : (bounded ? 0.5 : 0.45)
        presentation_group = _preview_legend_group(placed, legend_group, milliken)
        identity = get!(identities, presentation_group) do
            (;
                label = _preview_legend_label(presentation_group, legend_labels),
                group = Symbol("preview_", presentation_group),
            )
        end
        primitive = LineCableModels.resolve(offset, placed.primitive)
        resolved = LineCableModels.PlacedRegion(
            placed.source,
            primitive,
            placed.terminal,
            (patterns = placed.placement.patterns,),
            placed.paths,
        )
        for shape in DataModel.preview_shapes(resolved)
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

# The zoom factor of a system preview: a finite positive real number.
function validate(zoom_factor::Real, ::typeof(LineCableModels.preview))
    isfinite(zoom_factor) && zoom_factor > 0 || throw(ArgumentError(
        "zoom_factor must be finite and greater than zero",
    ))
    return zoom_factor
end
function validate(zoom_factor, ::typeof(LineCableModels.preview))
    throw(ArgumentError(
        "zoom_factor must be a positive real value; got $(typeof(zoom_factor))"))
end

function _native_system_limits(system, zoom_factor)
    zoom_factor === nothing || validate(zoom_factor, LineCableModels.preview)
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
