function _system_limits(system, zoom_factor)
    if zoom_factor !== nothing
        zoom_factor isa Real ||
            throw(ArgumentError("zoom_factor must be a positive real value"))
        isfinite(zoom_factor) && zoom_factor > 0 || throw(
            ArgumentError("zoom_factor must be finite and greater than zero"),
        )
    end
    horizontal = Float64[nominal(cable.horz) for cable in system.cables]
    vertical = Float64[nominal(cable.vert) for cable in system.cables]
    radii = Float64[max(
                        nominal(last(cable.design_data.components).conductor_group.r_ex),
                        nominal(last(cable.design_data.components).insulator_group.r_ex)
                    ) for cable in system.cables]
    center_x = isempty(horizontal) ? 0.0 : mean(horizontal)
    center_y = isempty(vertical) ? -1.0 : mean(vertical)
    half_x = isempty(horizontal) ? 1.0 :
             max(maximum(horizontal .+ radii) - center_x, center_x -
                                                          minimum(horizontal .- radii))
    half_y = isempty(vertical) ? 1.0 :
             max(maximum(vertical .+ radii) - center_y, center_y -
                                                        minimum(vertical .- radii))
    base_halfspan = max(half_x, half_y, eps(Float64))
    halfspan = base_halfspan * 1.05 * (zoom_factor === nothing ? 1.5 : Float64(zoom_factor))
    return (
        (center_x - halfspan, center_x + halfspan),
        (center_y - halfspan, center_y + halfspan)
    )
end

function _earth_colorbars(earth_model)
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
    return _colorbar_specs(extrema(resistivities), extrema(permeabilities),
        extrema(permittivities); alpha_value = 0.25)
end

function _system_shapes(system, earth_model, limits, display_legend)
    polygons = PreviewPolygon[]
    references = (
        PreviewReferenceLine([0.0], :air_earth, :black, 1.5),
    )
    if earth_model !== nothing && !earth_model.vertical_layers
        cumulative_depth = 0.0
        fill_minimum = limits[2][1] - 5.0
        fill_horizontal = (limits[1][1] - 5.0, limits[1][2] + 5.0)
        for (index, layer) in enumerate(earth_model.layers[2:end])
            top = cumulative_depth
            bottom = if isinf(layer.thickness)
                fill_minimum
            else
                cumulative_depth -= nominal(layer.thickness)
            end
            material = (;
                rho = layer.rho, eps_r = layer.eps_r, mu_r = layer.mu_r)
            geometry = Point2f[
                (fill_horizontal[1], top),
                (fill_horizontal[2], top),
                (fill_horizontal[2], bottom),
                (fill_horizontal[1], bottom)
            ]
            push!(
                polygons,
                _preview_polygon(
                    geometry,
                    display_legend ? "Earth layer $index" : nothing,
                    Symbol("earth_$index"),
                    _material_color(material; alpha_value = 0.25);
                    stroke = :transparent,
                    width = 0.0
                )
            )
        end
    end
    for cable in system.cables
        append!(
            polygons,
            _design_shapes(
                cable.design_data,
                nominal(cable.horz),
                nominal(cable.vert);
                display_legend = false
            )
        )
    end
    return polygons, references
end

PlotBuilder.dispatch_on(::Type{SystemPreviewPlotDefinition}) = LineCableSystem
function PlotBuilder.input_kwargs(::Type{SystemPreviewPlotDefinition})
    (
        :earth_model,
        :zoom_factor,
        :display_legend,
        :display_id,
        :display_colorbars
    )
end
PlotBuilder.renderer_kwargs(::Type{SystemPreviewPlotDefinition}) = (:size,)
function PlotBuilder.input_defaults(::Type{SystemPreviewPlotDefinition}, ::LineCableSystem)
    (;
        earth_model = nothing,
        zoom_factor = nothing,
        display_legend = true,
        display_id = false,
        display_colorbars = true
    )
end
function PlotBuilder.renderer_defaults(::Type{SystemPreviewPlotDefinition}, ::LineCableSystem)
    (; size = (900, 700))
end

function PlotBuilder.resolve(::Type{SystemPreviewPlotDefinition}, recipe::PlotBuilder.PlotRecipe)
    zoom_factor = recipe.input.zoom_factor
    if zoom_factor !== nothing
        zoom_factor isa Real || throw(
            ArgumentError("zoom_factor must be a positive real value"),
        )
        isfinite(zoom_factor) && zoom_factor > 0 || throw(
            ArgumentError("zoom_factor must be finite and greater than zero"),
        )
    end
    all(name -> getproperty(recipe.input, name) isa Bool,
        (:display_legend, :display_id, :display_colorbars)) || throw(
        ArgumentError("display_legend, display_id, and display_colorbars must be Bool"),
    )
    recipe.renderer.size isa Tuple{Int, Int} || throw(
        ArgumentError("size must be a tuple of two integers"),
    )
    return recipe
end

function PlotBuilder.fetch(::Type{SystemPreviewPlotDefinition}, recipe::PlotBuilder.PlotRecipe)
    system = recipe.object
    limits = _system_limits(system, recipe.input.zoom_factor)
    polygons, references = _system_shapes(
        system,
        recipe.input.earth_model,
        limits,
        recipe.input.display_legend
    )
    title = _system_title(Val(recipe.input.display_id), system)
    colorbars = _system_colorbars(
        Val(recipe.input.display_colorbars),
        recipe.input.earth_model
    )
    identity = (; kind = :system, id = system.system_id)
    payload = PreviewPayload(
        polygons,
        references,
        title,
        limits,
        colorbars,
        PlotBuilder.LegendDefinition(enabled = recipe.input.display_legend),
        identity,
        PlotBuilder.ExportDefinition(
            theme = recipe.renderer.export_theme,
            name = system.system_id,
            open_file = recipe.renderer.open_export
        )
    )
    return PlotBuilder.PlotRecipe(
        SystemPreviewPlotDefinition,
        system,
        merge(recipe.input, (; payload)),
        recipe.renderer
    )
end

function PlotBuilder._composition(
        ::Type{SystemPreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe
)
    return Val(:empty)
end

_system_title(::Val{false}, system) = "Cable system cross-section"
_system_title(::Val{true}, system) = "Cable system cross-section: $(system.system_id)"

_system_colorbars(::Val{false}, earth_model) = ()
_system_colorbars(::Val{true}, earth_model) = _earth_colorbars(earth_model)

function PlotBuilder.default_title(
        ::Type{SystemPreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key,
        view_key
)
    return recipe.input.payload.title
end

function PlotBuilder.default_figsize(
        ::Type{SystemPreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return recipe.renderer.size
end

function PlotBuilder.layout_spec(
        ::Type{SystemPreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return :preview
end

function PlotBuilder.colorbar_specs(
        ::Type{SystemPreviewPlotDefinition},
        mode::Val,
        recipe::PlotBuilder.PlotRecipe,
        page_key
)
    return recipe.input.payload.colorbars
end

function PlotBuilder.legend_spec(
        ::Type{SystemPreviewPlotDefinition}, mode::Val,
        recipe::PlotBuilder.PlotRecipe, page_key)
    return recipe.input.payload.legend
end

function PlotBuilder.page_identity(
        ::Type{SystemPreviewPlotDefinition}, mode::Val,
        recipe::PlotBuilder.PlotRecipe, page_key)
    return recipe.input.payload.key
end

function PlotBuilder.export_spec(
        ::Type{SystemPreviewPlotDefinition}, mode::Val,
        recipe::PlotBuilder.PlotRecipe, page_key, title::AbstractString)
    return recipe.input.payload.export_definition
end
