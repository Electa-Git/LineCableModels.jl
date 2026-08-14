struct CablePreviewPlotSpec <: PlotBuilder.AbstractPlotSpec end
struct SystemPreviewPlotSpec <: PlotBuilder.AbstractPlotSpec end
struct MaterialScalePlotSpec <: PlotBuilder.AbstractPlotSpec end

const _RHO_MIN = 1.0e-9
const _RHO_METAL_MAX = 1.0e-6
const _RHO_SEMIMETAL_MAX = 1.0e-4
const _RHO_SEMICON_MAX = 1.0e3
const _RHO_LEAKY_MAX = 1.0e8
const _RHO_MAX = 1.0e10

const _METAL_COLORS = [
    RGB(0.92, 0.90, 0.86),
    RGB(0.89, 0.89, 0.89),
    RGB(0.86, 0.89, 0.92),
    RGB(0.70, 0.72, 0.75)
]
const _SEMIMETAL_COLORS = [RGB(0.70, 0.72, 0.75), RGB(0.80, 0.75, 0.65)]
const _SEMICON_COLORS = [RGB(1.00, 0.83, 0.40), RGB(0.85, 0.55, 0.18)]
const _LEAKY_COLORS = [RGB(0.42, 0.55, 0.15), RGB(0.13, 0.13, 0.13)]
const _INSULATOR_COLORS = [RGB(0.07, 0.07, 0.07), RGB(0.00, 0.00, 0.00)]
const _MU_COLORS = [RGB(0.20, 0.50, 0.95), RGB(0.56, 0.00, 0.91)]
const _EPS_COLORS = [RGB(0.00, 0.85, 0.70), RGB(0.00, 0.55, 0.90)]

function _gradient(colors, value::Real)
    t = clamp(Float64(value), 0.0, 1.0)
    position = t * (length(colors) - 1)
    index = clamp(floor(Int, position) + 1, 1, length(colors) - 1)
    fraction = position - (index - 1)
    first_color = RGB(colors[index])
    second_color = RGB(colors[index + 1])
    return RGB(
        (1 - fraction) * red(first_color) + fraction * red(second_color),
        (1 - fraction) * green(first_color) + fraction * green(second_color),
        (1 - fraction) * blue(first_color) + fraction * blue(second_color)
    )
end

function _log_fraction(value, lower, upper)
    clamped = clamp(Float64(value), Float64(lower), Float64(upper))
    return (log10(clamped) - log10(lower)) / (log10(upper) - log10(lower))
end

function _minimum_lightness(color::RGB, minimum::Float64 = 0.07)
    hsl = HSL(color)
    return RGB(HSL(hsl.h, hsl.s, max(hsl.l, minimum)))
end

function _base_material_color(resistivity::Real)
    if !isfinite(resistivity)
        return _INSULATOR_COLORS[end]
    elseif resistivity <= _RHO_METAL_MAX
        return _gradient(_METAL_COLORS, _log_fraction(max(resistivity, 1.0e-8), 1.0e-8, _RHO_METAL_MAX))
    elseif resistivity <= _RHO_SEMIMETAL_MAX
        return _gradient(
            _SEMIMETAL_COLORS,
            _log_fraction(resistivity, _RHO_METAL_MAX, _RHO_SEMIMETAL_MAX)
        )
    elseif resistivity <= _RHO_SEMICON_MAX
        return _gradient(
            _SEMICON_COLORS,
            _log_fraction(resistivity, _RHO_SEMIMETAL_MAX, _RHO_SEMICON_MAX)
        )
    elseif resistivity <= _RHO_LEAKY_MAX
        return _gradient(
            _LEAKY_COLORS,
            _log_fraction(resistivity, _RHO_SEMICON_MAX, _RHO_LEAKY_MAX)
        )
    end
    return _minimum_lightness(
        _gradient(
        _INSULATOR_COLORS,
        _log_fraction(min(resistivity, _RHO_MAX), _RHO_LEAKY_MAX, _RHO_MAX)
    ),
    )
end

function _alpha_composite(background::RGBA, foreground::RGBA)
    output_alpha = alpha(foreground) + alpha(background) * (1 - alpha(foreground))
    iszero(output_alpha) && return RGBA(0, 0, 0, 0)
    return RGBA(
        (red(foreground) * alpha(foreground) +
         red(background) * alpha(background) * (1 - alpha(foreground))) / output_alpha,
        (green(foreground) * alpha(foreground) +
         green(background) * alpha(background) * (1 - alpha(foreground))) / output_alpha,
        (blue(foreground) * alpha(foreground) +
         blue(background) * alpha(background) * (1 - alpha(foreground))) / output_alpha,
        output_alpha
    )
end

function _material_color(material; alpha_value::Real = 1.0)
    resistivity = to_nominal(material.rho)
    relative_permittivity = to_nominal(material.eps_r)
    relative_permeability = to_nominal(material.mu_r)
    base = _minimum_lightness(_base_material_color(resistivity))
    lightness = HSL(base).l

    mu_fraction = clamp(_log_fraction(max(relative_permeability, 1.0), 1.0, 300.0), 0, 1)
    mu_tint = _gradient(_MU_COLORS, mu_fraction)
    mu_alpha = 0.50 * mu_fraction * (0.6 + 0.4 * (1 - lightness))

    eps_fraction = clamp(_log_fraction(max(relative_permittivity, 1.0), 1.0, 1000.0), 0, 1)
    eps_tint = _gradient(_EPS_COLORS, eps_fraction)
    band_weight = resistivity > _RHO_SEMICON_MAX ? 1.0 :
                  (resistivity > _RHO_METAL_MAX ? 0.6 : 0.35)
    eps_alpha = (0.20 + 0.40 * band_weight) * eps_fraction * (0.55 + 0.45 * (1 - lightness))

    color = _alpha_composite(RGBA(base.r, base.g, base.b, 1.0), RGBA(mu_tint.r, mu_tint.g, mu_tint.b, mu_alpha))
    color = _alpha_composite(color, RGBA(eps_tint.r, eps_tint.g, eps_tint.b, eps_alpha))
    return RGBA(red(color), green(color), blue(color), alpha_value)
end

_finite_point(point) = isfinite(point[1]) && isfinite(point[2])

function _circle_points(radius, xcenter, ycenter; count::Int = 128)
    angles = range(0, 2π; length = count)
    return filter(
        _finite_point,
        Point2f.(xcenter .+ radius .* cos.(angles), ycenter .+ radius .* sin.(angles))
    )
end

function _annulus_polygon(inner_radius, outer_radius, xcenter, ycenter; count::Int = 256)
    outer = _circle_points(outer_radius, xcenter, ycenter; count)
    iszero(inner_radius) && return Polygon(outer)
    inner = reverse(_circle_points(inner_radius, xcenter, ycenter; count))
    return Polygon(outer, [inner])
end

function _radial_wedge(
        inner_radius, outer_radius, width, center_angle, xcenter, ycenter; count::Int = 32)
    angle_width = iszero(inner_radius) ? 0.0 : width / inner_radius
    outer_angles = range(center_angle - angle_width / 2, center_angle + angle_width / 2; length = count)
    inner_angles = reverse(outer_angles)
    xvalues = vcat(
        xcenter .+ outer_radius .* cos.(outer_angles),
        xcenter .+ inner_radius .* cos.(inner_angles)
    )
    yvalues = vcat(
        ycenter .+ outer_radius .* sin.(outer_angles),
        ycenter .+ inner_radius .* sin.(inner_angles)
    )
    return Point2f.(xvalues, yvalues)
end

function _polygon_series(geometry, label, group, color; stroke = :black, width = 0.5)
    return PlotBuilder.SeriesSpec(
        :polygon,
        nothing,
        nothing,
        geometry,
        label;
        attributes = (; color, strokecolor = stroke, strokewidth = width, group)
    )
end

function _layer_series!(series, layer, label, group, xcenter, ycenter; include_label = true)
    color = _material_color(layer.material_props)
    first_label = include_label ? label : nothing
    if layer isa CircStrands
        wire_radius = to_nominal(layer.radius_wire)
        lay_radius = layer.num_wires == 1 ? 0.0 : to_nominal(layer.r_in)
        coordinates = calc_circstrands_coords(
            layer.num_wires,
            wire_radius,
            lay_radius;
            C = (xcenter, ycenter)
        )
        for (index, (xvalue, yvalue)) in enumerate(coordinates)
            push!(
                series,
                _polygon_series(
                    _circle_points(wire_radius, xvalue, yvalue),
                    index == 1 ? first_label : nothing,
                    group,
                    color
                )
            )
        end
    elseif layer isa RectStrands
        for index in 1:layer.num_wires
            angle = (index - 1) * 2π / layer.num_wires
            geometry = _radial_wedge(
                to_nominal(layer.r_in),
                to_nominal(layer.r_ex),
                to_nominal(layer.width),
                angle,
                xcenter,
                ycenter
            )
            push!(series, _polygon_series(geometry, index == 1 ? first_label : nothing, group, color))
        end
    elseif layer isa Union{Strip, Tubular, Semicon, Insulator}
        geometry = _annulus_polygon(
            to_nominal(layer.r_in),
            to_nominal(layer.r_ex),
            xcenter,
            ycenter
        )
        push!(series, _polygon_series(
            geometry, first_label, group, color; stroke = :transparent, width = 0.0))
    elseif layer isa ConductorGroup
        for (index, nested) in enumerate(layer.layers)
            _layer_series!(series, nested, label, group, xcenter, ycenter;
                include_label = index == 1 && include_label)
        end
    elseif layer isa Sector
        geometry = Point2f[(vertex[1] + xcenter, vertex[2] + ycenter)
                           for vertex in layer.vertices]
        push!(series, _polygon_series(geometry, first_label, group, color))
    elseif layer isa SectorInsulator
        outer = Point2f[(vertex[1] + xcenter, vertex[2] + ycenter)
                        for vertex in layer.outer_vertices]
        inner = reverse(Point2f[(vertex[1] + xcenter, vertex[2] + ycenter)
                                for vertex in layer.inner_sector.vertices])
        push!(series, _polygon_series(Polygon(outer, [inner]), first_label, group, color))
    else
        @warn "unsupported cable-preview layer" layer_type = typeof(layer)
    end
    return series
end

function _design_series(design, xcenter, ycenter; display_legend::Bool)
    series = PlotBuilder.SeriesSpec[]
    outer_radius = try
        to_nominal(design.components[end].insulator_group.r_ex)
    catch
        NaN
    end
    if isfinite(outer_radius) && outer_radius > 0
        push!(series,
            _polygon_series(_circle_points(outer_radius, xcenter, ycenter), nothing,
                :background, :white; stroke = :transparent, width = 0.0))
    end
    function preview_identity(layer, layer_name)
        hasproperty(layer, :material_props) || return layer_name, Symbol(layer_name)
        material = layer.material_props
        rho = to_nominal(material.rho)
        mu_r = to_nominal(material.mu_r)
        eps_r = to_nominal(material.eps_r)
        displayed_rho = isfinite(rho) ? round(Float64(rho); sigdigits = 2) : rho
        label = "$layer_name ρ=$displayed_rho"
        key = replace("$(layer_name)_$(rho)_$(mu_r)_$(eps_r)", r"[^0-9A-Za-z]+" => "_")
        return label, Symbol(key)
    end
    for component in design.components
        for layer in component.conductor_group.layers
            layer_name = lowercase(string(nameof(typeof(layer))))
            label, group = preview_identity(layer, layer_name)
            _layer_series!(series, layer, label, group, xcenter,
                ycenter; include_label = display_legend)
        end
        for layer in component.insulator_group.layers
            layer_name = lowercase(string(nameof(typeof(layer))))
            label, group = preview_identity(layer, layer_name)
            _layer_series!(series, layer, label, group, xcenter,
                ycenter; include_label = display_legend)
        end
    end
    return series
end

function _each_material(callback, design)
    function visit(layer)
        if layer isa ConductorGroup
            foreach(visit, layer.layers)
        elseif hasproperty(layer, :material_props)
            callback(layer.material_props)
        end
    end
    for component in design.components
        foreach(visit, component.conductor_group.layers)
        foreach(visit, component.insulator_group.layers)
    end
    return nothing
end

function _property_ranges(design)
    resistivities = Float64[]
    permeabilities = Float64[]
    permittivities = Float64[]
    _each_material(design) do material
        resistivity = to_nominal(material.rho)
        permeability = to_nominal(material.mu_r)
        permittivity = to_nominal(material.eps_r)
        isfinite(resistivity) && push!(resistivities, resistivity)
        isfinite(permeability) && push!(permeabilities, permeability)
        isfinite(permittivity) && push!(permittivities, permittivity)
    end
    return (
        isempty(resistivities) ? (_RHO_MIN, _RHO_MAX) : extrema(resistivities),
        isempty(permeabilities) ? (1.0, 300.0) : extrema(permeabilities),
        isempty(permittivities) ? (1.0, 1000.0) : extrema(permittivities)
    )
end

function _color_samples(function_value, lower, upper; count::Int = 256)
    lower_value = max(Float64(lower), floatmin(Float64))
    upper_value = max(Float64(upper), nextfloat(lower_value))
    values = 10.0 .^ range(log10(lower_value), log10(upper_value); length = count)
    return [function_value(value) for value in values]
end

function _compact_number(value::Real)
    numeric = Float64(value)
    isfinite(numeric) || return string(numeric)
    return @sprintf("%.3g", numeric)
end

function _colorbar_range(range)
    lower, upper = Float64.(range)
    lower_label = _compact_number(lower)
    upper_label = _compact_number(upper)
    if lower_label == upper_label
        representative = lower > 0 && upper > 0 ? sqrt(lower * upper) : (lower + upper) / 2
        return representative, nextfloat(representative), ([0.5], [lower_label])
    end
    return lower, upper,
    ([0.0, 1.0], [lower_label, upper_label])
end

function _colorbar_specs(rho_range, mu_range, eps_range; alpha_value = 1.0)
    rho_lower, rho_upper, rho_ticks = _colorbar_range(rho_range)
    mu_lower, mu_upper, mu_ticks = _colorbar_range(max.(mu_range, 1.0))
    eps_lower, eps_upper, eps_ticks = _colorbar_range(max.(eps_range, 1.0))
    gray = RGB(0.5, 0.5, 0.5)
    dark = RGB(0.1, 0.1, 0.1)
    return (
        (
            label = "ρ [Ω·m]",
            colormap = _color_samples(_base_material_color, rho_lower, rho_upper),
            limits = (0.0, 1.0),
            ticks = rho_ticks
        ),
        (
            label = "μᵣ",
            colormap = _color_samples(
                value -> begin
                    fraction = _log_fraction(value, 1.0, 300.0)
                    tint = _gradient(_MU_COLORS, fraction)
                    overlay = RGBA(tint.r, tint.g, tint.b, 0.5 * fraction)
                    color = _alpha_composite(RGBA(gray.r, gray.g, gray.b, 1.0), overlay)
                    RGBA(red(color), green(color), blue(color), alpha_value)
                end,
                mu_lower,
                mu_upper
            ),
            limits = (0.0, 1.0),
            ticks = mu_ticks
        ),
        (
            label = "εᵣ",
            colormap = _color_samples(
                value -> begin
                    fraction = _log_fraction(value, 1.0, 1000.0)
                    tint = _gradient(_EPS_COLORS, fraction)
                    overlay = RGBA(tint.r, tint.g, tint.b, 0.6 * fraction)
                    color = _alpha_composite(RGBA(dark.r, dark.g, dark.b, 1.0), overlay)
                    RGBA(red(color), green(color), blue(color), alpha_value)
                end,
                eps_lower,
                eps_upper
            ),
            limits = (0.0, 1.0),
            ticks = eps_ticks
        )
    )
end

function _distance_axes()
    quantity = UnitHandler.QuantityTag{:distance}()
    unit = UnitHandler.units(:base, :meter)
    return (
        PlotBuilder.AxisSpec(:x, quantity, unit, "y [m]", :linear),
        PlotBuilder.AxisSpec(:y, quantity, unit, "z [m]", :linear)
    )
end

PlotBuilder.dispatch_on(::Type{CablePreviewPlotSpec}) = CableDesign
function PlotBuilder.input_kwargs(::Type{CablePreviewPlotSpec})
    (
        :x_offset,
        :y_offset,
        :display_legend,
        :display_id,
        :display_colorbars
    )
end
PlotBuilder.renderer_kwargs(::Type{CablePreviewPlotSpec}) = (:size,)
function PlotBuilder.input_defaults(::Type{CablePreviewPlotSpec}, ::CableDesign)
    (;
        x_offset = 0.0,
        y_offset = 0.0,
        display_legend = true,
        display_id = false,
        display_colorbars = true
    )
end
function PlotBuilder.renderer_defaults(::Type{CablePreviewPlotSpec}, ::CableDesign)
    (; size = (900, 700))
end

function PlotBuilder.resolve_input(::Type{CablePreviewPlotSpec}, recipe::NamedTuple)
    recipe.input.x_offset isa Real || throw(ArgumentError("x_offset must be real"))
    recipe.input.y_offset isa Real || throw(ArgumentError("y_offset must be real"))
    all(name -> getproperty(recipe.input, name) isa Bool,
        (:display_legend, :display_id, :display_colorbars)) || throw(
        ArgumentError("display_legend, display_id, and display_colorbars must be Bool"),
    )
    recipe.renderer.size isa Tuple{Int, Int} || throw(
        ArgumentError("size must be a tuple of two integers"),
    )
    return recipe
end

function PlotBuilder.make_axes(
        ::Type{CablePreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
)
    xaxis, yaxis = _distance_axes()
    return (; xaxis, yaxis, zaxis = nothing)
end

function PlotBuilder.make_series(
        ::Type{CablePreviewPlotSpec},
        mode::Val,
        grouping::Val,
        recipe::NamedTuple,
        page_key,
        view_key,
        axes::NamedTuple
)
    return _design_series(
        recipe.object,
        recipe.input.x_offset,
        recipe.input.y_offset;
        display_legend = recipe.input.display_legend
    )
end

_cable_title(::Val{false}, design) = "Cable design preview"
_cable_title(::Val{true}, design) = "Cable design preview: $(design.cable_id)"

_cable_colorbars(::Val{false}, design) = ()
_cable_colorbars(::Val{true}, design) = _colorbar_specs(_property_ranges(design)...)

function PlotBuilder.default_title(
        ::Type{CablePreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
)
    return _cable_title(Val(recipe.input.display_id), recipe.object)
end

function PlotBuilder.view_key(
        ::Type{CablePreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
)
    return (; kind = :cable)
end

function PlotBuilder.view_attributes(
        ::Type{CablePreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
)
    return (; aspect = :data)
end

function PlotBuilder.default_figsize(
        ::Type{CablePreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key
)
    return recipe.renderer.size
end

function PlotBuilder.figure_layout(
        ::Type{CablePreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key
)
    return :preview
end

function PlotBuilder.page_kwargs(
        ::Type{CablePreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        views::Vector{PlotBuilder.ViewSpec}
)
    return (;
        colorbars = _cable_colorbars(Val(recipe.input.display_colorbars), recipe.object),
        display_legend = recipe.input.display_legend,
        export_name = recipe.object.cable_id,
        configuration = (;
            x_offset = recipe.input.x_offset,
            y_offset = recipe.input.y_offset,
            display_id = recipe.input.display_id,
            display_colorbars = recipe.input.display_colorbars
        )
    )
end

function _system_limits(system, zoom_factor)
    if zoom_factor !== nothing
        zoom_factor isa Real ||
            throw(ArgumentError("zoom_factor must be a positive real value"))
        isfinite(zoom_factor) && zoom_factor > 0 || throw(
            ArgumentError("zoom_factor must be finite and greater than zero"),
        )
    end
    horizontal = Float64[to_nominal(cable.horz) for cable in system.cables]
    vertical = Float64[to_nominal(cable.vert) for cable in system.cables]
    radii = Float64[max(
                        to_nominal(last(cable.design_data.components).conductor_group.r_ex),
                        to_nominal(last(cable.design_data.components).insulator_group.r_ex)
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
        push!(resistivities, to_nominal(layer.base_rho_g))
        push!(permeabilities, to_nominal(layer.base_mur_g))
        push!(permittivities, to_nominal(layer.base_epsr_g))
    end
    isempty(resistivities) && return ()
    return _colorbar_specs(extrema(resistivities), extrema(permeabilities),
        extrema(permittivities); alpha_value = 0.25)
end

function _system_series(system, earth_model, limits, display_legend)
    series = PlotBuilder.SeriesSpec[
        PlotBuilder.SeriesSpec(
        :hline,
        nothing,
        [0.0],
        nothing,
        nothing;
        attributes = (; color = :black, linewidth = 1.5, group = :air_earth)
    ),
    ]
    if earth_model !== nothing && !earth_model.vertical_layers
        cumulative_depth = 0.0
        fill_minimum = limits[2][1] - 5.0
        fill_horizontal = (limits[1][1] - 5.0, limits[1][2] + 5.0)
        for (index, layer) in enumerate(earth_model.layers[2:end])
            top = cumulative_depth
            bottom = if isinf(layer.t)
                fill_minimum
            else
                cumulative_depth -= to_nominal(layer.t)
            end
            material = (;
                rho = layer.base_rho_g, eps_r = layer.base_epsr_g, mu_r = layer.base_mur_g)
            geometry = Point2f[
                (fill_horizontal[1], top),
                (fill_horizontal[2], top),
                (fill_horizontal[2], bottom),
                (fill_horizontal[1], bottom)
            ]
            push!(
                series,
                _polygon_series(
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
            series,
            _design_series(
                cable.design_data,
                to_nominal(cable.horz),
                to_nominal(cable.vert);
                display_legend = false
            )
        )
    end
    return series
end

PlotBuilder.dispatch_on(::Type{SystemPreviewPlotSpec}) = LineCableSystem
function PlotBuilder.input_kwargs(::Type{SystemPreviewPlotSpec})
    (
        :earth_model,
        :zoom_factor,
        :display_legend,
        :display_id,
        :display_colorbars
    )
end
PlotBuilder.renderer_kwargs(::Type{SystemPreviewPlotSpec}) = (:size,)
function PlotBuilder.input_defaults(::Type{SystemPreviewPlotSpec}, ::LineCableSystem)
    (;
        earth_model = nothing,
        zoom_factor = nothing,
        display_legend = true,
        display_id = false,
        display_colorbars = true
    )
end
function PlotBuilder.renderer_defaults(::Type{SystemPreviewPlotSpec}, ::LineCableSystem)
    (; size = (900, 700))
end

function PlotBuilder.resolve_input(::Type{SystemPreviewPlotSpec}, recipe::NamedTuple)
    _system_limits(recipe.object, recipe.input.zoom_factor)
    all(name -> getproperty(recipe.input, name) isa Bool,
        (:display_legend, :display_id, :display_colorbars)) || throw(
        ArgumentError("display_legend, display_id, and display_colorbars must be Bool"),
    )
    recipe.renderer.size isa Tuple{Int, Int} || throw(
        ArgumentError("size must be a tuple of two integers"),
    )
    return recipe
end

function PlotBuilder.make_axes(
        ::Type{SystemPreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
)
    xaxis, yaxis = _distance_axes()
    return (; xaxis, yaxis, zaxis = nothing)
end

function PlotBuilder.make_series(
        ::Type{SystemPreviewPlotSpec},
        mode::Val,
        grouping::Val,
        recipe::NamedTuple,
        page_key,
        view_key,
        axes::NamedTuple
)
    limits = _system_limits(recipe.object, recipe.input.zoom_factor)
    return _system_series(
        recipe.object,
        recipe.input.earth_model,
        limits,
        recipe.input.display_legend
    )
end

_system_title(::Val{false}, system) = "Cable system cross-section"
_system_title(::Val{true}, system) = "Cable system cross-section: $(system.system_id)"

_system_colorbars(::Val{false}, earth_model) = ()
_system_colorbars(::Val{true}, earth_model) = _earth_colorbars(earth_model)

function PlotBuilder.default_title(
        ::Type{SystemPreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
)
    return _system_title(Val(recipe.input.display_id), recipe.object)
end

function PlotBuilder.view_key(
        ::Type{SystemPreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
)
    return (; kind = :system)
end

function PlotBuilder.view_attributes(
        ::Type{SystemPreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
)
    return (;
        aspect = :data,
        limits = _system_limits(recipe.object, recipe.input.zoom_factor)
    )
end

function PlotBuilder.default_figsize(
        ::Type{SystemPreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key
)
    return recipe.renderer.size
end

function PlotBuilder.figure_layout(
        ::Type{SystemPreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key
)
    return :preview
end

function PlotBuilder.page_kwargs(
        ::Type{SystemPreviewPlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        views::Vector{PlotBuilder.ViewSpec}
)
    return (;
        colorbars = _system_colorbars(
            Val(recipe.input.display_colorbars),
            recipe.input.earth_model
        ),
        display_legend = recipe.input.display_legend,
        export_name = recipe.object.system_id,
        configuration = (;
            zoom_factor = recipe.input.zoom_factor,
            display_id = recipe.input.display_id,
            display_colorbars = recipe.input.display_colorbars
        )
    )
end

PlotBuilder.dispatch_on(::Type{MaterialScalePlotSpec}) = Nothing
PlotBuilder.renderer_kwargs(::Type{MaterialScalePlotSpec}) = (:size,)
function PlotBuilder.renderer_defaults(::Type{MaterialScalePlotSpec}, ::Nothing)
    (; size = (800, 400))
end

function PlotBuilder.resolve_input(::Type{MaterialScalePlotSpec}, recipe::NamedTuple)
    recipe.renderer.size isa Tuple{Int, Int} || throw(
        ArgumentError("size must be a tuple of two integers"),
    )
    return recipe
end

function PlotBuilder.grouping_mode(
        ::Type{MaterialScalePlotSpec},
        mode::Val,
        recipe::NamedTuple
)
    return Val(:empty)
end

function PlotBuilder.default_title(
        ::Type{MaterialScalePlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        view_key
)
    return "Material property color scale"
end

function PlotBuilder.default_figsize(
        ::Type{MaterialScalePlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key
)
    return recipe.renderer.size
end

function PlotBuilder.figure_layout(
        ::Type{MaterialScalePlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key
)
    return :material_scale
end

function PlotBuilder.page_kwargs(
        ::Type{MaterialScalePlotSpec},
        mode::Val,
        recipe::NamedTuple,
        page_key,
        views::Vector{PlotBuilder.ViewSpec}
)
    return (;
        colorbars = _colorbar_specs(
            (_RHO_MIN, _RHO_MAX),
            (1.0, 300.0),
            (1.0, 1000.0)
        ),
        display_legend = false,
        export_name = "material_scale",
        controls = PlotBuilder.control_definitions(
            reset = false,
            xlog = false,
            ylog = false,
            legend = false,
            visibility = false,
            zoom = false
        )
    )
end
