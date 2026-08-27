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
    resistivity = nominal(material.rho)
    relative_permittivity = nominal(material.eps_r)
    relative_permeability = nominal(material.mu_r)
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
    iszero(inner_radius) && return GeometryBasics.Polygon(outer)
    inner = reverse(_circle_points(inner_radius, xcenter, ycenter; count))
    return GeometryBasics.Polygon(outer, [inner])
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

function _preview_polygon(geometry, label, group, color; stroke = :black, width = 0.5)
    resolved_label = label === nothing ? nothing : String(label)
    return PreviewPolygon(
        geometry,
        resolved_label,
        group,
        color,
        stroke,
        Float64(width)
    )
end

function _annular_preview_shapes(layer, context)
    geometry = _annulus_polygon(
        nominal(layer.r_in),
        nominal(layer.r_ex),
        context.xcenter,
        context.ycenter
    )
    label = context.include_label ? context.label : nothing
    return PreviewPolygon[
        _preview_polygon(
            geometry,
            label,
            context.group,
            _material_color(layer.material_props);
            stroke = :transparent,
            width = 0.0
        ),
    ]
end

_preview_layer_name(layer) = lowercase(string(nameof(typeof(layer))))

function _preview_layer_identities(component_id, layers, role::Symbol)
    names = _preview_layer_name.(layers)
    totals = Dict{String, Int}()
    for name in names
        totals[name] = get(totals, name, 0) + 1
    end

    seen = Dict{String, Int}()
    component_name = uppercasefirst(replace(String(component_id), '_' => ' '))
    identities = Tuple{String, Symbol}[]
    for (index, name) in enumerate(names)
        occurrence = get(seen, name, 0) + 1
        seen[name] = occurrence
        suffix = totals[name] == 1 ? "" : " $occurrence"
        label = "$component_name: $name$suffix"
        key = replace(
            "preview_$(component_id)_$(role)_$index",
            r"[^0-9A-Za-z]+" => "_"
        )
        push!(identities, (label, Symbol(key)))
    end
    return identities
end

function _design_shapes(design, xcenter, ycenter; display_legend::Bool)
    shapes = PreviewPolygon[]
    outer_radius = try
        nominal(design.components[end].insulator_group.r_ex)
    catch
        NaN
    end
    if isfinite(outer_radius) && outer_radius > 0
        push!(shapes,
            _preview_polygon(_circle_points(outer_radius, xcenter, ycenter), nothing,
                :background, :white; stroke = :transparent, width = 0.0))
    end
    for component in design.components
        conductor_layers = component.conductor_group.layers
        conductor_identities = _preview_layer_identities(
            component.id,
            conductor_layers,
            :conductor
        )
        for (layer, (label, group)) in zip(conductor_layers, conductor_identities)
            append!(shapes, preview_shapes(layer, (;
                label,
                group,
                xcenter,
                ycenter,
                include_label = display_legend
            )))
        end
        insulator_layers = component.insulator_group.layers
        insulator_identities = _preview_layer_identities(
            component.id,
            insulator_layers,
            :insulator
        )
        for (layer, (label, group)) in zip(insulator_layers, insulator_identities)
            append!(shapes, preview_shapes(layer, (;
                label,
                group,
                xcenter,
                ycenter,
                include_label = display_legend
            )))
        end
    end
    return shapes
end

function _each_material(callback, design::CableDesign)
    for component in design.components
        for layer in component.conductor_group.layers
            foreach(callback, preview_materials(layer))
        end
        for layer in component.insulator_group.layers
            foreach(callback, preview_materials(layer))
        end
    end
    return nothing
end

function _each_material(callback, designs::AbstractVector{<:CableDesign})
    for design in designs
        _each_material(callback, design)
    end
    return nothing
end

function _property_ranges(design)
    resistivities = Float64[]
    permeabilities = Float64[]
    permittivities = Float64[]
    _each_material(design) do material
        resistivity = nominal(material.rho)
        permeability = nominal(material.mu_r)
        permittivity = nominal(material.eps_r)
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
    function permeability_color(value)
        fraction = _log_fraction(value, 1.0, 300.0)
        tint = _gradient(_MU_COLORS, fraction)
        overlay = RGBA(tint.r, tint.g, tint.b, 0.5 * fraction)
        color = _alpha_composite(RGBA(gray.r, gray.g, gray.b, 1.0), overlay)
        return RGBA(red(color), green(color), blue(color), alpha_value)
    end
    function permittivity_color(value)
        fraction = _log_fraction(value, 1.0, 1000.0)
        tint = _gradient(_EPS_COLORS, fraction)
        overlay = RGBA(tint.r, tint.g, tint.b, 0.6 * fraction)
        color = _alpha_composite(RGBA(dark.r, dark.g, dark.b, 1.0), overlay)
        return RGBA(red(color), green(color), blue(color), alpha_value)
    end
    return (
        PlotBuilder.ColorbarDefinition(
            "ρ [Ω·m]",
            _color_samples(_base_material_color, rho_lower, rho_upper),
            (0.0, 1.0),
            rho_ticks
        ),
        PlotBuilder.ColorbarDefinition(
            "μᵣ",
            _color_samples(permeability_color, mu_lower, mu_upper),
            (0.0, 1.0),
            mu_ticks
        ),
        PlotBuilder.ColorbarDefinition(
            "εᵣ",
            _color_samples(permittivity_color, eps_lower, eps_upper),
            (0.0, 1.0),
            eps_ticks
        ),
    )
end
