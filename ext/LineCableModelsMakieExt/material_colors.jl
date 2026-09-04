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
const _LEAKY_COLORS = [RGB(0.42, 0.55, 0.15), RGB(0.20, 0.24, 0.25)]
const _INSULATOR_COLORS = [RGB(0.20, 0.24, 0.25), RGB(0.78, 0.86, 0.88)]
const _PERFECT_INSULATOR_COLOR = RGB(0.90, 0.94, 0.96)
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
        return _PERFECT_INSULATOR_COLOR
    elseif resistivity <= _RHO_METAL_MAX
        return _gradient(
            _METAL_COLORS,
            _log_fraction(max(resistivity, 1.0e-8), 1.0e-8, _RHO_METAL_MAX)
        )
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
    return _minimum_lightness(_gradient(
        _INSULATOR_COLORS,
        _log_fraction(min(resistivity, _RHO_MAX), _RHO_LEAKY_MAX, _RHO_MAX)
    ))
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

function _material_color(material; alpha::Real = 1.0)
    resistivity = LineCableModels.nominal(material.rho)
    relative_permittivity = LineCableModels.nominal(material.eps_r)
    relative_permeability = LineCableModels.nominal(material.mu_r)
    base = _minimum_lightness(_base_material_color(resistivity))
    lightness = HSL(base).l

    mu_fraction = clamp(
        _log_fraction(max(relative_permeability, 1.0), 1.0, 300.0), 0, 1)
    mu_tint = _gradient(_MU_COLORS, mu_fraction)
    mu_alpha = 0.50 * mu_fraction * (0.6 + 0.4 * (1 - lightness))

    eps_fraction = clamp(
        _log_fraction(max(relative_permittivity, 1.0), 1.0, 1000.0), 0, 1)
    eps_tint = _gradient(_EPS_COLORS, eps_fraction)
    band_weight = resistivity > _RHO_SEMICON_MAX ? 1.0 :
                  (resistivity > _RHO_METAL_MAX ? 0.6 : 0.35)
    eps_alpha = (0.20 + 0.40 * band_weight) * eps_fraction *
                (0.55 + 0.45 * (1 - lightness))

    color = _alpha_composite(
        RGBA(base.r, base.g, base.b, 1.0),
        RGBA(mu_tint.r, mu_tint.g, mu_tint.b, mu_alpha)
    )
    color = _alpha_composite(
        color,
        RGBA(eps_tint.r, eps_tint.g, eps_tint.b, eps_alpha)
    )
    return RGBA(red(color), green(color), blue(color), alpha)
end

function _compact_number(value::Real)
    numeric = Float64(value)
    isfinite(numeric) || return string(numeric)
    compact = @sprintf("%.3g", numeric)
    occursin('e', compact) || return compact
    coefficient, exponent = split(compact, 'e'; limit = 2)
    superscript = replace(
        string(parse(Int, exponent)),
        "-" => "⁻", "0" => "⁰", "1" => "¹", "2" => "²", "3" => "³",
        "4" => "⁴", "5" => "⁵", "6" => "⁶", "7" => "⁷", "8" => "⁸",
        "9" => "⁹"
    )
    return coefficient == "1" ? "10$superscript" : "$coefficient × 10$superscript"
end

function _colorbar_range(property_range)
    lower, upper = Float64.(property_range)
    lower_label = _compact_number(lower)
    upper_label = _compact_number(upper)
    if lower_label == upper_label
        representative = lower > 0 && upper > 0 ? sqrt(lower * upper) : (lower + upper) / 2
        return representative, nextfloat(representative), ([0.5], [lower_label])
    end
    return lower, upper, ([0.0, 1.0], [lower_label, upper_label])
end

function _color_samples(color, lower, upper; count::Int = 256)
    lower_value = max(Float64(lower), floatmin(Float64))
    upper_value = max(Float64(upper), nextfloat(lower_value))
    values = 10.0 .^ range(log10(lower_value), log10(upper_value); length = count)
    return [color(value) for value in values]
end

function LineCableModels.materialcolors(
        property::Symbol,
        property_range = nothing;
        alpha::Real = 1.0
)
    property in (:rho, :mu_r, :eps_r) || throw(ArgumentError(
        "material property must be :rho, :mu_r, or :eps_r",
    ))
    0 <= alpha <= 1 || throw(ArgumentError(
        "alpha must be between zero and one",
    ))
    ranges = LineCableModels.DataModel.material_property_ranges()
    resolved_range = property_range === nothing ? getproperty(ranges, property) :
                     property_range
    resolved_range isa Tuple && length(resolved_range) == 2 &&
    all(value -> value isa Real, resolved_range) || throw(ArgumentError(
        "a material property range must be a tuple of two real values",
    ))

    gray = RGB(0.5, 0.5, 0.5)
    dark = RGB(0.1, 0.1, 0.1)
    permeability_color(value) = begin
        fraction = _log_fraction(value, 1.0, 300.0)
        tint = _gradient(_MU_COLORS, fraction)
        _alpha_composite(
            RGBA(gray.r, gray.g, gray.b, 1.0),
            RGBA(tint.r, tint.g, tint.b, 0.5 * fraction)
        )
    end
    permittivity_color(value) = begin
        fraction = _log_fraction(value, 1.0, 1000.0)
        tint = _gradient(_EPS_COLORS, fraction)
        _alpha_composite(
            RGBA(dark.r, dark.g, dark.b, 1.0),
            RGBA(tint.r, tint.g, tint.b, 0.6 * fraction)
        )
    end
    label, color, bounded_range = if property === :rho
        "ρ [Ω·m]", _base_material_color, resolved_range
    elseif property === :mu_r
        "μᵣ", permeability_color, max.(resolved_range, 1.0)
    else
        "εᵣ", permittivity_color, max.(resolved_range, 1.0)
    end
    lower, upper, ticks = _colorbar_range(bounded_range)
    colors = _color_samples(color, lower, upper)
    colors = RGBA[RGBA(red(value), green(value), blue(value), alpha) for value in colors]
    return (; label, colormap = colors, limits = (0.0, 1.0), ticks)
end

function _material_schemes(ranges; kwargs...)
    return Tuple(
        LineCableModels.materialcolors(
            property,
            getproperty(ranges, property);
            kwargs...
        ) for property in (:rho, :mu_r, :eps_r)
    )
end
