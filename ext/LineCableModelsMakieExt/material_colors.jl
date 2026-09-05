const _CONDUCTOR_RHO_RANGE = (1.72e-8, 2.5e-7)
const _SEMICON_RHO_RANGE = (1.0e-2, 1.0e4)
const _DIELECTRIC_EPS_RANGE = (1.0, 10.0)
const _MAGNETIC_MU_RANGE = (1.0, 300.0)
const _EARTH_RHO_RANGE = (1.0e-1, 1.0e4)

const _CONDUCTOR_COLORS = [
    RGB(221 / 255, 228 / 255, 236 / 255),
    RGB(203 / 255, 185 / 255, 160 / 255),
    RGB(198 / 255, 134 / 255, 66 / 255)
]
const _DIELECTRIC_COLORS = [
    RGB(1, 1, 1),
    RGB(248 / 255, 244 / 255, 232 / 255),
    RGB(242 / 255, 230 / 255, 200 / 255),
    RGB(234 / 255, 207 / 255, 154 / 255)
]
const _SEMICON_COLORS = [
    RGB(0.349, 0.392, 0.404),
    RGB(0.604, 0.467, 0.349)
]
const _MAGNETIC_COLORS = [
    RGB(89 / 255, 102 / 255, 216 / 255),
    RGB(214 / 255, 79 / 255, 216 / 255)
]
const _EARTH_COLORS = [
    RGB(0.760, 0.704, 0.590),
    RGB(0.357, 0.286, 0.220)
]

function _oklab_blend(first_color, second_color, value::Real)
    t = clamp(Float64(value), 0.0, 1.0)
    first_oklab = convert(Oklab{Float64}, convert(RGB{Float64}, first_color))
    second_oklab = convert(Oklab{Float64}, convert(RGB{Float64}, second_color))
    mixed = Oklab(
        (1 - t) * first_oklab.l + t * second_oklab.l,
        (1 - t) * first_oklab.a + t * second_oklab.a,
        (1 - t) * first_oklab.b + t * second_oklab.b,
    )
    rgb = convert(RGB{Float64}, mixed)
    return RGB(
        clamp(red(rgb), 0, 1),
        clamp(green(rgb), 0, 1),
        clamp(blue(rgb), 0, 1),
    )
end

function _gradient(colors, value::Real)
    t = clamp(Float64(value), 0.0, 1.0)
    position = t * (length(colors) - 1)
    index = clamp(floor(Int, position) + 1, 1, length(colors) - 1)
    fraction = position - (index - 1)
    return _oklab_blend(colors[index], colors[index + 1], fraction)
end

function _log_fraction(value, lower, upper)
    clamped = clamp(Float64(value), Float64(lower), Float64(upper))
    return (log10(clamped) - log10(lower)) / (log10(upper) - log10(lower))
end

function _linear_fraction(value, lower, upper)
    clamped = clamp(Float64(value), Float64(lower), Float64(upper))
    return (clamped - lower) / (upper - lower)
end

_conductor_color(resistivity::Real) = _gradient(
    _CONDUCTOR_COLORS,
    _log_fraction(resistivity, _CONDUCTOR_RHO_RANGE...)
)

_dielectric_color(relative_permittivity::Real) = _gradient(
    _DIELECTRIC_COLORS,
    _linear_fraction(max(relative_permittivity, 1.0), _DIELECTRIC_EPS_RANGE...)
)

_semicon_color(resistivity::Real) = _gradient(
    _SEMICON_COLORS,
    _log_fraction(resistivity, _SEMICON_RHO_RANGE...)
)

_earth_color(resistivity::Real) = _gradient(
    _EARTH_COLORS,
    _log_fraction(resistivity, _EARTH_RHO_RANGE...)
)

_material_base(material::Material) = _material_base(Val(material.kind), material)
_material_base(layer::EarthLayer) = _earth_color(nominal(layer.rho))
_material_base(::Val{:conductor}, material::Material) = _conductor_color(nominal(material.rho))
_material_base(::Val{:insulator}, material::Material) = _dielectric_color(nominal(material.eps_r))
_material_base(::Val{:semicon}, material::Material) = _semicon_color(nominal(material.rho))
function _material_base(::Val{Kind}, material) where {Kind}
    throw(ArgumentError(
        "preview coloring is not defined for material kind :$Kind"
    ))
end

function _magnetic_overlay(base, relative_permeability::Real)
    fraction = _log_fraction(
        max(relative_permeability, 1.0), _MAGNETIC_MU_RANGE...
    )
    iszero(fraction) && return RGB(base)
    tint = _gradient(_MAGNETIC_COLORS, fraction)
    return _oklab_blend(base, tint, 0.35fraction^0.7)
end

function _material_color(material::Union{Material, EarthLayer}; alpha::Real = 1.0)
    relative_permeability = LineCableModels.nominal(material.mu_r)
    base = _material_base(material)
    color = _magnetic_overlay(base, relative_permeability)
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

function _color_samples(
        color,
        lower,
        upper;
        count::Int = 256,
        logarithmic::Bool = true,
)
    lower_value = max(Float64(lower), floatmin(Float64))
    upper_value = max(Float64(upper), nextfloat(lower_value))
    values = logarithmic ?
             10.0 .^ range(log10(lower_value), log10(upper_value); length = count) :
             range(lower_value, upper_value; length = count)
    return [color(value) for value in values]
end

function _bounded_property_range(property_range, supported_range)
    lower, upper = Float64.(property_range)
    lower <= upper || throw(ArgumentError(
        "a material property range must be ordered from lower to upper"
    ))
    return (
        clamp(lower, supported_range...),
        clamp(upper, supported_range...)
    )
end

function _material_scheme(
        label,
        color,
        property_range,
        supported_range,
        alpha;
        logarithmic::Bool = true,
)
    bounded_range = _bounded_property_range(property_range, supported_range)
    lower, upper, ticks = _colorbar_range(bounded_range)
    colors = _color_samples(color, lower, upper; logarithmic)
    colors = RGBA[RGBA(red(value), green(value), blue(value), alpha) for value in colors]
    return (; label, colormap = colors, limits = (0.0, 1.0), ticks)
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

    neutral = RGB(0.94, 0.94, 0.94)
    permeability_color(value) = _magnetic_overlay(neutral, value)
    label, color, supported_range, logarithmic = if property === :rho
        "ρ [Ω·m]", _conductor_color, _CONDUCTOR_RHO_RANGE, true
    elseif property === :mu_r
        "μᵣ", permeability_color, _MAGNETIC_MU_RANGE, true
    else
        "εᵣ", _dielectric_color, _DIELECTRIC_EPS_RANGE, false
    end
    return _material_scheme(
        label,
        color,
        resolved_range,
        supported_range,
        alpha;
        logarithmic,
    )
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


function _earth_material_schemes(ranges; alpha::Real = 0.25)
    return (
        _material_scheme(
            "ρ [Ω·m]", _earth_color, ranges.rho, _EARTH_RHO_RANGE, alpha
        ),
        LineCableModels.materialcolors(:mu_r, ranges.mu_r; alpha),
        LineCableModels.materialcolors(:eps_r, ranges.eps_r; alpha)
    )
end
