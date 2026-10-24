_wire_area(diameter::Real) = (one(diameter) * pi / 4) * diameter^2

const _AWG_BASE = 92
const _D0_MM = 0.127
const _AREA0_MM2 = 0.012668

function awg_to_d_mm(number::Real)
    oftype(float(number), _D0_MM) *
    oftype(float(number), _AWG_BASE) ^
    ((oftype(float(number), 36) - number) / oftype(float(number), 39))
end
function awg_to_area_mm2(number::Real)
    oftype(float(number), _AREA0_MM2) *
    oftype(float(number), _AWG_BASE) ^
    ((oftype(float(number), 36) - number) / oftype(float(number), 19.5))
end

"""Return the AWG number corresponding to `diameter` \\[mm\\]."""
function d_mm_to_awg(diameter::Real)
    oftype(float(diameter), 36) -
    oftype(float(diameter), 39) *
    log(diameter / oftype(float(diameter), _D0_MM)) /
    log(oftype(float(diameter), _AWG_BASE))
end
"""Return the AWG number corresponding to solid metal `area` \\[mm²\\]."""
function area_mm2_to_awg(area::Real)
    oftype(float(area), 36) -
    oftype(float(area), 19.5) *
    log(area / oftype(float(area), _AREA0_MM2)) /
    log(oftype(float(area), _AWG_BASE))
end

function awg_label(number::Integer)
    number == -3 && return "0000 (4/0)"
    number == -2 && return "000 (3/0)"
    number == -1 && return "00 (2/0)"
    number == 0 && return "0 (1/0)"
    return string(number)
end

function awg_sizes(::Type{T}, awg_min::Integer = -3, awg_max::Integer = 40) where {T <: Real}
    awg_min <= awg_max || throw(ArgumentError("awg_min must not exceed awg_max"))
    return [(awg_label(number), convert(T, awg_to_d_mm(number)) / T(1000))
            for number in awg_min:awg_max]
end

awg_sizes(awg_min::Integer = -3, awg_max::Integer = 40) = awg_sizes(Float64, awg_min, awg_max)

"""
Apply a fill factor to solid area to approximate stranded metallic area.
"""
function stranded_area_mm2(number::Real; fill_factor::Real = 0.94)
    factor, area = promote(float(fill_factor), float(awg_to_area_mm2(number)))
    return factor * area
end

const _WIRE_RULES = Tuple{Int, Int, Union{Int, Nothing}}[
    (10, 6, 7), (16, 6, 7), (25, 6, 7), (35, 6, 7),
    (50, 6, 19), (70, 12, 19), (95, 15, 19), (120, 15, 37),
    (150, 15, 37), (185, 30, 37), (240, 30, 37), (300, 30, 61),
    (400, 53, 61), (500, 53, 61), (630, 53, 91), (800, 53, 91),
    (1000, 53, 91)
]

_hex_wires(layers::Int) = 1 + 3layers * (layers - 1)

"""Return wire-count bounds for the requested metal `target_area` \\[mm²\\]."""
function _allowed_wires(target_area::Real)
    for (threshold, minimum, maximum) in _WIRE_RULES
        target_area <= threshold && return minimum, maximum
    end
    return 53, nothing
end

function _allowed_wires(wires::Int, minimum::Int, maximum::Union{Int, Nothing})
    return maximum === nothing ? wires >= minimum : minimum <= wires <= maximum
end
