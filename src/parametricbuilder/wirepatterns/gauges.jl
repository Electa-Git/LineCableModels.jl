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

function d_mm_to_awg(diameter_mm::Real)
    oftype(float(diameter_mm), 36) -
    oftype(float(diameter_mm), 39) *
    log(diameter_mm / oftype(float(diameter_mm), _D0_MM)) /
    log(oftype(float(diameter_mm), _AWG_BASE))
end
function area_mm2_to_awg(area_mm2::Real)
    oftype(float(area_mm2), 36) -
    oftype(float(area_mm2), 19.5) *
    log(area_mm2 / oftype(float(area_mm2), _AREA0_MM2)) /
    log(oftype(float(area_mm2), _AWG_BASE))
end

function awg_label(number::Integer)
    number == -3 && return "0000 (4/0)"
    number == -2 && return "000 (3/0)"
    number == -1 && return "00 (2/0)"
    number == 0 && return "0 (1/0)"
    return string(number)
end

function awg_sizes(::Type{T}, nmin::Integer = -3, nmax::Integer = 40) where {T <: Real}
    nmin <= nmax || throw(ArgumentError("nmin must not exceed nmax"))
    return [(awg_label(number), convert(T, awg_to_d_mm(number)) / T(1000))
            for number in nmin:nmax]
end

awg_sizes(nmin::Integer = -3, nmax::Integer = 40) = awg_sizes(Float64, nmin, nmax)

"Apply a fill factor to solid area to approximate stranded metallic area."
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

function _allowed_wires(target_mm2::Real)
    for (threshold, minimum, maximum) in _WIRE_RULES
        target_mm2 <= threshold && return minimum, maximum
    end
    return 53, nothing
end

function _allowed_wires(wires::Int, minimum::Int, maximum::Union{Int, Nothing})
    return maximum === nothing ? wires >= minimum : minimum <= wires <= maximum
end
