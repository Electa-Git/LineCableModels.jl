const LP_FIG_SIZE = (800, 400)

const METRIC_PREFIX_EXPONENT = Dict(
    :yocto => -24,
    :zepto => -21,
    :atto => -18,
    :femto => -15,
    :pico => -12,
    :nano => -9,
    :micro => -6,
    :milli => -3,
    :centi => -2,
    :deci => -1,
    :base => 0,
    :deca => 1,
    :hecto => 2,
    :kilo => 3,
    :mega => 6,
    :giga => 9,
    :tera => 12,
    :peta => 15,
    :exa => 18,
    :zetta => 21,
    :yotta => 24
)

const METRIC_PREFIX_SYMBOL = Dict(
    :yocto => "y",
    :zepto => "z",
    :atto => "a",
    :femto => "f",
    :pico => "p",
    :nano => "n",
    :micro => "μ",
    :milli => "m",
    :centi => "c",
    :deci => "d",
    :base => "",
    :deca => "da",
    :hecto => "h",
    :kilo => "k",
    :mega => "M",
    :giga => "G",
    :tera => "T",
    :peta => "P",
    :exa => "E",
    :zetta => "Z",
    :yotta => "Y"
)

const DEFAULT_QUANTITY_UNITS = Dict(
    :impedance => :base,
    :admittance => :base,
    :resistance => :base,
    :inductance => :milli,
    :conductance => :base,
    :capacitance => :micro,
    :angle => :base
)

struct UnitSpec
    symbol::String
    per_length::Bool
end

struct ComponentMetadata
    component::Symbol
    quantity::Symbol
    symbol::String
    title::String
    axis_label::String
    unit::UnitSpec
end

function get_description(::SeriesImpedance)
    (
        impedance = "Series impedance",
        resistance = "Series resistance",
        inductance = "Series inductance"
    )
end

get_symbol(::SeriesImpedance) = (
    impedance = "Z",
    resistance = "R",
    inductance = "L"
)

get_unit_symbol(::SeriesImpedance) = (
    impedance = "Ω",
    resistance = "Ω",
    inductance = "H"
)

function get_description(::ShuntAdmittance)
    (
        admittance = "Shunt admittance",
        conductance = "Shunt conductance",
        capacitance = "Shunt capacitance"
    )
end

get_symbol(::ShuntAdmittance) = (
    admittance = "Y",
    conductance = "G",
    capacitance = "C"
)

function get_unit_symbol(::ShuntAdmittance)
    (
        admittance = "S",
        conductance = "S",
        capacitance = "F"
    )
end

parent_kind(::SeriesImpedance) = :series_impedance
parent_kind(::ShuntAdmittance) = :shunt_admittance

metric_exponent(prefix::Symbol) = get(METRIC_PREFIX_EXPONENT, prefix) do
    Base.error("Unsupported metric prefix :$(prefix)")
end

prefix_symbol(prefix::Symbol) = get(METRIC_PREFIX_SYMBOL, prefix) do
    Base.error("Unsupported metric prefix :$(prefix)")
end

quantity_scale(prefix::Symbol) = 10.0 ^ (-metric_exponent(prefix))
length_scale(prefix::Symbol) = 10.0 ^ (metric_exponent(prefix))
frequency_scale(prefix::Symbol) = quantity_scale(prefix)

function unit_text(quantity_prefix::Symbol, base_unit::String)
    ps = prefix_symbol(quantity_prefix)
    return isempty(ps) ? base_unit : string(ps, base_unit)
end

function length_unit_text(prefix::Symbol)
    ps = prefix_symbol(prefix)
    return isempty(ps) ? "m" : string(ps, "m")
end

function composite_unit(
        quantity_prefix::Symbol,
        base_unit::String,
        per_length::Bool,
        length_prefix::Symbol
)
    numerator = unit_text(quantity_prefix, base_unit)
    if per_length
        denominator = length_unit_text(length_prefix)
        return string(numerator, "/", denominator)
    else
        return numerator
    end
end

function frequency_axis_label(prefix::Symbol)
    unit = unit_text(prefix, "Hz")
    return string("frequency [", unit, "]")
end

function normalize_quantity_units(units)
    table = Dict(DEFAULT_QUANTITY_UNITS)
    if units isa Symbol
        for key in keys(table)
            table[key] = units
        end
    elseif units isa NamedTuple
        for (key, val) in pairs(units)
            table[key] = val
        end
    elseif units isa AbstractDict
        for (key, val) in units
            table[key] = val
        end
    elseif units === nothing
        return table
    else
        Base.error("Unsupported quantity unit specification $(typeof(units))")
    end
    return table
end

function resolve_quantity_prefix(quantity::Symbol, units::AbstractDict{Symbol, Symbol})
    return get(units, quantity, get(DEFAULT_QUANTITY_UNITS, quantity, :base))
end

function resolve_conductors(data_dims::NTuple{3, Int}, con)
    nrows, ncols, _ = data_dims
    if con === nothing
        return collect(1:nrows), collect(1:ncols)
    elseif con isa Tuple && length(con) == 2
        isel = collect_indices(con[1], nrows)
        jsel = collect_indices(con[2], ncols)
        return isel, jsel
    else
        Base.error("Conductor selector must be a tuple (i_sel, j_sel)")
    end
end

function collect_indices(sel, n)
    if sel === nothing
        return collect(1:n)
    elseif sel isa Integer
        (1 <= sel <= n) ||
            Base.error("Index $(sel) out of bounds for dimension of size $(n)")
        return [sel]
    elseif sel isa AbstractVector
        indices = collect(Int, sel)
        for idx in indices
            (1 <= idx <= n) ||
                Base.error("Index $(idx) out of bounds for dimension of size $(n)")
        end
        return indices
    elseif sel isa AbstractRange
        indices = collect(sel)
        for idx in indices
            (1 <= idx <= n) ||
                Base.error("Index $(idx) out of bounds for dimension of size $(n)")
        end
        return indices
    elseif sel isa Colon
        return collect(1:n)
    else
        Base.error("Unsupported selector $(sel)")
    end
end

function components_for(
        obj::SeriesImpedance,
        mode::Symbol,
        coord::Symbol;
        per_length::Bool = true
)
    desc = get_description(obj)
    sym = get_symbol(obj)
    units = get_unit_symbol(obj)
    if mode == :ZY
        coord in (:cart, :polar) || Base.error("Unsupported coordinate system $(coord)")
        if coord == :cart
            return ComponentMetadata[
                ComponentMetadata(:real, :impedance, sym.impedance,
                    string(desc.impedance, " – real part"),
                    string("real(", sym.impedance, ")"),
                    UnitSpec(units.impedance, per_length)),
                ComponentMetadata(:imag, :impedance, sym.impedance,
                    string(desc.impedance, " – imaginary part"),
                    string("imag(", sym.impedance, ")"),
                    UnitSpec(units.impedance, per_length))
            ]
        else
            return ComponentMetadata[
                ComponentMetadata(:magnitude, :impedance, sym.impedance,
                    string(desc.impedance, " – magnitude"),
                    string("|", sym.impedance, "|"),
                    UnitSpec(units.impedance, per_length)),
                ComponentMetadata(:angle, :angle, sym.impedance,
                    string(desc.impedance, " – angle"),
                    string("angle(", sym.impedance, ")"),
                    UnitSpec("deg", false))
            ]
        end
    elseif mode == :RLCG
        return ComponentMetadata[
            ComponentMetadata(:resistance, :resistance, sym.resistance,
                desc.resistance,
                sym.resistance,
                UnitSpec(units.resistance, per_length)),
            ComponentMetadata(:inductance, :inductance, sym.inductance,
                desc.inductance,
                sym.inductance,
                UnitSpec(units.inductance, per_length))
        ]
    else
        Base.error("Unsupported mode $(mode)")
    end
end

function components_for(
        obj::ShuntAdmittance,
        mode::Symbol,
        coord::Symbol;
        per_length::Bool = true
)
    desc = get_description(obj)
    sym = get_symbol(obj)
    units = get_unit_symbol(obj)
    if mode == :ZY
        coord in (:cart, :polar) || Base.error("Unsupported coordinate system $(coord)")
        if coord == :cart
            return ComponentMetadata[
                ComponentMetadata(:real, :admittance, sym.admittance,
                    string(desc.admittance, " – real part"),
                    string("real(", sym.admittance, ")"),
                    UnitSpec(units.admittance, per_length)),
                ComponentMetadata(:imag, :admittance, sym.admittance,
                    string(desc.admittance, " – imaginary part"),
                    string("imag(", sym.admittance, ")"),
                    UnitSpec(units.admittance, per_length))
            ]
        else
            return ComponentMetadata[
                ComponentMetadata(:magnitude, :admittance, sym.admittance,
                    string(desc.admittance, " – magnitude"),
                    string("|", sym.admittance, "|"),
                    UnitSpec(units.admittance, per_length)),
                ComponentMetadata(:angle, :angle, sym.admittance,
                    string(desc.admittance, " – angle"),
                    string("angle(", sym.admittance, ")"),
                    UnitSpec("deg", false))
            ]
        end
    elseif mode == :RLCG
        if (coord == :cart || coord == :polar)
            @warn "Ignoring argument :$(coord) for RLCG parameters"
        end
        return ComponentMetadata[
            ComponentMetadata(:conductance, :conductance, sym.conductance,
                desc.conductance,
                sym.conductance,
                UnitSpec(units.conductance, per_length)),
            ComponentMetadata(:capacitance, :capacitance, sym.capacitance,
                desc.capacitance,
                sym.capacitance,
                UnitSpec(units.capacitance, per_length))
        ]
    else
        Base.error("Unsupported mode $(mode)")
    end
end

function component_values(component::Symbol, slice, freqs::Vector{<:Real})
    data = collect(slice)
    if component === :real
        return (real.(data))
    elseif component === :imag
        return (imag.(data))
    elseif component === :magnitude
        return (abs.(data))
    elseif component === :angle
        return rad2deg.((angle.(data)))
    elseif component === :resistance || component === :conductance
        return (real.(data))
    elseif component === :inductance
        imag_part = (imag.(data))
        return reactance_to_l(imag_part, freqs)
    elseif component === :capacitance
        imag_part = (imag.(data))
        return reactance_to_c(imag_part, freqs)
    else
        Base.error("Unsupported component $(component)")
    end
end

function reactance_to_l(imag_part::Vector{<:Real}, freqs::Vector{<:Real})
    result = similar(freqs, promote_type(eltype(imag_part), eltype(freqs)))
    two_pi = 2π
    for idx in eachindex(freqs)
        f = freqs[idx]
        if iszero(f)
            result[idx] = NaN
        else
            result[idx] = imag_part[idx] / (two_pi * f)
        end
    end
    return result
end

function reactance_to_c(imag_part::Vector{<:Real}, freqs::Vector{<:Real})
    result = similar(freqs, promote_type(eltype(imag_part), eltype(freqs)))
    two_pi = 2π
    for idx in eachindex(freqs)
        f = freqs[idx]
        if iszero(f)
            result[idx] = NaN
        else
            result[idx] = imag_part[idx] / (two_pi * f)
        end
    end
    return result
end

function legend_label(symbol::String, i::Int, j::Int)
    return string(symbol, "(", i, ",", j, ")")
end
