# --------------------------------------------------------------------------
# Metric prefixes
# --------------------------------------------------------------------------

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

const UNIT_SYMBOL = Dict(
    :ohm => "Ω",
    :henry => "H",
    :farad => "F",
    :siemens => "S",
    :meter => "m",
    :hertz => "Hz",
    :degree => "°",
    :dimensionless => ""
)

@inline function _prefix_exp(prefix::Symbol)
    return get(METRIC_PREFIX_EXPONENT, prefix) do
        throw(ArgumentError("unsupported metric prefix :$prefix"))
    end
end

@inline function _prefix_symbol(prefix::Symbol)
    return get(METRIC_PREFIX_SYMBOL, prefix) do
        throw(ArgumentError("unsupported metric prefix :$prefix"))
    end
end

# --------------------------------------------------------------------------
# Core types
# --------------------------------------------------------------------------

"""
$(TYPEDEF)

Describe one physical unit and its metric prefix.

$(TYPEDFIELDS)
"""
@kwdef struct Unit
    "Unit name, such as `:ohm`, `:meter`, or `:dimensionless`."
    name::Symbol = :dimensionless
    "Metric prefix name, such as `:base`, `:milli`, or `:kilo`."
    prefix::Symbol = :base
end

"""
$(TYPEDEF)

Describe numerator and denominator units, such as Ω/km.

$(TYPEDFIELDS)
"""
@kwdef struct Units
    "Numerator units."
    base::Vector{Unit} = [Unit()]  # dimensionless by default
    "Denominator units."
    per::Vector{Unit} = Unit[]
end

# --------------------------------------------------------------------------
# Construct a composite unit from one numerator unit.
# --------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Construct a single-numerator [`Units`](@ref) value with an optional
denominator.

# Arguments

- `prefix`: Numerator metric prefix.
- `name`: Numerator unit name.

# Keywords

- `per`: Optional `(prefix, name)` denominator. Default: `nothing`.

# Returns

- A [`Units`](@ref) description.
"""
function units(
        prefix::Symbol,
        name::Symbol;
        per::Union{Nothing, Tuple{Symbol, Symbol}} = nothing
)
    b = Unit(name = name, prefix = prefix)
    if per === nothing
        return Units(base = [b], per = Unit[])
    else
        pfx2, name2 = per
        d = Unit(name = name2, prefix = pfx2)
        return Units(base = [b], per = [d])
    end
end

# --------------------------------------------------------------------------
# Formatting
# --------------------------------------------------------------------------

"""
$(TYPEDSIGNATURES)

Return the metric prefix and SI symbol for `u`.
"""
function get_label(u::Unit)
    prefix_str = _prefix_symbol(u.prefix)
    unit_str = get(UNIT_SYMBOL, u.name, String(u.name))
    return string(prefix_str, unit_str)
end

"""
$(TYPEDSIGNATURES)

Return a display label for composite units. A dot joins units within the
numerator or denominator. Parentheses enclose a denominator with several
units. Dimensionless entries are omitted.
"""
function get_label(u::Units)
    base_units = [x for x in u.base if x.name != :dimensionless]
    per_units = [x for x in u.per if x.name != :dimensionless]

    base_strs = [get_label(b) for b in base_units if get_label(b) != ""]
    per_strs = [get_label(p) for p in per_units if get_label(p) != ""]

    if isempty(base_strs) && isempty(per_strs)
        return ""
    end

    base_str = isempty(base_strs) ? "" : join(base_strs, ".")

    if isempty(per_strs)
        return base_str
    else
        if isempty(base_str)
            base_str = "1"
        end
        if length(per_strs) == 1
            return string(base_str, "/", per_strs[1])
        else
            return string(base_str, "/(", join(per_strs, "."), ")")
        end
    end
end
