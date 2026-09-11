const _METRIC_PREFIX_EXPONENT = (
    yocto = -24, zepto = -21, atto = -18, femto = -15, pico = -12,
    nano = -9, micro = -6, milli = -3, centi = -2, deci = -1,
    base = 0, deca = 1, hecto = 2, kilo = 3, mega = 6, giga = 9,
    tera = 12, peta = 15, exa = 18, zetta = 21, yotta = 24
)

const _METRIC_PREFIX_SYMBOL = (
    yocto = "y", zepto = "z", atto = "a", femto = "f", pico = "p",
    nano = "n", micro = "μ", milli = "m", centi = "c", deci = "d",
    base = "", deca = "da", hecto = "h", kilo = "k", mega = "M",
    giga = "G", tera = "T", peta = "P", exa = "E", zetta = "Z", yotta = "Y"
)

const _UNIT_SYMBOL = (
    dimensionless = "", ohm = "Ω", henry = "H", farad = "F",
    siemens = "S", meter = "m", hertz = "Hz", second = "s",
    radian = "rad", degree = "°"
)

@inline function _prefix_exponent(prefix::Symbol)
    hasproperty(_METRIC_PREFIX_EXPONENT, prefix) || throw(
        ArgumentError("unsupported metric prefix :$prefix"),
    )
    return getproperty(_METRIC_PREFIX_EXPONENT, prefix)
end

@inline function _prefix_symbol(prefix::Symbol)
    hasproperty(_METRIC_PREFIX_SYMBOL, prefix) || throw(
        ArgumentError("unsupported metric prefix :$prefix"),
    )
    return getproperty(_METRIC_PREFIX_SYMBOL, prefix)
end

@inline function _unit_symbol(name::Symbol)
    hasproperty(_UNIT_SYMBOL, name) || throw(ArgumentError("unsupported unit :$name"))
    return getproperty(_UNIT_SYMBOL, name)
end

"""
$(TYPEDEF)

Describe one supported physical unit and its metric prefix.

$(TYPEDFIELDS)
"""
struct Unit
    "Unit name."
    name::Symbol
    "Metric prefix name."
    prefix::Symbol

    function Unit(name::Symbol, prefix::Symbol)
        _unit_symbol(name)
        _prefix_exponent(prefix)
        name in (:dimensionless, :radian, :degree) && prefix !== :base && throw(
            ArgumentError("unit :$name does not accept metric prefix :$prefix"),
        )
        return new(name, prefix)
    end
end

Unit(name::Symbol) = Unit(name, :base)
Unit(; name::Symbol = :dimensionless, prefix::Symbol = :base) = Unit(name, prefix)

"""
$(TYPEDEF)

Describe numerator and denominator units, such as Ω/km.

$(TYPEDFIELDS)
"""
struct UnitExpr{N <: Tuple, D <: Tuple}
    "Numerator units."
    numerator::N
    "Denominator units."
    denominator::D

    function UnitExpr(numerator::N, denominator::D) where {N <: Tuple, D <: Tuple}
        all(unit -> unit isa Unit, numerator) || throw(
            ArgumentError("unit-expression numerators must contain Unit values"),
        )
        all(unit -> unit isa Unit, denominator) || throw(
            ArgumentError("unit-expression denominators must contain Unit values"),
        )
        return new{N, D}(numerator, denominator)
    end
end

function UnitExpr(; numerator::Tuple = (Unit(),), denominator::Tuple = ())
    UnitExpr(numerator, denominator)
end

"""
$(TYPEDSIGNATURES)

Construct a unit expression with one numerator and an optional denominator.

# Arguments

- `prefix`: Numerator metric prefix.
- `name`: Numerator unit name.

# Keywords

- `per`: Optional `(prefix, name)` denominator.

# Returns

- A validated [`UnitExpr`](@ref).
"""
function units(
        prefix::Symbol,
        name::Symbol;
        per::Union{Nothing, Tuple{Symbol, Symbol}} = nothing
)
    numerator = (Unit(name, prefix),)
    denominator = per === nothing ? () : (Unit(last(per), first(per)),)
    return UnitExpr(numerator, denominator)
end

"Return the SI symbol, including its metric prefix, for `unit`."
label(unit::Unit) = string(_prefix_symbol(unit.prefix), _unit_symbol(unit.name))

function _unit_side_label(side::Tuple)
    labels = map(label, filter(unit -> unit.name !== :dimensionless, side))
    return join(labels, ".")
end

"""
$(TYPEDSIGNATURES)

Format a unit expression using SI symbols.
"""
function label(expression::UnitExpr)
    numerator = _unit_side_label(expression.numerator)
    denominator = _unit_side_label(expression.denominator)
    isempty(denominator) && return numerator
    isempty(numerator) && (numerator = "1")
    denominator = length(filter(unit -> unit.name !== :dimensionless,
        expression.denominator)) > 1 ? "($denominator)" : denominator
    return "$numerator/$denominator"
end

Base.inv(expression::UnitExpr) = UnitExpr(expression.denominator, expression.numerator)
