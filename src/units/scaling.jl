@inline _dimension(name::Symbol) = name in (:radian, :degree) ? :angle : name
function _dimensions(expression::UnitExpr)
    return (
        map(unit -> _dimension(unit.name), expression.numerator),
        map(unit -> _dimension(unit.name), expression.denominator)
    )
end

"""
$(TYPEDSIGNATURES)

Return the numeric factor that converts values between compatible unit
expressions.

# Errors

- Throws `ArgumentError` when the expressions represent different dimensions.
"""
function scale_factor(from::UnitExpr, to::UnitExpr)
    return scale_factor(from,to,Float64)
end

"""
$(TYPEDSIGNATURES)

Convert compatible units using constants evaluated in scalar precision `T`.
Conversion between degrees and radians preserves this precision instead of
rounding π through `Float64`.
"""
function scale_factor(from::UnitExpr,to::UnitExpr,::Type{T}) where {T<:AbstractFloat}
    _dimensions(from) == _dimensions(to) || throw(
        ArgumentError("cannot convert $(label(from)) to incompatible unit $(label(to))"),
    )
    exponent(expression)=sum(unit -> _prefix_exponent(unit.prefix),expression.numerator;init=0) -
        sum(unit -> _prefix_exponent(unit.prefix),expression.denominator;init=0)
    degrees(expression)=count(unit -> unit.name===:degree,expression.numerator) -
        count(unit -> unit.name===:degree,expression.denominator)
    return T(10)^(exponent(from)-exponent(to)) *
        (T(π)/T(180))^(degrees(from)-degrees(to))
end

"""
$(TYPEDSIGNATURES)

Resolve a display unit for a quantity from one optional override.

The value `nothing` retains the registered display unit, with `length_prefix` as the
metric prefix of a per-length denominator. A `Symbol` also replaces its numerator metric
prefix, and a `UnitExpr` is returned unchanged.

# Errors

- Throws `ArgumentError` for unsupported override values or invalid metric
  prefixes.
"""
function display_unit(quantity::Quantity, override; length_prefix::Symbol = :kilo)
    throw(ArgumentError(
        "display-unit overrides must be a metric prefix, UnitExpr, or nothing",
    ))
end
function display_unit(quantity::Quantity, ::Nothing; length_prefix::Symbol = :kilo)
    expression = display_unit(quantity)
    has_length = any(unit -> unit.name === :meter, expression.denominator)
    retained = filter(unit -> unit.name !== :meter, expression.denominator)
    denominator = has_length ? (retained..., Unit(:meter, length_prefix)) : retained
    return UnitExpr(expression.numerator, denominator)
end
function display_unit(quantity::Quantity, override::UnitExpr; length_prefix::Symbol = :kilo)
    return override
end
function display_unit(quantity::Quantity, override::Symbol; length_prefix::Symbol = :kilo)
    expression = display_unit(quantity, nothing; length_prefix)
    isempty(expression.numerator) && return expression
    numerator = (Unit(first(expression.numerator).name, override),
        Base.tail(expression.numerator)...)
    return UnitExpr(numerator, expression.denominator)
end

function scale_factor(quantity::Quantity, target::UnitExpr)
    scale_factor(native_unit(quantity), target)
end
