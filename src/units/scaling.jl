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

# Per unit length keeps the length in the denominator, with `length_prefix`. The total
# line length removes it.
function _with_basis(expression::UnitExpr, ::PerUnitLength; length_prefix::Symbol = :base)
    has_length = any(unit -> unit.name === :meter, expression.denominator)
    retained = filter(unit -> unit.name !== :meter, expression.denominator)
    denominator = has_length ? (retained..., Unit(:meter, length_prefix)) : retained
    return UnitExpr(expression.numerator, denominator)
end
function _with_basis(expression::UnitExpr, ::TotalLength; length_prefix::Symbol = :base)
    retained = filter(unit -> unit.name !== :meter, expression.denominator)
    return UnitExpr(expression.numerator, retained)
end

function native_unit(quantity::Quantity, basis::Union{Symbol, LengthBasis})
    _with_basis(native_unit(quantity), LengthBasis(basis))
end
function display_unit(
        quantity::Quantity,
        basis::Union{Symbol, LengthBasis};
        length_prefix::Symbol = :kilo,
        prefix::Union{Nothing, Symbol} = nothing
)
    expression = _with_basis(display_unit(quantity), LengthBasis(basis); length_prefix)
    prefix === nothing && return expression
    isempty(expression.numerator) && return expression
    numerator = (
        Unit(first(expression.numerator).name, prefix),
        Base.tail(expression.numerator)...
    )
    return UnitExpr(numerator, expression.denominator)
end

"""
$(TYPEDSIGNATURES)

Resolve a display unit for a quantity and basis from one optional override.

The value `nothing` retains the registered display unit, a `Symbol` replaces
its numerator metric prefix, and a `UnitExpr` is returned unchanged.

# Errors

- Throws `ArgumentError` for unsupported override values or invalid metric
  prefixes.
"""
function display_unit(
        quantity::Quantity,
        basis::Union{Symbol, LengthBasis},
        override;
        length_prefix::Symbol = :kilo
)
    throw(ArgumentError(
        "display-unit overrides must be a metric prefix, UnitExpr, or nothing",
    ))
end
function display_unit(quantity::Quantity, basis::Union{Symbol, LengthBasis}, ::Nothing;
        length_prefix::Symbol = :kilo)
    return display_unit(quantity, basis; length_prefix)
end
function display_unit(quantity::Quantity, basis::Union{Symbol, LengthBasis}, override::UnitExpr;
        length_prefix::Symbol = :kilo)
    return override
end
function display_unit(quantity::Quantity, basis::Union{Symbol, LengthBasis}, override::Symbol;
        length_prefix::Symbol = :kilo)
    return display_unit(quantity, basis; length_prefix, prefix = override)
end

function scale_factor(quantity::Quantity, target::UnitExpr)
    scale_factor(native_unit(quantity), target)
end
function scale_factor(quantity::Quantity, basis::Union{Symbol, LengthBasis}, target::UnitExpr)
    scale_factor(native_unit(quantity, basis), target)
end
