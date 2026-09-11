@inline _dimension(name::Symbol) = name in (:radian, :degree) ? :angle : name
@inline _unit_scale(unit::Unit) = 10.0^_prefix_exponent(unit.prefix) *
                                  (unit.name === :degree ? π / 180 : 1.0)

function _dimensions(expression::UnitExpr)
    return (
        map(unit -> _dimension(unit.name), expression.numerator),
        map(unit -> _dimension(unit.name), expression.denominator)
    )
end

function _expression_scale(expression::UnitExpr)
    numerator = prod(_unit_scale, expression.numerator; init = 1.0)
    denominator = prod(_unit_scale, expression.denominator; init = 1.0)
    return numerator / denominator
end

"""
$(TYPEDSIGNATURES)

Return the numeric factor that converts values between compatible unit
expressions.

# Errors

- Throws `ArgumentError` when the expressions represent different dimensions.
"""
function scale_factor(from::UnitExpr, to::UnitExpr)
    _dimensions(from) == _dimensions(to) || throw(
        ArgumentError("cannot convert $(label(from)) to incompatible unit $(label(to))"),
    )
    return _expression_scale(from) / _expression_scale(to)
end

function _with_basis(expression::UnitExpr, basis::Symbol; length_prefix::Symbol = :base)
    basis in (:pul, :total) || throw(
        ArgumentError("basis must be :pul or :total"),
    )
    has_length = any(unit -> unit.name === :meter, expression.denominator)
    retained = filter(unit -> unit.name !== :meter, expression.denominator)
    denominator = basis === :pul && has_length ?
                  (retained..., Unit(:meter, length_prefix)) : retained
    return UnitExpr(expression.numerator, denominator)
end

function native_unit(quantity::Quantity, basis::Symbol)
    _with_basis(native_unit(quantity), basis)
end
function display_unit(
        quantity::Quantity,
        basis::Symbol;
        length_prefix::Symbol = :kilo,
        prefix::Union{Nothing, Symbol} = nothing
)
    expression = _with_basis(display_unit(quantity), basis; length_prefix)
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
        basis::Symbol,
        override;
        length_prefix::Symbol = :kilo
)
    override === nothing &&
        return display_unit(quantity, basis; length_prefix)
    override isa UnitExpr && return override
    override isa Symbol || throw(ArgumentError(
        "display-unit overrides must be a metric prefix, UnitExpr, or nothing",
    ))
    return display_unit(quantity, basis; length_prefix, prefix = override)
end

function scale_factor(quantity::Quantity, target::UnitExpr)
    scale_factor(native_unit(quantity), target)
end
function scale_factor(quantity::Quantity, basis::Symbol, target::UnitExpr)
    scale_factor(native_unit(quantity, basis), target)
end
