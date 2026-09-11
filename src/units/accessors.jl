# Shared physical selectors extended independently by DataModel and Engine.
quantity(::typeof(R)) = Quantity{:series_resistance}()
quantity(::typeof(L)) = Quantity{:series_inductance}()
quantity(::typeof(C)) = Quantity{:shunt_capacitance}()

"Return the scientific label registered for `selector`."
label(selector::Function) = label(quantity(selector))

"Return the scientific label registered for `selector` and `transform`."
label(selector::Function, transform::Function) = label(quantity(selector, transform))

"Return the scientific symbol registered for `selector`."
symbol(selector::Function) = symbol(quantity(selector))

"Return the scientific symbol registered for `selector` and `transform`."
symbol(selector::Function, transform::Function) = symbol(quantity(selector, transform))

"Return the native unit registered for `selector`."
native_unit(selector::Function) = native_unit(quantity(selector))

"Return the native unit registered for `selector` at `basis`."
native_unit(selector::Function, basis::Symbol) = native_unit(quantity(selector), basis)

"Return the native unit registered for `selector` and `transform`."
function native_unit(selector::Function, transform::Function)
    native_unit(quantity(selector, transform))
end

"Return the native unit registered for `selector` and `transform` at `basis`."
function native_unit(selector::Function, transform::Function, basis::Symbol)
    native_unit(quantity(selector, transform), basis)
end

"Return the display unit registered for `selector`."
display_unit(selector::Function) = display_unit(quantity(selector))

"Return the display unit registered for `selector` at `basis`."
function display_unit(selector::Function, basis::Symbol; kwargs...)
    return display_unit(quantity(selector), basis; kwargs...)
end

"Return the display unit registered for `selector` and `transform`."
function display_unit(selector::Function, transform::Function)
    display_unit(quantity(selector, transform))
end

"Return the display unit registered for `selector` and `transform` at `basis`."
function display_unit(
        selector::Function,
        transform::Function,
        basis::Symbol;
        kwargs...
)
    return display_unit(quantity(selector, transform), basis; kwargs...)
end

"Return the factor from the native unit of `selector` to `target`."
function scale_factor(selector::Function, target::UnitExpr)
    scale_factor(quantity(selector), target)
end

"Return the factor from the native unit of `selector` at `basis` to `target`."
function scale_factor(selector::Function, basis::Symbol, target::UnitExpr)
    scale_factor(quantity(selector), basis, target)
end

"Return the factor from the native unit of `selector` and `transform` to `target`."
function scale_factor(selector::Function, transform::Function, target::UnitExpr)
    scale_factor(quantity(selector, transform), target)
end

"Return the factor from the native unit of `selector` and `transform` at `basis` to `target`."
function scale_factor(
        selector::Function,
        transform::Function,
        basis::Symbol,
        target::UnitExpr
)
    return scale_factor(quantity(selector, transform), basis, target)
end
