# A magnitude, angle, real or imaginary part is a component of the selector's quantity.
function quantity(selector::Function,
        transform::Union{typeof(abs), typeof(angle), typeof(real), typeof(imag)})
    return quantity(quantity(selector), transform)
end

# A diagonal keeps the selector's quantity. `observables` declares where it applies.
quantity(selector::Function, ::typeof(diag)) = quantity(selector)

"""
Return the scientific label registered for `selector`.
"""
label(selector::Function) = label(quantity(selector))

"""
Return the scientific label registered for `selector` and `transform`.
"""
label(selector::Function, transform::Function) = label(quantity(selector, transform))

"""
Return the scientific symbol registered for `selector`.
"""
symbol(selector::Function) = symbol(quantity(selector))

"""
Return the scientific symbol registered for `selector` and `transform`.
"""
symbol(selector::Function, transform::Function) = symbol(quantity(selector, transform))

"""
Return the native unit registered for `selector`.
"""
native_unit(selector::Function) = native_unit(quantity(selector))

"""
Return the native unit registered for `selector` and `transform`.
"""
function native_unit(selector::Function, transform::Function)
    native_unit(quantity(selector, transform))
end

"""
Return the display unit registered for `selector`.
"""
display_unit(selector::Function) = display_unit(quantity(selector))

"""
Return the display unit registered for `selector` and `transform`.
"""
function display_unit(selector::Function, transform::Function)
    display_unit(quantity(selector, transform))
end

"""
Return the factor from the native unit of `selector` to `target`.
"""
function scale_factor(selector::Function, target::UnitExpr)
    scale_factor(quantity(selector), target)
end

"""
Return the factor from the native unit of `selector` and `transform` to `target`.
"""
function scale_factor(selector::Function, transform::Function, target::UnitExpr)
    scale_factor(quantity(selector, transform), target)
end
