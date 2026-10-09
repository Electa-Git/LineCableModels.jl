"""
Return the nominal value of a deterministic or uncertain quantity.
"""
nominal(value) = value
nominal(value::Complex) = complex(nominal(real(value)), nominal(imag(value)))
nominal(values::AbstractArray) = nominal.(values)

"""
Return the standard uncertainty of a quantity. A deterministic number returns the zero
of its own type.
"""
uncertainty(value::Number) = zero(value)
uncertainty(value::Complex) = complex(uncertainty(real(value)), uncertainty(imag(value)))
uncertainty(values::AbstractArray) = uncertainty.(values)
