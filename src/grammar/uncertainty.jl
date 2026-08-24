"Return the nominal value of a deterministic or uncertain quantity."
nominal(value) = value
nominal(value::Complex) = complex(nominal(real(value)), nominal(imag(value)))
nominal(values::AbstractArray) = nominal.(values)

"Return the standard uncertainty of a quantity. Deterministic numbers return zero."
standard_uncertainty(value::Number) = zero(nominal(value))
standard_uncertainty(::Any) = 0.0
