# The quantities without a selector.
native_unit(::Quantity{Frequency, Nothing}) = units(:base, :hertz)
display_unit(::Quantity{Frequency, Nothing}) = units(:base, :hertz)
label(::Quantity{Frequency, Nothing}) = "Frequency"
symbol(::Quantity{Frequency, Nothing}) = "f"

native_unit(::Quantity{Dimensionless, Nothing}) = units(:base, :dimensionless)
display_unit(::Quantity{Dimensionless, Nothing}) = units(:base, :dimensionless)
label(::Quantity{Dimensionless, Nothing}) = "Dimensionless"
symbol(::Quantity{Dimensionless, Nothing}) = ""

# Derived-unit rules. A component reads its unit, label and symbol from its base
# quantity. A magnitude, real part or imaginary part keeps the base unit, and an angle
# is in radians, shown in degrees.
function native_unit(
        ::Quantity{Name, <:Union{typeof(abs), typeof(real), typeof(imag)}}) where {Name}
    return native_unit(Quantity{Name, Nothing}())
end
function display_unit(
        ::Quantity{Name, <:Union{typeof(abs), typeof(real), typeof(imag)}}) where {Name}
    return display_unit(Quantity{Name, Nothing}())
end
native_unit(::Quantity{Name, typeof(angle)}) where {Name} = units(:base, :radian)
display_unit(::Quantity{Name, typeof(angle)}) where {Name} = units(:base, :degree)

label(::Quantity{Name, typeof(abs)}) where {Name} =
    label(Quantity{Name, Nothing}()) * " magnitude"
label(::Quantity{Name, typeof(angle)}) where {Name} =
    label(Quantity{Name, Nothing}()) * " angle"
label(::Quantity{Name, typeof(real)}) where {Name} =
    "Real " * lowercase(label(Quantity{Name, Nothing}()))
label(::Quantity{Name, typeof(imag)}) where {Name} =
    "Imaginary " * lowercase(label(Quantity{Name, Nothing}()))
symbol(::Quantity{Name, typeof(abs)}) where {Name} =
    "|" * symbol(Quantity{Name, Nothing}()) * "|"
symbol(::Quantity{Name, typeof(angle)}) where {Name} =
    "∠" * symbol(Quantity{Name, Nothing}())
symbol(::Quantity{Name, typeof(real)}) where {Name} =
    "Re(" * symbol(Quantity{Name, Nothing}()) * ")"
symbol(::Quantity{Name, typeof(imag)}) where {Name} =
    "Im(" * symbol(Quantity{Name, Nothing}()) * ")"

family(::Quantity{Name, Nothing}) where {Name} = nothing
family(::Quantity{Name, <:Function}) where {Name} = family(Quantity{Name, Nothing}())

function label(quantity::Quantity, unit::UnitExpr)
    unit_text = label(unit)
    return isempty(unit_text) ? label(quantity) : "$(label(quantity)) [$unit_text]"
end
