native_unit(::Quantity{:frequency, Nothing}) = units(:base, :hertz)
display_unit(::Quantity{:frequency, Nothing}) = units(:base, :hertz)
label(::Quantity{:frequency, Nothing}) = "Frequency"
symbol(::Quantity{:frequency, Nothing}) = "f"

native_unit(::Quantity{:series_resistance, Nothing}) =
    units(:base, :ohm; per = (:base, :meter))
display_unit(::Quantity{:series_resistance, Nothing}) =
    units(:base, :ohm; per = (:kilo, :meter))
label(::Quantity{:series_resistance, Nothing}) = "Series resistance"
symbol(::Quantity{:series_resistance, Nothing}) = "R"

native_unit(::Quantity{:series_reactance, Nothing}) =
    units(:base, :ohm; per = (:base, :meter))
display_unit(::Quantity{:series_reactance, Nothing}) =
    units(:base, :ohm; per = (:kilo, :meter))
label(::Quantity{:series_reactance, Nothing}) = "Series reactance"
symbol(::Quantity{:series_reactance, Nothing}) = "X"

native_unit(::Quantity{:series_inductance, Nothing}) =
    units(:base, :henry; per = (:base, :meter))
display_unit(::Quantity{:series_inductance, Nothing}) =
    units(:milli, :henry; per = (:kilo, :meter))
label(::Quantity{:series_inductance, Nothing}) = "Series inductance"
symbol(::Quantity{:series_inductance, Nothing}) = "L"

native_unit(::Quantity{:series_impedance, Nothing}) =
    units(:base, :ohm; per = (:base, :meter))
display_unit(::Quantity{:series_impedance, Nothing}) =
    units(:base, :ohm; per = (:kilo, :meter))
label(::Quantity{:series_impedance, Nothing}) = "Series impedance"
symbol(::Quantity{:series_impedance, Nothing}) = "Z"

native_unit(::Quantity{:shunt_conductance, Nothing}) =
    units(:base, :siemens; per = (:base, :meter))
display_unit(::Quantity{:shunt_conductance, Nothing}) =
    units(:base, :siemens; per = (:kilo, :meter))
label(::Quantity{:shunt_conductance, Nothing}) = "Shunt conductance"
symbol(::Quantity{:shunt_conductance, Nothing}) = "G"

native_unit(::Quantity{:shunt_susceptance, Nothing}) =
    units(:base, :siemens; per = (:base, :meter))
display_unit(::Quantity{:shunt_susceptance, Nothing}) =
    units(:base, :siemens; per = (:kilo, :meter))
label(::Quantity{:shunt_susceptance, Nothing}) = "Shunt susceptance"
symbol(::Quantity{:shunt_susceptance, Nothing}) = "B"

native_unit(::Quantity{:shunt_capacitance, Nothing}) =
    units(:base, :farad; per = (:base, :meter))
display_unit(::Quantity{:shunt_capacitance, Nothing}) =
    units(:micro, :farad; per = (:kilo, :meter))
label(::Quantity{:shunt_capacitance, Nothing}) = "Shunt capacitance"
symbol(::Quantity{:shunt_capacitance, Nothing}) = "C"

native_unit(::Quantity{:shunt_admittance, Nothing}) =
    units(:base, :siemens; per = (:base, :meter))
display_unit(::Quantity{:shunt_admittance, Nothing}) =
    units(:base, :siemens; per = (:kilo, :meter))
label(::Quantity{:shunt_admittance, Nothing}) = "Shunt admittance"
symbol(::Quantity{:shunt_admittance, Nothing}) = "Y"

native_unit(::Quantity{:phase_angle, Nothing}) = units(:base, :radian)
display_unit(::Quantity{:phase_angle, Nothing}) = units(:base, :degree)
label(::Quantity{:phase_angle, Nothing}) = "Phase angle"
symbol(::Quantity{:phase_angle, Nothing}) = "∠"

native_unit(::Quantity{:dimensionless, Nothing}) = units(:base, :dimensionless)
display_unit(::Quantity{:dimensionless, Nothing}) = units(:base, :dimensionless)
label(::Quantity{:dimensionless, Nothing}) = "Dimensionless"
symbol(::Quantity{:dimensionless, Nothing}) = ""

# Derived-unit rules. A component reads its unit, label and symbol from its base
# quantity. A magnitude, real part or imaginary part keeps the base unit, and an angle
# takes the unit of the phase angle.
function native_unit(
        ::Quantity{Name, <:Union{typeof(abs), typeof(real), typeof(imag)}}) where {Name}
    return native_unit(Quantity{Name, Nothing}())
end
function display_unit(
        ::Quantity{Name, <:Union{typeof(abs), typeof(real), typeof(imag)}}) where {Name}
    return display_unit(Quantity{Name, Nothing}())
end
native_unit(::Quantity{Name, typeof(angle)}) where {Name} =
    native_unit(Quantity{:phase_angle, Nothing}())
display_unit(::Quantity{Name, typeof(angle)}) where {Name} =
    display_unit(Quantity{:phase_angle, Nothing}())

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

const _SeriesQuantity = Union{
    Quantity{:series_resistance, Nothing},
    Quantity{:series_reactance, Nothing},
    Quantity{:series_inductance, Nothing},
    Quantity{:series_impedance, Nothing}
}

const _ShuntQuantity = Union{
    Quantity{:shunt_conductance, Nothing},
    Quantity{:shunt_susceptance, Nothing},
    Quantity{:shunt_capacitance, Nothing},
    Quantity{:shunt_admittance, Nothing}
}

family(::_SeriesQuantity) = Val(:series)
family(::_ShuntQuantity) = Val(:shunt)
family(::Quantity{Name, Nothing}) where {Name} = nothing
family(::Quantity{Name, <:Function}) where {Name} = family(Quantity{Name, Nothing}())

function label(quantity::Quantity, unit::UnitExpr)
    unit_text = label(unit)
    return isempty(unit_text) ? label(quantity) : "$(label(quantity)) [$unit_text]"
end
