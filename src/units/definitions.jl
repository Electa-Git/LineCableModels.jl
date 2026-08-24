native_unit(::QuantityTag{:frequency}) = units(:base, :hertz)
display_unit(::QuantityTag{:frequency}) = units(:base, :hertz)
label(::QuantityTag{:frequency}) = "Frequency"
symbol(::QuantityTag{:frequency}) = "f"

native_unit(::QuantityTag{:series_resistance}) = units(:base, :ohm; per = (:base, :meter))
display_unit(::QuantityTag{:series_resistance}) = units(:base, :ohm; per = (:kilo, :meter))
label(::QuantityTag{:series_resistance}) = "Series resistance"
symbol(::QuantityTag{:series_resistance}) = "R"

native_unit(::QuantityTag{:series_reactance}) = units(:base, :ohm; per = (:base, :meter))
display_unit(::QuantityTag{:series_reactance}) = units(:base, :ohm; per = (:kilo, :meter))
label(::QuantityTag{:series_reactance}) = "Series reactance"
symbol(::QuantityTag{:series_reactance}) = "X"

native_unit(::QuantityTag{:series_inductance}) = units(:base, :henry; per = (:base, :meter))
display_unit(::QuantityTag{:series_inductance}) = units(:milli, :henry; per = (:kilo, :meter))
label(::QuantityTag{:series_inductance}) = "Series inductance"
symbol(::QuantityTag{:series_inductance}) = "L"

native_unit(::QuantityTag{:series_impedance}) = units(:base, :ohm; per = (:base, :meter))
display_unit(::QuantityTag{:series_impedance}) = units(:base, :ohm; per = (:kilo, :meter))
label(::QuantityTag{:series_impedance}) = "Series impedance"
symbol(::QuantityTag{:series_impedance}) = "Z"

native_unit(::QuantityTag{:shunt_conductance}) = units(:base, :siemens; per = (:base, :meter))
display_unit(::QuantityTag{:shunt_conductance}) = units(:base, :siemens; per = (:kilo, :meter))
label(::QuantityTag{:shunt_conductance}) = "Shunt conductance"
symbol(::QuantityTag{:shunt_conductance}) = "G"

native_unit(::QuantityTag{:shunt_susceptance}) = units(:base, :siemens; per = (:base, :meter))
display_unit(::QuantityTag{:shunt_susceptance}) = units(:base, :siemens; per = (:kilo, :meter))
label(::QuantityTag{:shunt_susceptance}) = "Shunt susceptance"
symbol(::QuantityTag{:shunt_susceptance}) = "B"

native_unit(::QuantityTag{:shunt_capacitance}) = units(:base, :farad; per = (:base, :meter))
display_unit(::QuantityTag{:shunt_capacitance}) = units(:micro, :farad; per = (:kilo, :meter))
label(::QuantityTag{:shunt_capacitance}) = "Shunt capacitance"
symbol(::QuantityTag{:shunt_capacitance}) = "C"

native_unit(::QuantityTag{:shunt_admittance}) = units(:base, :siemens; per = (:base, :meter))
display_unit(::QuantityTag{:shunt_admittance}) = units(:base, :siemens; per = (:kilo, :meter))
label(::QuantityTag{:shunt_admittance}) = "Shunt admittance"
symbol(::QuantityTag{:shunt_admittance}) = "Y"

native_unit(::QuantityTag{:phase_angle}) = units(:base, :radian)
display_unit(::QuantityTag{:phase_angle}) = units(:base, :degree)
label(::QuantityTag{:phase_angle}) = "Phase angle"
symbol(::QuantityTag{:phase_angle}) = "∠"

native_unit(::QuantityTag{:distance}) = units(:base, :meter)
display_unit(::QuantityTag{:distance}) = units(:base, :meter)
label(::QuantityTag{:distance}) = "Distance"
symbol(::QuantityTag{:distance}) = "d"

native_unit(::QuantityTag{:dimensionless}) = units(:base, :dimensionless)
display_unit(::QuantityTag{:dimensionless}) = units(:base, :dimensionless)
label(::QuantityTag{:dimensionless}) = "Dimensionless"
symbol(::QuantityTag{:dimensionless}) = ""

native_unit(::QuantityTag{(:series_impedance, :magnitude)}) =
    native_unit(QuantityTag{:series_impedance}())
display_unit(::QuantityTag{(:series_impedance, :magnitude)}) =
    display_unit(QuantityTag{:series_impedance}())
native_unit(::QuantityTag{(:shunt_admittance, :magnitude)}) =
    native_unit(QuantityTag{:shunt_admittance}())
display_unit(::QuantityTag{(:shunt_admittance, :magnitude)}) =
    display_unit(QuantityTag{:shunt_admittance}())
native_unit(::QuantityTag{(:series_impedance, :phase_angle)}) =
    native_unit(QuantityTag{:phase_angle}())
display_unit(::QuantityTag{(:series_impedance, :phase_angle)}) =
    display_unit(QuantityTag{:phase_angle}())
native_unit(::QuantityTag{(:shunt_admittance, :phase_angle)}) =
    native_unit(QuantityTag{:phase_angle}())
display_unit(::QuantityTag{(:shunt_admittance, :phase_angle)}) =
    display_unit(QuantityTag{:phase_angle}())

label(::QuantityTag{(:series_impedance, :magnitude)}) = "Series impedance magnitude"
symbol(::QuantityTag{(:series_impedance, :magnitude)}) = "|Z|"
label(::QuantityTag{(:series_impedance, :phase_angle)}) = "Series impedance phase angle"
symbol(::QuantityTag{(:series_impedance, :phase_angle)}) = "∠Z"
label(::QuantityTag{(:shunt_admittance, :magnitude)}) = "Shunt admittance magnitude"
symbol(::QuantityTag{(:shunt_admittance, :magnitude)}) = "|Y|"
label(::QuantityTag{(:shunt_admittance, :phase_angle)}) = "Shunt admittance phase angle"
symbol(::QuantityTag{(:shunt_admittance, :phase_angle)}) = "∠Y"

native_unit(::QuantityTag{:series_impedance_absolute_error}) =
    native_unit(QuantityTag{:series_impedance}())
display_unit(::QuantityTag{:series_impedance_absolute_error}) =
    display_unit(QuantityTag{:series_impedance}())
label(::QuantityTag{:series_impedance_absolute_error}) = "Series impedance absolute error"
symbol(::QuantityTag{:series_impedance_absolute_error}) = "ΔZ"

native_unit(::QuantityTag{:shunt_admittance_absolute_error}) =
    native_unit(QuantityTag{:shunt_admittance}())
display_unit(::QuantityTag{:shunt_admittance_absolute_error}) =
    display_unit(QuantityTag{:shunt_admittance}())
label(::QuantityTag{:shunt_admittance_absolute_error}) = "Shunt admittance absolute error"
symbol(::QuantityTag{:shunt_admittance_absolute_error}) = "ΔY"

native_unit(::QuantityTag{:series_impedance_relative_error}) = units(:base, :dimensionless)
display_unit(::QuantityTag{:series_impedance_relative_error}) = units(:base, :dimensionless)
label(::QuantityTag{:series_impedance_relative_error}) = "Series impedance relative error"
symbol(::QuantityTag{:series_impedance_relative_error}) = "εZ"

native_unit(::QuantityTag{:shunt_admittance_relative_error}) = units(:base, :dimensionless)
display_unit(::QuantityTag{:shunt_admittance_relative_error}) = units(:base, :dimensionless)
label(::QuantityTag{:shunt_admittance_relative_error}) = "Shunt admittance relative error"
symbol(::QuantityTag{:shunt_admittance_relative_error}) = "εY"

function label(quantity::QuantityTag, unit::UnitExpr)
    unit_text = label(unit)
    return isempty(unit_text) ? label(quantity) : "$(label(quantity)) [$unit_text]"
end
