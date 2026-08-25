native_unit(::Quantity{:frequency}) = units(:base, :hertz)
display_unit(::Quantity{:frequency}) = units(:base, :hertz)
label(::Quantity{:frequency}) = "Frequency"
symbol(::Quantity{:frequency}) = "f"

native_unit(::Quantity{:series_resistance}) = units(:base, :ohm; per = (:base, :meter))
display_unit(::Quantity{:series_resistance}) = units(:base, :ohm; per = (:kilo, :meter))
label(::Quantity{:series_resistance}) = "Series resistance"
symbol(::Quantity{:series_resistance}) = "R"

native_unit(::Quantity{:series_reactance}) = units(:base, :ohm; per = (:base, :meter))
display_unit(::Quantity{:series_reactance}) = units(:base, :ohm; per = (:kilo, :meter))
label(::Quantity{:series_reactance}) = "Series reactance"
symbol(::Quantity{:series_reactance}) = "X"

native_unit(::Quantity{:series_inductance}) = units(:base, :henry; per = (:base, :meter))
display_unit(::Quantity{:series_inductance}) = units(:milli, :henry; per = (:kilo, :meter))
label(::Quantity{:series_inductance}) = "Series inductance"
symbol(::Quantity{:series_inductance}) = "L"

native_unit(::Quantity{:series_impedance}) = units(:base, :ohm; per = (:base, :meter))
display_unit(::Quantity{:series_impedance}) = units(:base, :ohm; per = (:kilo, :meter))
label(::Quantity{:series_impedance}) = "Series impedance"
symbol(::Quantity{:series_impedance}) = "Z"

native_unit(::Quantity{:shunt_conductance}) = units(:base, :siemens; per = (:base, :meter))
display_unit(::Quantity{:shunt_conductance}) = units(:base, :siemens; per = (:kilo, :meter))
label(::Quantity{:shunt_conductance}) = "Shunt conductance"
symbol(::Quantity{:shunt_conductance}) = "G"

native_unit(::Quantity{:shunt_susceptance}) = units(:base, :siemens; per = (:base, :meter))
display_unit(::Quantity{:shunt_susceptance}) = units(:base, :siemens; per = (:kilo, :meter))
label(::Quantity{:shunt_susceptance}) = "Shunt susceptance"
symbol(::Quantity{:shunt_susceptance}) = "B"

native_unit(::Quantity{:shunt_capacitance}) = units(:base, :farad; per = (:base, :meter))
display_unit(::Quantity{:shunt_capacitance}) = units(:micro, :farad; per = (:kilo, :meter))
label(::Quantity{:shunt_capacitance}) = "Shunt capacitance"
symbol(::Quantity{:shunt_capacitance}) = "C"

native_unit(::Quantity{:shunt_admittance}) = units(:base, :siemens; per = (:base, :meter))
display_unit(::Quantity{:shunt_admittance}) = units(:base, :siemens; per = (:kilo, :meter))
label(::Quantity{:shunt_admittance}) = "Shunt admittance"
symbol(::Quantity{:shunt_admittance}) = "Y"

native_unit(::Quantity{:phase_angle}) = units(:base, :radian)
display_unit(::Quantity{:phase_angle}) = units(:base, :degree)
label(::Quantity{:phase_angle}) = "Phase angle"
symbol(::Quantity{:phase_angle}) = "∠"

native_unit(::Quantity{:distance}) = units(:base, :meter)
display_unit(::Quantity{:distance}) = units(:base, :meter)
label(::Quantity{:distance}) = "Distance"
symbol(::Quantity{:distance}) = "d"

native_unit(::Quantity{:dimensionless}) = units(:base, :dimensionless)
display_unit(::Quantity{:dimensionless}) = units(:base, :dimensionless)
label(::Quantity{:dimensionless}) = "Dimensionless"
symbol(::Quantity{:dimensionless}) = ""

native_unit(::Quantity{(:series_impedance, :magnitude)}) =
    native_unit(Quantity{:series_impedance}())
display_unit(::Quantity{(:series_impedance, :magnitude)}) =
    display_unit(Quantity{:series_impedance}())
native_unit(::Quantity{(:shunt_admittance, :magnitude)}) =
    native_unit(Quantity{:shunt_admittance}())
display_unit(::Quantity{(:shunt_admittance, :magnitude)}) =
    display_unit(Quantity{:shunt_admittance}())
native_unit(::Quantity{(:series_impedance, :phase_angle)}) =
    native_unit(Quantity{:phase_angle}())
display_unit(::Quantity{(:series_impedance, :phase_angle)}) =
    display_unit(Quantity{:phase_angle}())
native_unit(::Quantity{(:shunt_admittance, :phase_angle)}) =
    native_unit(Quantity{:phase_angle}())
display_unit(::Quantity{(:shunt_admittance, :phase_angle)}) =
    display_unit(Quantity{:phase_angle}())

label(::Quantity{(:series_impedance, :magnitude)}) = "Series impedance magnitude"
symbol(::Quantity{(:series_impedance, :magnitude)}) = "|Z|"
label(::Quantity{(:series_impedance, :phase_angle)}) = "Series impedance phase angle"
symbol(::Quantity{(:series_impedance, :phase_angle)}) = "∠Z"
label(::Quantity{(:shunt_admittance, :magnitude)}) = "Shunt admittance magnitude"
symbol(::Quantity{(:shunt_admittance, :magnitude)}) = "|Y|"
label(::Quantity{(:shunt_admittance, :phase_angle)}) = "Shunt admittance phase angle"
symbol(::Quantity{(:shunt_admittance, :phase_angle)}) = "∠Y"

native_unit(::Quantity{:series_impedance_absolute_error}) =
    native_unit(Quantity{:series_impedance}())
display_unit(::Quantity{:series_impedance_absolute_error}) =
    display_unit(Quantity{:series_impedance}())
label(::Quantity{:series_impedance_absolute_error}) = "Series impedance absolute error"
symbol(::Quantity{:series_impedance_absolute_error}) = "ΔZ"

native_unit(::Quantity{:shunt_admittance_absolute_error}) =
    native_unit(Quantity{:shunt_admittance}())
display_unit(::Quantity{:shunt_admittance_absolute_error}) =
    display_unit(Quantity{:shunt_admittance}())
label(::Quantity{:shunt_admittance_absolute_error}) = "Shunt admittance absolute error"
symbol(::Quantity{:shunt_admittance_absolute_error}) = "ΔY"

native_unit(::Quantity{:series_impedance_relative_error}) = units(:base, :dimensionless)
display_unit(::Quantity{:series_impedance_relative_error}) = units(:base, :dimensionless)
label(::Quantity{:series_impedance_relative_error}) = "Series impedance relative error"
symbol(::Quantity{:series_impedance_relative_error}) = "εZ"

native_unit(::Quantity{:shunt_admittance_relative_error}) = units(:base, :dimensionless)
display_unit(::Quantity{:shunt_admittance_relative_error}) = units(:base, :dimensionless)
label(::Quantity{:shunt_admittance_relative_error}) = "Shunt admittance relative error"
symbol(::Quantity{:shunt_admittance_relative_error}) = "εY"

function label(quantity::Quantity, unit::UnitExpr)
    unit_text = label(unit)
    return isempty(unit_text) ? label(quantity) : "$(label(quantity)) [$unit_text]"
end
