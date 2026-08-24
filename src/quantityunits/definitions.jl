# --------------------------------------------------------------------------
# Base quantities
# --------------------------------------------------------------------------

# Frequency
default_unit(::QuantityTag{:frequency}) = units(:base, :hertz)
display_unit(::QuantityTag{:frequency}) = units(:base, :hertz)
get_label(::QuantityTag{:frequency}) = "Frequency"
get_symbol(::QuantityTag{:frequency}) = "f"

# Series resistance
default_unit(::QuantityTag{:resistance}) = units(:base, :ohm; per = (:base, :meter))    # Ω/m
display_unit(::QuantityTag{:resistance}) = units(:base, :ohm; per = (:kilo, :meter))    # Ω/km
get_label(::QuantityTag{:resistance}) = "Series resistance"
get_symbol(::QuantityTag{:resistance}) = "R"

# Series inductance
default_unit(::QuantityTag{:inductance}) = units(:base, :henry; per = (:base, :meter))
display_unit(::QuantityTag{:inductance}) = units(:milli, :henry; per = (:kilo, :meter))
get_label(::QuantityTag{:inductance}) = "Series inductance"
get_symbol(::QuantityTag{:inductance}) = "L"

# Shunt capacitance
default_unit(::QuantityTag{:capacitance}) = units(:base, :farad; per = (:base, :meter))
display_unit(::QuantityTag{:capacitance}) = units(:micro, :farad; per = (:kilo, :meter))
get_label(::QuantityTag{:capacitance}) = "Shunt capacitance"
get_symbol(::QuantityTag{:capacitance}) = "C"

# Shunt conductance
default_unit(::QuantityTag{:conductance}) = units(:base, :siemens; per = (:base, :meter))
display_unit(::QuantityTag{:conductance}) = units(:base, :siemens; per = (:kilo, :meter))
get_label(::QuantityTag{:conductance}) = "Shunt conductance"
get_symbol(::QuantityTag{:conductance}) = "G"

# Series impedance
default_unit(::QuantityTag{:series_impedance}) = units(:base, :ohm; per = (:base, :meter))
display_unit(::QuantityTag{:series_impedance}) = units(:base, :ohm; per = (:kilo, :meter))
get_label(::QuantityTag{:series_impedance}) = "Series impedance"
get_symbol(::QuantityTag{:series_impedance}) = "Z"

# Shunt admittance
function default_unit(::QuantityTag{:shunt_admittance})
    units(:base, :siemens; per = (:base, :meter))
end
function display_unit(::QuantityTag{:shunt_admittance})
    units(:base, :siemens; per = (:kilo, :meter))
end
get_label(::QuantityTag{:shunt_admittance}) = "Shunt admittance"
get_symbol(::QuantityTag{:shunt_admittance}) = "Y"

# Inductive reactance
default_unit(::QuantityTag{:reactance}) = units(:base, :ohm; per = (:base, :meter))
display_unit(::QuantityTag{:reactance}) = units(:base, :ohm; per = (:kilo, :meter))
get_label(::QuantityTag{:reactance}) = "Inductive reactance"
get_symbol(::QuantityTag{:reactance}) = "X"

# Capacitive susceptance
default_unit(::QuantityTag{:susceptance}) = units(:base, :siemens; per = (:base, :meter))
display_unit(::QuantityTag{:susceptance}) = units(:base, :siemens; per = (:kilo, :meter))
get_label(::QuantityTag{:susceptance}) = "Capacitive susceptance"
get_symbol(::QuantityTag{:susceptance}) = "B"

# Angle
default_unit(::QuantityTag{:angle}) = units(:base, :degree)
display_unit(::QuantityTag{:angle}) = units(:base, :degree)
get_label(::QuantityTag{:angle}) = "Angle"
get_symbol(::QuantityTag{:angle}) = "∠"

# magnitude uses same unit as base quantity
function default_unit(::QuantityTag{(:series_impedance, :re)})
    default_unit(QuantityTag{:resistance}())
end
function default_unit(::QuantityTag{(:series_impedance, :im)})
    default_unit(QuantityTag{:reactance}())
end
function default_unit(::QuantityTag{(:series_impedance, :abs)})
    default_unit(QuantityTag{:series_impedance}())
end
function default_unit(::QuantityTag{(:shunt_admittance, :re)})
    default_unit(QuantityTag{:conductance}())
end
function default_unit(::QuantityTag{(:shunt_admittance, :im)})
    default_unit(QuantityTag{:susceptance}())
end
function default_unit(::QuantityTag{(:shunt_admittance, :abs)})
    default_unit(QuantityTag{:shunt_admittance}())
end

# angle unit
function default_unit(::QuantityTag{(:series_impedance, :angle)})
    default_unit(QuantityTag{:angle}())
end
function default_unit(::QuantityTag{(:shunt_admittance, :angle)})
    default_unit(QuantityTag{:angle}())
end

# labels
get_label(::QuantityTag{(:series_impedance, :abs)}) = "Series impedance magnitude"
get_label(::QuantityTag{(:series_impedance, :angle)}) = "Series impedance angle"
get_label(::QuantityTag{(:shunt_admittance, :abs)}) = "Shunt admittance magnitude"
get_label(::QuantityTag{(:shunt_admittance, :angle)}) = "Shunt admittance angle"

function default_unit(::QuantityTag{:series_impedance_absolute_error})
    default_unit(QuantityTag{:series_impedance}())
end
function display_unit(::QuantityTag{:series_impedance_absolute_error})
    display_unit(QuantityTag{:series_impedance}())
end
function get_label(::QuantityTag{:series_impedance_absolute_error})
    "Series impedance absolute error"
end
get_symbol(::QuantityTag{:series_impedance_absolute_error}) = "ΔZ"

function default_unit(::QuantityTag{:shunt_admittance_absolute_error})
    default_unit(QuantityTag{:shunt_admittance}())
end
function display_unit(::QuantityTag{:shunt_admittance_absolute_error})
    display_unit(QuantityTag{:shunt_admittance}())
end
function get_label(::QuantityTag{:shunt_admittance_absolute_error})
    "Shunt admittance absolute error"
end
get_symbol(::QuantityTag{:shunt_admittance_absolute_error}) = "ΔY"

default_unit(::QuantityTag{:series_impedance_relative_error}) = units(:base, :dimensionless)
display_unit(::QuantityTag{:series_impedance_relative_error}) = units(:base, :dimensionless)
function get_label(::QuantityTag{:series_impedance_relative_error})
    "Series impedance relative error"
end
get_symbol(::QuantityTag{:series_impedance_relative_error}) = "εZ"

default_unit(::QuantityTag{:shunt_admittance_relative_error}) = units(:base, :dimensionless)
display_unit(::QuantityTag{:shunt_admittance_relative_error}) = units(:base, :dimensionless)
function get_label(::QuantityTag{:shunt_admittance_relative_error})
    "Shunt admittance relative error"
end
get_symbol(::QuantityTag{:shunt_admittance_relative_error}) = "εY"
