# Engine-owned line-parameter quantity and unit accessors.
Units.quantity(::typeof(frequencies)) = Units.Quantity{:frequency}()
Units.quantity(::typeof(Z)) = Units.Quantity{:series_impedance}()
Units.quantity(::typeof(Y)) = Units.Quantity{:shunt_admittance}()
Units.quantity(::typeof(X)) = Units.Quantity{:series_reactance}()
Units.quantity(::typeof(G)) = Units.Quantity{:shunt_conductance}()
Units.quantity(::typeof(B)) = Units.Quantity{:shunt_susceptance}()

Units.quantity(::typeof(Z), ::typeof(abs)) =
    Units.Quantity{(:series_impedance, :magnitude)}()
Units.quantity(::typeof(Z), ::typeof(angle)) =
    Units.Quantity{(:series_impedance, :phase_angle)}()
Units.quantity(::typeof(Y), ::typeof(abs)) =
    Units.Quantity{(:shunt_admittance, :magnitude)}()
Units.quantity(::typeof(Y), ::typeof(angle)) =
    Units.Quantity{(:shunt_admittance, :phase_angle)}()

Units.quantity(::typeof(Z), ::typeof(absolute_error)) =
    Units.Quantity{:series_impedance_absolute_error}()
Units.quantity(::typeof(Z), ::typeof(relative_error)) =
    Units.Quantity{:series_impedance_relative_error}()
Units.quantity(::typeof(Y), ::typeof(absolute_error)) =
    Units.Quantity{:shunt_admittance_absolute_error}()
Units.quantity(::typeof(Y), ::typeof(relative_error)) =
    Units.Quantity{:shunt_admittance_relative_error}()
