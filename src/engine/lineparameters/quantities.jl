# Engine-owned line-parameter quantity and unit accessors.
Units.quantity(::typeof(frequencies)) = Units.QuantityTag{:frequency}()
Units.quantity(::typeof(Z)) = Units.QuantityTag{:series_impedance}()
Units.quantity(::typeof(Y)) = Units.QuantityTag{:shunt_admittance}()
Units.quantity(::typeof(X)) = Units.QuantityTag{:series_reactance}()
Units.quantity(::typeof(G)) = Units.QuantityTag{:shunt_conductance}()
Units.quantity(::typeof(B)) = Units.QuantityTag{:shunt_susceptance}()

Units.quantity(::typeof(Z), ::typeof(abs)) =
    Units.QuantityTag{(:series_impedance, :magnitude)}()
Units.quantity(::typeof(Z), ::typeof(angle)) =
    Units.QuantityTag{(:series_impedance, :phase_angle)}()
Units.quantity(::typeof(Y), ::typeof(abs)) =
    Units.QuantityTag{(:shunt_admittance, :magnitude)}()
Units.quantity(::typeof(Y), ::typeof(angle)) =
    Units.QuantityTag{(:shunt_admittance, :phase_angle)}()

Units.quantity(::typeof(Z), ::typeof(absolute_error)) =
    Units.QuantityTag{:series_impedance_absolute_error}()
Units.quantity(::typeof(Z), ::typeof(relative_error)) =
    Units.QuantityTag{:series_impedance_relative_error}()
Units.quantity(::typeof(Y), ::typeof(absolute_error)) =
    Units.QuantityTag{:shunt_admittance_absolute_error}()
Units.quantity(::typeof(Y), ::typeof(relative_error)) =
    Units.QuantityTag{:shunt_admittance_relative_error}()
