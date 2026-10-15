# Engine-owned line-parameter quantity and unit accessors.
Units.quantity(::typeof(frequencies)) = Units.Quantity{:frequency, Nothing}()
Units.quantity(::typeof(Z)) = Units.Quantity{:series_impedance, Nothing}()
Units.quantity(::typeof(Y)) = Units.Quantity{:shunt_admittance, Nothing}()
Units.quantity(::typeof(X)) = Units.Quantity{:series_reactance, Nothing}()
Units.quantity(::typeof(G)) = Units.Quantity{:shunt_conductance, Nothing}()
Units.quantity(::typeof(B)) = Units.Quantity{:shunt_susceptance, Nothing}()

# Error rules. An error is a component of the selector's quantity. An absolute error
# keeps the base unit, and a relative error is dimensionless.
function Units.quantity(selector::Function,
        transform::Union{typeof(absolute_error), typeof(relative_error)})
    return Units.quantity(Units.quantity(selector), transform)
end
Units.native_unit(::Units.Quantity{Name, typeof(absolute_error)}) where {Name} =
    Units.native_unit(Units.Quantity{Name, Nothing}())
Units.display_unit(::Units.Quantity{Name, typeof(absolute_error)}) where {Name} =
    Units.display_unit(Units.Quantity{Name, Nothing}())
Units.native_unit(::Units.Quantity{Name, typeof(relative_error)}) where {Name} =
    Units.units(:base, :dimensionless)
Units.display_unit(::Units.Quantity{Name, typeof(relative_error)}) where {Name} =
    Units.units(:base, :dimensionless)
Units.label(::Units.Quantity{Name, typeof(absolute_error)}) where {Name} =
    Units.label(Units.Quantity{Name, Nothing}()) * " absolute error"
Units.label(::Units.Quantity{Name, typeof(relative_error)}) where {Name} =
    Units.label(Units.Quantity{Name, Nothing}()) * " relative error"
Units.symbol(::Units.Quantity{Name, typeof(absolute_error)}) where {Name} =
    "Δ" * Units.symbol(Units.Quantity{Name, Nothing}())
Units.symbol(::Units.Quantity{Name, typeof(relative_error)}) where {Name} =
    "ε" * Units.symbol(Units.Quantity{Name, Nothing}())
