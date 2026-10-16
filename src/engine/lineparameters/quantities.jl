# Engine-owned line-parameter quantity and unit accessors.
Units.quantity(::typeof(frequencies)) = Units.Quantity{Units.Frequency, Nothing}()
Units.quantity(selector::Union{typeof(Z), typeof(Y), typeof(X), typeof(G), typeof(B)}) =
    Units.Quantity{typeof(selector), Nothing}()

Units.native_unit(::Units.Quantity{typeof(Z), Nothing}) =
    Units.units(:base, :ohm; per = (:base, :meter))
Units.display_unit(::Units.Quantity{typeof(Z), Nothing}) =
    Units.units(:base, :ohm; per = (:kilo, :meter))
Units.label(::Units.Quantity{typeof(Z), Nothing}) = "Series impedance"
Units.symbol(::Units.Quantity{typeof(Z), Nothing}) = "Z"

Units.native_unit(::Units.Quantity{typeof(X), Nothing}) =
    Units.units(:base, :ohm; per = (:base, :meter))
Units.display_unit(::Units.Quantity{typeof(X), Nothing}) =
    Units.units(:base, :ohm; per = (:kilo, :meter))
Units.label(::Units.Quantity{typeof(X), Nothing}) = "Series reactance"
Units.symbol(::Units.Quantity{typeof(X), Nothing}) = "X"

Units.native_unit(::Units.Quantity{typeof(Y), Nothing}) =
    Units.units(:base, :siemens; per = (:base, :meter))
Units.display_unit(::Units.Quantity{typeof(Y), Nothing}) =
    Units.units(:base, :siemens; per = (:kilo, :meter))
Units.label(::Units.Quantity{typeof(Y), Nothing}) = "Shunt admittance"
Units.symbol(::Units.Quantity{typeof(Y), Nothing}) = "Y"

Units.native_unit(::Units.Quantity{typeof(G), Nothing}) =
    Units.units(:base, :siemens; per = (:base, :meter))
Units.display_unit(::Units.Quantity{typeof(G), Nothing}) =
    Units.units(:base, :siemens; per = (:kilo, :meter))
Units.label(::Units.Quantity{typeof(G), Nothing}) = "Shunt conductance"
Units.symbol(::Units.Quantity{typeof(G), Nothing}) = "G"

Units.native_unit(::Units.Quantity{typeof(B), Nothing}) =
    Units.units(:base, :siemens; per = (:base, :meter))
Units.display_unit(::Units.Quantity{typeof(B), Nothing}) =
    Units.units(:base, :siemens; per = (:kilo, :meter))
Units.label(::Units.Quantity{typeof(B), Nothing}) = "Shunt susceptance"
Units.symbol(::Units.Quantity{typeof(B), Nothing}) = "B"

Units.family(::Units.Quantity{<:Union{typeof(Z), typeof(X)}, Nothing}) = Val(:series)
Units.family(::Units.Quantity{<:Union{typeof(Y), typeof(G), typeof(B)}, Nothing}) = Val(:shunt)

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
