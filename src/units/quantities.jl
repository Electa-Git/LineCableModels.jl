"""
$(TYPEDEF)

Represent the scientific meaning selected by `Name` and `Component`. It holds no
numerical values or presentation metadata. `Name` is the type of the quantity's
selector, such as `typeof(Z)`, or a type that its module declares where no selector
exists, such as [`Frequency`](@ref). `Component` is the type of the transform that
derives the quantity from its base quantity, such as `typeof(abs)`, and `Nothing` for
the base quantity itself.
"""
struct Quantity{Name, Component} end

"""
$(TYPEDEF)

The name of the frequency quantity, `Quantity{Frequency, Nothing}`. Commons and the
reports build it, and they load before Engine, which owns the `frequencies` selector.
"""
abstract type Frequency end

"""
$(TYPEDEF)

The name of a quantity without unit or symbol, such as a categorical plot axis.
"""
abstract type Dimensionless end

"""
    quantity(selector::Function)
    quantity(selector::Function, transform::Function)
    quantity(base::Quantity, transform::Function)

Return the typed [`Quantity`](@ref) identity registered for a scientific
selector or selector-transform pair, or the component of a base quantity under
`transform`.

A selector without a registered quantity has no method.
"""
function quantity end

quantity(::Quantity{Name, Nothing}, ::Component) where {Name, Component <: Function} =
    Quantity{Name, Component}()
