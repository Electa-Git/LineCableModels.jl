# Engine-owned line-parameter quantity and unit accessors.
QuantityUnits.quantity(::typeof(Z)) = QuantityUnits.QuantityTag{:series_impedance}()
QuantityUnits.quantity(::typeof(Y)) = QuantityUnits.QuantityTag{:shunt_admittance}()
QuantityUnits.quantity(::typeof(X)) = QuantityUnits.QuantityTag{:reactance}()
QuantityUnits.quantity(::typeof(G)) = QuantityUnits.QuantityTag{:conductance}()
QuantityUnits.quantity(::typeof(B)) = QuantityUnits.QuantityTag{:susceptance}()

QuantityUnits.quantity(::typeof(Z), ::Val{:re}) = QuantityUnits.quantity(R)
QuantityUnits.quantity(::typeof(Z), ::Val{:im}) = QuantityUnits.quantity(X)
function QuantityUnits.quantity(::typeof(Z), ::Val{:abs})
    QuantityUnits.QuantityTag{(:series_impedance, :abs)}()
end
function QuantityUnits.quantity(::typeof(Z), ::Val{:angle})
    QuantityUnits.QuantityTag{(:series_impedance, :angle)}()
end
QuantityUnits.quantity(::typeof(Y), ::Val{:re}) = QuantityUnits.quantity(G)
QuantityUnits.quantity(::typeof(Y), ::Val{:im}) = QuantityUnits.quantity(B)
function QuantityUnits.quantity(::typeof(Y), ::Val{:abs})
    QuantityUnits.QuantityTag{(:shunt_admittance, :abs)}()
end
function QuantityUnits.quantity(::typeof(Y), ::Val{:angle})
    QuantityUnits.QuantityTag{(:shunt_admittance, :angle)}()
end
