# Shared physical selectors extended independently by DataModel and Engine.
quantity(::typeof(R)) = QuantityTag{:series_resistance}()
quantity(::typeof(L)) = QuantityTag{:series_inductance}()
quantity(::typeof(C)) = QuantityTag{:shunt_capacitance}()
