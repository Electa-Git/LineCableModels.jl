# Shared physical selectors extended independently by DataModel and Engine.
quantity(::typeof(R)) = Quantity{:series_resistance}()
quantity(::typeof(L)) = Quantity{:series_inductance}()
quantity(::typeof(C)) = Quantity{:shunt_capacitance}()
