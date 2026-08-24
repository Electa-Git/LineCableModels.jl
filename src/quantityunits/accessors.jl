# Shared physical accessors extended independently by DataModel and Engine.
quantity(::typeof(R)) = QuantityTag{:resistance}()
quantity(::typeof(L)) = QuantityTag{:inductance}()
quantity(::typeof(C)) = QuantityTag{:capacitance}()
