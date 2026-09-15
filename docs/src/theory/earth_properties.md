# Earth properties

Frequency-dependent soil relations determine material conductivity and
permittivity. Equivalent homogeneous-earth models (EHEM) instead reduce a
layered earth to effective material values for a separately selected
earth-return formula. Neither operation is itself an impedance or
admittance matrix.

## Frequency-dependent soil relations

Each relation retains its reference resistivity, units, fitted frequency
range, and source assumptions. `:default` is a routing alias for the explicit
`:constant` pass-through.

- [Default frequency-dependent earth material](earth-properties/frequency-dependent/default.md)
- [Registered frequency-dependent soil relations](earth-properties/frequency-dependent/formulas.md)

## Equivalent homogeneous earth

These reductions replace the layered response by a homogeneous
approximation. They must be paired with the earth-return formula and
material reconstruction stated by the source.

- [Default equivalent homogeneous-earth rule](earth-properties/equivalent-homogeneous/default.md)

[Back to Contents](contents.md)
