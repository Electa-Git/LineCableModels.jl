# Dielectrics

```@meta
CurrentModule = LineCableModels.Engine
```

The dielectric layers between a conductor and the next conductor outward form its interval.
The magnetic field in the interval adds a series term to the loop impedance of the conductor
[Ametani1980](@cite), with the relative permeability `mu_r_ins` of the interval. For a
stranded conductor, that permeability includes the solenoid factor of the helical current of
its wires, as [Cable data model](../data-model.md) describes.

Each layer also has an admittivity ``\kappa``, which the shunt formula reads to build the
potential coefficients inside the cable by the series chain of the layers, which the
[matrix formulation](matrix.md) page gives. An insulator takes the law of
`insulation_admittance` and a semicon the law of `semicon_admittance`, by the class of its
material.
The one law of this release, `:lossless`, leaves out conduction and polarization loss in both.

## Insulation impedance

```@autodocs
Modules = [InsulationImpedance]
Order = [:module, :constant, :type, :function]
Public = true
Private = false
```

## Insulation admittance

```@autodocs
Modules = [InsulationAdmittance]
Order = [:module, :constant, :type, :function]
Public = true
Private = false
```

## Semicon admittance

```@autodocs
Modules = [SemiconAdmittance]
Order = [:module, :constant, :type, :function]
Public = true
Private = false
```

[Back to Contents](contents.md)
