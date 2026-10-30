# Conductors

```@meta
CurrentModule = LineCableModels.Engine
```

A coaxial formulation reads each terminal of a cable as a tubular conductor of inner radius
``r_\mathrm{in}`` and outer radius ``r_\mathrm{ex}``, a solid conductor when
``r_\mathrm{in} = 0``. The internal impedance of the tube is the electromagnetic response
inside the metal, given by three surface terms per unit length
[Schelkunoff1934](@cite): the inner surface impedance, the outer surface impedance and the
transfer impedance through the wall. The loop recursion of the
[matrix formulation](matrix.md) adds the outer surface impedance of a conductor and the inner
surface impedance of the next conductor outward to the loop impedance of the conductor. The
loop impedance less the transfer impedance adds to the entries between the conductor and each
conductor inside it.

A stranded conductor or an armor enters as the equivalent tube of its terminal, which the data
model builds: the tube spans the conductive regions and takes their direct-current resistance
and geometric mean radius. A strand on a helix is longer than the cable by its overlength, so
the resistance of the tube includes it, as [Cable data model](../data-model.md) describes.

The resistivity of each conductor follows from its material at the operating temperature of
the problem by the resistivity law, before the surface impedances are computed.

## Internal impedance

```@autodocs
Modules = [InternalImpedance]
Order = [:module, :constant, :type, :function]
Public = true
Private = false
```

## Resistivity

```@autodocs
Modules = [Resistivity]
Order = [:module, :constant, :type, :function]
Public = true
Private = false
```

[Back to Contents](contents.md)
