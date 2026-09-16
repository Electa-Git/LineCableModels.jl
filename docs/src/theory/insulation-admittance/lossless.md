# Lossless cable-layer admittivity

## Identification

| Field | Value |
| --- | --- |
| Family | Insulation admittance and semicon admittance |
| Formula identifier | `:lossless` |
| Compact description | `Lossless` |

**Description.** Lossless frequency-domain dielectric admittivity retaining
displacement current while suppressing material conduction and polarization
loss.

```math
\\kappa=j\\omega\\varepsilon_0\\varepsilon_r.
```

The common radial annulus operator converts this material coefficient to
``Y=2\\pi\\kappa/\\ln(b/a)``. The `:default` selection is only a routing alias
to this equation.

## Scope

The same explicit equation is available in the insulation and semicon
admittance boundaries. The selected boundary still determines which material
layers receive it.

[Back to insulation parameters](../insulation_parameters.md)
