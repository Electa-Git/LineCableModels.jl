# Default cable-layer admittivity

## Identification and source

| Field | Value |
| --- | --- |
| Family | Insulation admittance |
| Formula identifier | `:default` |
| Route | `:lossless` |
| Documentation status | Routing alias; the equation is documented separately. |

**Description.** `:default` routes to the explicit `:lossless` constitutive
relation. It is retained as the package-level selection convention and does
not define a second equation.

**Assumptions.**

_To be documented._

**Expression.** See [Lossless cable-layer admittivity](lossless.md).

**Approximation.** Lossless dielectric material: conduction and polarization
loss are suppressed.

**Limitations.** Select `:lossy` to retain material conduction and polarization
loss.

**Reference.** Standard frequency-domain dielectric constitutive relation.

[Back to the relevant theory overview](../insulation_parameters.md)
