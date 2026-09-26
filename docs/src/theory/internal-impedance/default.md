# Default cylindrical-conductor surface impedances

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Formula identifier | `:schelkunoff1934`; `:default` routes to this implementation |
| Explicit literature identifier | `:schelkunoff1934` |
| Documentation status | Registered and documented. |

**Description.** Exact cylindrical surface impedances for solid and hollow
round conductors. The package default and `:schelkunoff1934` use the same
Schelkunoff route; the latter exposes the author-year identity.

**Assumptions.**

The conductor is homogeneous, concentric, and linear. The conductor
permeability and resistivity are finite and positive.

**Expression.**

For inner radius ``a``, outer radius ``b``, and
``m=\sqrt{j\omega\mu/\rho}``, the surface terms use modified Bessel
functions and the denominator
``D=I_1(mb)K_1(ma)-K_1(mb)I_1(ma)``. The solid-cylinder limit is used when
``a=0``.

**Approximation.**

The radial solution is exact for the cylindrical conductor model and is
assembled recursively for concentric conductive terminals.

**Limitations.**

It does not model arbitrary proximity-induced angular current redistribution.

**Reference.**

S. A. Schelkunoff, *The Electromagnetic Theory of Coaxial Transmission Lines
and Cylindrical Shields*, 1934; A. Ametani, *A General Formulation of
Impedance and Admittance of Cables*, 1980.

[Back to the relevant theory overview](../internal_impedance.md)
