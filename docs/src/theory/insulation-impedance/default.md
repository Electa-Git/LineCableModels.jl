# Default coaxial-insulation magnetic series impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | Insulation impedance |
| Formula identifier | `:default` |
| Explicit literature identifier | `:ametani1980` |
| Documentation status | Registered and documented. |

**Description.** Longitudinal magnetic impedance of one concentric insulation
region. The package default and `:ametani1980` use the same Ametani route; the
latter exposes the author-year identity.

**Assumptions.**

The region is a homogeneous annulus with finite positive relative
permeability.

**Expression.**

For inner and outer radii ``a`` and ``b``,
``Z_{ins}=j\omega\mu_0\mu_r\ln(b/a)/(2\pi)``.

**Approximation.**

The term is assembled with conductor surface impedances in the cable series
impedance matrix.

**Limitations.**

The term vanishes for zero-thickness or zero-inner-radius regions under the
package's boundary convention.

**Reference.**

 A. Ametani, *A General Formulation of Impedance and Admittance of Cables*,
1980.

[Back to the relevant theory overview](../insulation_parameters.md)
