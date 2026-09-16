# Default frequency-dependent earth material

## Identification and source

| Field | Value |
| --- | --- |
| Family | Frequency-dependent earth properties |
| Formula identifier | `:default` |
| Route | `:constant` |
| Documentation status | Routing alias; the pass-through is documented separately. |

**Description.** `:default` routes to the explicit `:constant` static
earth-material pass-through. The registered frequency-dependent relations are
listed in [Registered frequency-dependent soil relations](formulas.md).

**Assumptions.**

The input material is homogeneous and its properties are already represented
at the requested frequency.

**Expression.**

``\rho(f)=\rho_0``, ``\varepsilon_r(f)=\varepsilon_{r,0}``, and
``\mu_r(f)=\mu_{r,0}``.

**Approximation.**

No frequency-dependent constitutive correction is applied.

**Limitations.**

Selecting `:default` does not model soil dispersion. Select `:constant` for the
same pass-through explicitly, or select one of the registered empirical laws.

**Reference.**

Package routing default; no author equation is claimed by this registration.

[Back to the relevant theory overview](../../earth_properties.md)
