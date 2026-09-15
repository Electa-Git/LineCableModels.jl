# Default equivalent homogeneous-earth rule

## Identification and source

| Field | Value |
| --- | --- |
| Family | Equivalent homogeneous earth |
| Formula identifier | `:bottommost`; `:default` routes to this rule |
| Documentation status | Registered and documented. |

**Description.** `:default` routes to the explicit `:bottommost` rule.
The bottommost soil layer supplies the equivalent
homogeneous resistivity, relative permittivity, and relative permeability.
This is a package policy, not an author-named implementation.

**Assumptions.**

The layer property vectors include air at index 1 and soil layers at indices
2 through N.

**Expression.**

``(\rho,\varepsilon_r,\mu_r)_{equivalent}=
 (\rho_N,\varepsilon_{r,N},\mu_{r,N})``.

**Approximation.**

The selected material is passed to the consuming homogeneous-earth formula.

Martins-Britto et al. reported that deep-layer conductivity can predominate in
the cases they studied, but this registration does not implement their
equivalent-conductivity equation.

**Limitations.**

The approximation may be unsuitable for strong conductivity contrasts or
frequency ranges outside the supporting study, and no accuracy claim is made
for selecting permittivity or permeability from the deepest layer.

**Reference.**

A. G. Martins-Britto, F. V. Lopes, and S. R. M. J. Rondineau, *Multilayer
Earth Structure Approximation by a Homogeneous Conductivity Soil for Ground
Return Impedance Calculations*, 2020 (context only; not the implemented
equation).

[Back to the relevant theory overview](../../earth_properties.md)
