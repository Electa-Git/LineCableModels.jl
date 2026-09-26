# Default modal decomposition and eigenvalue tracking

## Identification and source

| Field | Value |
| --- | --- |
| Family | Modal transformation |
| Formula identifier | `:default` |
| Explicit literature identifier | `:chrysochos2014` |
| Documentation status | Registered and documented. |

**Description.** Levenberg–Marquardt tracking of complex modal eigenpairs,
initialized from the preceding frequency. The package default and
`:chrysochos2014` use the same route; the latter exposes the author-year
identity.

**Assumptions.**

The phase-domain impedance and admittance matrices are fully coupled and the
frequency samples are ordered.

**Expression.**

The scaled matrix is
``\widetilde{S}=YZ/(-\omega^2\mu_0\varepsilon_0)-I``. Each eigenpair is
tracked by a real least-squares residual with ``t^Tt=1``.

**Approximation.**

An analytic real Jacobian and damped normal-equation step track the modes;
matched conventional eigensolutions are retained when iteration fails.

**Limitations.**

The route depends on frequency ordering and its iteration controls. It is not
a fixed symmetrical-component transform.

**Reference.**

A. I. Chrysochos, T. A. Papadopoulos, and G. K. Papagiannis, *Robust
Calculation of Frequency-Dependent Transmission-Line Transformation Matrices
Using the Levenberg–Marquardt Method*, 2014.

[Back to the relevant theory overview](../modal_decomposition.md)
