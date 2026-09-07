# Maaouni–Amri two-half-space potential/admittance approximation

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Parallel thin bare overhead wires; scalar self/mutual coefficients; no coating. |
| Calculated quantities | Thin-wire potential coefficient and qTEM admittance near a lossy interface |
| Earth structure | Homogeneous earth below lossless air. |
| Model and approximation | qTEM limit of the air transverse constant, followed by the exponential approximation (23) used in (25)–(26); no independent full-wave modal solve. |
| Main source | A. Maaouni, A. Amri, and A. Zouhir (2001) |
| Citation key(s) | `:Maaouni2001` |
| Evidence status | Original publication page images checked |

**Description.** Potential/admittance companion of the source's two-half-space qTEM wire formulation.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | The reduction sets the air transverse constant ``τ≈0``, hence the axial wavenumber to ``k_0``; it does not set the axial wavenumber to zero. | Section 3. |
| Air propagation constant ``γ_air`` | The source uses the real air wavenumber ``k_0=ω√(μ_0ε_0)``. | Section 2. |
| Earth propagation constant ``γ_earth`` | Complex medium ratio ``n``. | Parent equations. |
| Earth permittivity and displacement current | Retained through complex material constants. | Definitions. |
| Range of validity | Infinite thin overhead wires; both ``Re[z√(1−n²)]`` and ``Re[z̄√(1−n²)]`` must be positive. | Sections 2–3. |
| Earth permeability ``μ_earth`` | Medium ratio retained. | Definitions. |
| Arrangement | Scalar self/mutual wire coefficients. | Formulation. |
| Earth structure | Homogeneous earth below lossless air. | Figure 1. |
| Conductor and insulation geometry | Thin bare wire; no coating. | Model. |
| Constitutive and field assumptions | Linear isotropic harmonic media. | Derivation. |
| Conventions | Matrix admittance requires inversion after coefficient assembly. | Source TL assembly. |

**Expression.** The source's potential kernel uses the same analytical helper

```math
G(X,Y)\simeq\frac{n^2}{2(n^4-1)}[Q(bz)+Q(b\bar z)]
-\frac{P(b,z)+P(b,\bar z)-P(-b,z)-P(-b,\bar z)-n^2b[Q(-bz)+Q(-b\bar z)]}{2b(n^4-1)},
\qquad\text{(26)}
```

with ``z=k_0(Y+jX)``, ``b=j/\sqrt{1+n^2}``, ``Q(z)=e^{-z}E_1(-z)``, and ``P`` from (25). The assembled potential-coefficient matrix is inverted as a matrix to obtain the shunt parameters.

**Approximation.** qTEM/asymptotic analytical approximation; not a full-wave solved ``Γ``.

**Limitations.** Overhead thin wires above nonmagnetic homogeneous earth. Equation (26) is evaluated only within its stated transformed-argument domain; no buried or mixed coefficient is inferred.

**Reference.** [Maaouni2001](@cite), equations (2)–(7), (22)–(26).

**Transcription source.** Original page-image verification of the full ``G`` expression.

## Source transcription

The source uses ``e^{-j\omega t}``,
``n^2=\varepsilon_g/\varepsilon_0+j\sigma_g/(\omega\varepsilon_0)``,
and ``k_0=\omega\sqrt{\mu_0\varepsilon_0}``. Both half-spaces have
permeability ``\mu_0``; all wire axes are in air. The qTEM reduction
sets ``\tau=\sqrt{k_z^2-k_0^2}\simeq0``, not ``k_z=0``.
Numerical values under the survey's ``e^{j\omega t}`` convention are
the complex conjugates of the source expressions.

The analytical helper omitted from the earlier transcription is

```math
P(b,z)\simeq-\frac{1-n^2}{2b}
\left[
\ln\left(1+\frac{2}{z\sqrt{1-n^2}}\right)
+Q\left(bz+\frac{2b}{\sqrt{1-n^2}}\right)
\right]
+\frac{1}{z}+bQ(bz)\left(1+\frac{1-n^2}{2b^2}\right).
\qquad\text{(25)}
```

Equation (26) requires both
``\Re[z\sqrt{1-n^2}]>0`` and
``\Re[\bar z\sqrt{1-n^2}]>0``. The source justifies this condition
for its usual overhead geometry ``Y>X`` and large ``|n|``;
it is checked directly when evaluating the approximation.

Under ``e^{jωt}``, let ``G_+`` denote the conjugated value of (26). The potential coefficient is ``P_{mn}=[ln(ρ^*_{mn}/ρ_{mn})+G_+]/(2πε_0)``. It follows from (2) as ``τa_nK_1(τa_n)→1``. For a self term, the direct distance is the wire radius and the image distance is twice the axis height. The external potential matrix is assembled before inversion; ``G`` is neither an admittance entry nor the complete potential coefficient.

The unapproximated qTEM parent (7) is the overhead integral represented by [Wise1948](@cite). Equation (26), with the explicit helper (25), remains a distinct analytical approximation to that parent.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``G`` | unchanged | electric/potential correction kernel | source normalized |
| ``z,\bar z`` | unchanged | complex separation arguments | dimensionless |
| ``P,Q`` | unchanged | exponential-integral helpers | source convention |

## Evidence and approximation sources

The formulation was found during the broad lossy-dielectric/interface sweep.

## Limitations and discrepancies

The source's scalar exposition is not converted into entrywise inversion of matrix coefficients.
