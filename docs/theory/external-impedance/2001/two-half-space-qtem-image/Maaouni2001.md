# Maaouni–Amri two-half-space qTEM impedance approximation

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel thin bare overhead wires; lateral separation ``X`` and summed height ``Y``; no finite coating. |
| Calculated quantities | Thin-wire longitudinal impedance near a lossy interface and qTEM logarithmic image approximation |
| Earth structure | Homogeneous earth below lossless air. |
| Model and approximation | Equation (12) and the closed ``G`` form are qTEM/asymptotic reductions of the preceding exact spectral representation. |
| Main source | A. Maaouni, A. Amri, and A. Zouhir (2001) |
| Citation key(s) | `:Maaouni2001` |
| Evidence status | Original publication page images checked |

**Description.** Spectral parent and analytical qTEM image approximation for parallel thin overhead wires above homogeneous earth.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | The reduction sets the air transverse constant ``τ≈0``, hence the axial wavenumber to ``k_0``; it does not set the axial wavenumber to zero. | Section 3. |
| Air propagation constant ``γ_air`` | The source uses the real air wavenumber ``k_0=ω√(μ_0ε_0)``. | Section 2. |
| Earth propagation constant ``γ_earth`` | Ratio ``n`` represents the second medium. | (2)–(7). |
| Earth permittivity and displacement current | Complex material constants retained. | Constitutive definitions. |
| Range of validity | Thin infinite overhead wires. The image approximation uses ``τ≈0``; the companion closed potential kernel also requires both transformed separation arguments to have positive real parts. | Sections 2–3. |
| Earth permeability ``μ_earth`` | Both media have vacuum permeability. | Section 2 and Figure 1. |
| Arrangement | Self and mutual distances ``X,Y``. | (2)–(12). |
| Earth structure | Homogeneous earth below lossless air. | Figure 1. |
| Conductor and insulation geometry | Thin bare wire; no finite coating. | Model statement. |
| Constitutive and field assumptions | Linear isotropic media and harmonic fields. | Derivation. |
| Conventions | The source uses ``e^{-jωt}``; results are conjugated for ``e^{jωt}``. Decaying roots and principal exponential-integral branches are retained. | Section 2; (12), (25)–(26). |

**Expression.** The qTEM impedance image integral reduces as

```math
\begin{aligned}
J(X,Y)&\simeq\ln\frac{\rho_J^*}{\rho^*} \\
\rho_J^*&=\sqrt{X^2+(Y+Y_J)^2} \\
\rho^*&=\sqrt{X^2+Y^2} \\
Y_J&=\frac{2}{k_0\sqrt{1-n^2}}.
\end{aligned}\qquad\text{(12)}
```

The full closed evaluator for the companion kernel is

```math
G(X,Y)\simeq\frac{n^2}{2(n^4-1)}[Q(bz)+Q(b\bar z)]
-\frac{P(b,z)+P(b,\bar z)-P(-b,z)-P(-b,\bar z)-n^2b[Q(-bz)+Q(-b\bar z)]}{2b(n^4-1)},
\qquad\text{(26)}
```

where ``z=k_0(Y+jX)``, ``b=j/\sqrt{1+n^2}``, ``Q(z)=e^{-z}E_1(-z)``, and ``P`` is defined by (25).

**Approximation.** Equation (12) and the closed ``G`` form are qTEM/asymptotic reductions of the preceding exact spectral representation.

**Limitations.** Overhead thin wires above nonmagnetic homogeneous earth. Neither buried nor mixed conductor pairs follow from the extracted formulas.

**Reference.** [Maaouni2001](@cite), equations (2)–(26).

**Transcription source.** Original page images; conjugates, signs, and denominator ``n^4-1`` verified.

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

The longitudinal correction is ``J``. In the qTEM limit, the external series coefficient is ``jωμ_0[ln(ρ^*/ρ)+J]/(2π)`` under ``e^{jωt}``. Substituting (12) gives the air-referenced complex-image coefficient also written in [Ametani2014](@cite), (17)–(18). The historical image construction remains attributed to [Dubanton1969](@cite); the displacement-retaining witness here predates the 2014 restatement. The helper ``G`` belongs to the separate potential record and is not added to the series coefficient.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``X,Y`` | unchanged | lateral/normal separation variables | m |
| ``n`` | unchanged | medium parameter ratio | complex |
| ``J,G`` | unchanged | spectral correction kernels | source normalized |

## Evidence and approximation sources

The equation groups provide both impedance and potential/admittance outputs.

## Limitations and discrepancies

The source does not provide a general finite-layer recursion.
