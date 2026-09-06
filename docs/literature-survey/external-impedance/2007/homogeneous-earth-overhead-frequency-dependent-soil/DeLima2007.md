# de Lima–Portela overhead impedance with frequency-dependent soil

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel overhead conductors represented by heights and separation. No conductor radius or insulation parameter enters this ground correction. |
| Calculated quantities | Self and mutual overhead ground-return corrections per unit length |
| Earth structure | Linear homogeneous isotropic half-space with a planar air interface. |
| Model and approximation | Integral representation within the declared TEM/quasi-TEM model, with a separate empirical constitutive approximation (7)–(8). The authors evaluate the integrals by Gauss–Kronrod quadrature; the numerical evaluation is not a new physical kernel or an analytical asymptote. In the conduction-only limit the source sets ``\varepsilon_s=0`` and ``\sigma_s=\sigma_0``, hence ``\eta_s=\sqrt{j\omega\mu_0\sigma_0}``. |
| Main source | Antonio Carlos Siqueira de Lima and Carlos Portela (2007), extending Carson-type expressions to complex frequency-dependent soil parameters; the constitutive soil model is attributed to earlier Portela work |
| Citation key(s) | `:DeLima2007` |
| Evidence status | Original PDF equations checked against page images; longitudinal-reduction prescription unresolved |

**Description.** Self and mutual ground-return impedance corrections for parallel overhead conductors above homogeneous earth with coupled frequency-dependent conductivity and permittivity.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Appendix assumes ``\exp(-kz)`` with **unknown** ``k``. Final (1)–(2) retain no ``k``. The step eliminating it is not specified; this is not evidence that an impressed ``Γ=0`` was explicitly prescribed. | Stated — appendix opening, p. 497; equation-implied — (1)–(2), p. 493. |
| Air propagation constant ``γ_air`` | No air propagation constant is retained in (1)–(2); appendix air equation is ``\nabla^2E_a=0``. Its relation to the preceding unknown longitudinal ``k`` is unresolved. | Equation-implied — (1)–(2), p. 493; (15), p. 497. |
| Earth propagation constant ``γ_earth`` | ``\eta_s=\sqrt{j\omega\mu_0(\sigma_s+j\omega\varepsilon_s)}``. Square-root branch not separately stated. | Stated — definition after (2), p. 493. |
| Earth permittivity and displacement current | Retained. ``\sigma_s+j\omega\varepsilon_s`` is modeled jointly as a complex soil immittance ``\kappa'``; conductivity and dielectric response must not be varied independently of the chosen model. | Stated — (7)–(9) and section II-B, pp. 493–494. |
| Range of validity | Authors discuss acceptable TEM/quasi-TEM error up to about 1 MHz depending on soil, geometry, and line length; recommend other field approaches for very short lines or frequencies above some MHz. Soil-fit agreement up to 2 MHz is a separate measured-model claim. | Stated — section II opening, p. 493; section II-B, p. 494. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0`` in the final impedance and propagation constant. | Equation-implied — (1)–(2) and ``\eta_s`` definition. |
| Arrangement | Overhead, self and mutual; heights ``h_i,h_j`` and horizontal separation ``d_{ij}``. | Stated — section II-A, p. 493. |
| Earth structure | Linear homogeneous isotropic half-space with a planar air interface. | Stated — section II opening and appendix assumptions, pp. 493, 497. |
| Conductor and insulation geometry | Parallel overhead conductors represented by heights and separation. No conductor radius or insulation parameter enters this ground correction. | Equation-implied — (1)–(2); radius/internal and ideal external contributions are separate in (11), p. 494. |
| Constitutive and field assumptions | Linear, homogeneous, isotropic media; source TEM or quasi-TEM approximation; frequency-dependent scalar soil response. No additional approximation of the displayed integral kernel is introduced by the authors' quadrature method. | Stated — section II and appendix; quadrature paragraph, p. 493. |
| Conventions | ``\exp(j\omega t)`` time dependence and ``\exp(-kz)`` longitudinal dependence. Heights positive above the interface; ``\xi`` is the transverse spectral variable. Output is a per-unit-length ground correction. | Stated — appendix opening, p. 497; definitions on p. 493 and decomposition (11), p. 494. |

**Expression.** Source self and mutual correction integrals (1)–(2).

```math
z_{c_{ii}}=\frac{j\omega\mu_0}{\pi}
\int_0^\infty
\frac{\exp(-2h_i\xi)}{\xi+\sqrt{\xi^2+\eta_s^2}}\,d\xi,
\qquad\text{(1)}

z_{c_{ij}}=\frac{j\omega\mu_0}{\pi}
\int_0^\infty
\frac{\exp(-(h_i+h_j)\xi)}{\xi+\sqrt{\xi^2+\eta_s^2}}
\cos(d_{ij}\xi)\,d\xi,
\qquad\text{(2)}

\eta_s=\sqrt{j\omega\mu_0(\sigma_s+j\omega\varepsilon_s)}.
```

``h_i,h_j,d_{ij}`` are lengths in meters, ``\xi,\eta_s`` inverse lengths, ``\sigma_s`` soil conductivity in S/m, and ``\varepsilon_s`` soil permittivity in F/m. The source's joint soil model is

```math
\sigma_s+j\omega\varepsilon_s\simeq\kappa'
=\sigma_0+\delta_{\sigma_s}+j\delta_{\omega\varepsilon_s},
\qquad\text{(7)}

\delta_{\sigma_s}+j\delta_{\omega\varepsilon_s}
=\Delta_i\left(\frac{f}{10^6}\right)^\alpha
\left(\cot(\alpha\pi/2)+j\right).
\qquad\text{(8)}
```

``\sigma_0`` is low-frequency conductivity, ``\Delta_i`` is in S/m, and ``\alpha`` controls frequency dependence. ``f`` is frequency and ``\omega=2\pi f``. ``\Delta_i`` is the source's soil-fit parameter, not a conductor-indexed geometric quantity. The source-prescribed line assembly is

```math
Z=Z_i+Z_{\mathrm{ext}}+Z_g,
\qquad\text{(11)}
```

where the matrices are conductor internal, ideal external, and ground-return impedance per unit length. Equations (1)–(2) supply ``Z_g``; they are not the complete external impedance.

**Approximation.** Integral representation within the declared TEM/quasi-TEM model, with a separate empirical constitutive approximation (7)–(8). The authors evaluate the integrals by Gauss–Kronrod quadrature; the numerical evaluation is not a new physical kernel or an analytical asymptote. In the conduction-only limit the source sets ``\varepsilon_s=0`` and ``\sigma_s=\sigma_0``, hence ``\eta_s=\sqrt{j\omega\mu_0\sigma_0}``.

**Limitations.** No explicit derivation step fixes or eliminates the appendix's unknown longitudinal ``k``. Branch prescriptions are not separately stated. The source leaves the corresponding frequency-dependent ground-admittance extension to future research; it must not be inferred by applying the impedance's soil substitution to an unrelated shunt formula.

**Reference.** [DeLima2007](@cite), (1)–(2), (7)–(8), p. 493; (9), (11), p. 494; appendix opening, p. 497.

**Transcription source.** Original PDF images, pp. 493–494 and 497–498. The complete kernels, prefactors, integration domains, soil definition, signs, and model parameters were checked visually. PDF text extraction omits most equations and was used for navigation only.

## Source transcription

Equations (1)–(2), (7)–(8), and (11) above retain the source's notation and decomposition. The illustrative soil fit in (9) is also printed:

```math
\kappa'=84.16\,10^{-6}
+\omega^{0.71603}(0.057849+j0.12097)10^{-6}.
\qquad\text{(9)}
```

This particular sample uses ``\sigma_0=84.16\,\mu\mathrm S/\mathrm m``, ``\Delta_i=8.92028\,\mathrm{mS}/\mathrm m``, and ``\alpha=0.71603``. Equation (9) takes ``\omega`` in rad/s and returns ``\kappa'`` in S/m. It is an example of (7)–(8), not a universal earth law. The source credits earlier soil measurement/model work in references [3]–[5], [10]–[11]; no new priority is assigned to this supporting constitutive model.

## Notation map

Notation is unchanged.

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``z_{c_{ii}},z_{c_{ij}}`` | unchanged | Overhead self/mutual earth corrections | Ω/m |
| ``\eta_s`` | unchanged | Soil bulk propagation constant | inverse m |
| ``\xi`` | unchanged | Transverse integration variable | inverse m; nonnegative domain |
| ``h_i,h_j,d_{ij}`` | unchanged | Heights and horizontal separation | m |
| ``\sigma_s,\varepsilon_s,\mu_0`` | unchanged | Soil conductivity, permittivity, fixed permeability | S/m, F/m, H/m |
| ``\kappa',\sigma_0`` | unchanged | Complex soil immittance and low-frequency conductivity | S/m |
| ``\delta_{\sigma_s},\delta_{\omega\varepsilon_s},\Delta_i,\alpha`` | unchanged | Conductive/displacement increments, fitted amplitude, exponent | S/m for first three; exponent dimensionless |
| ``Z,Z_i,Z_{\mathrm{ext}},Z_g`` | unchanged | Complete/internal/ideal-external/ground impedance matrices | Ω/m |
| ``k,z,t,j,\omega,f`` | unchanged | Unknown longitudinal propagation function, axial coordinate, time, imaginary unit, angular frequency, frequency | inverse m, m, s, dimensionless, rad/s, Hz |

## Evidence and approximation sources

The exact material change is the joint complex soil response inside ``\eta_s``. The source explicitly compares this treatment with conventional constant ``\sigma_0`` and does not reinstate an air wave term in the final kernel. Model validity and experimental soil-fit range are distinct statements on pp. 493–494. The numerical quadrature scheme is described only to identify how the authors evaluated the physical integrals; no implementation is added.

## Limitations and discrepancies

- The appendix introduces an unknown longitudinal ``k`` but final transverse equations lack it; its reduction cannot be reconstructed from an assumed ``Γ=0``.
- The appendix definition after (16) prints ``\eta_g=j\omega\mu(\sigma+j\omega\varepsilon)`` without the square root present in main-text ``\eta_s``. This is a source inconsistency; the main-text definition is retained for (1)–(2).
- Appendix (14) omits a ``j`` inside its material parenthesis compared with (13) and the main-text propagation definition. No derived repair is supplied.
- The original article's DOI is recorded from the source identity; the existing BibTeX entry is copied without filling its missing DOI.
