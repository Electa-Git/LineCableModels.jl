# De Lima et al. quasi-full-wave single-conductor impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinitely long circular conductor of radius ``r`` at signed region-1 coordinate ``h``; buried application is bare. |
| Calculated quantities | Full-wave parent per-unit-length impedance and noniterative quasi-full-wave evaluation for one overhead or bare buried conductor parallel to a planar interface |
| Earth structure | Two homogeneous half-spaces separated by a plane interface. |
| Model and approximation | The qFW operation is the source's substitution of an image-approximation propagation constant ``\bar\gamma`` into ``u_i`` and ``\eta_1``. It avoids Newton iteration of the full-wave modal equation (7). It is distinct from qTEM, which sets ``\gamma\approx\gamma_1`` and uses a leading logarithmic term, and from closed-form image evaluation of the Sommerfeld integrals. |
| Main source | A. C. S. de Lima, A. P. C. Magalhães, P. E. D. Rocha, R. A. Meyberg, and M. T. C. de Barros (2018) |
| Citation key(s) | `:DeLima2018` |
| Evidence status | Original publication page images checked |

**Description.** Per-unit-length impedance of a single infinitely long thin conductor in either of two lossy homogeneous half-spaces, evaluated from a full-wave scalar/vector-potential formulation or with the source's noniterative quasi-full-wave longitudinal propagation estimate.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Source ``\gamma`` is an unknown modal propagation constant with ``e^{-\gamma z}``. Full wave solves (7) iteratively. The qFW formula instead prescribes ``\bar\gamma`` from an image approximation and evaluates spectral roots with it. | Stated — (1), (7)–(9), and (15)–(16), pp. 1874–1875. |
| Air propagation constant ``γ_air`` | Whichever region is air uses ``\gamma_i=\sqrt{j\omega\mu_i(\sigma_i+j\omega\varepsilon_i)}``; the displayed model permits both media to be lossy, rather than setting air to zero. | Stated — below (1), p. 1874. |
| Earth propagation constant ``γ_earth`` | Same indexed definition ``\gamma_i=\sqrt{j\omega\mu_i(\sigma_i+j\omega\varepsilon_i)}``; earth is medium 2 for overhead and medium 1 for buried use. | Stated — p. 1874 and §III-B, p. 1876. |
| Earth permittivity and displacement current | Retained in each medium through ``\sigma_i+j\omega\varepsilon_i``. | Stated — p. 1874. |
| Range of validity | Thin infinite conductor and one planar interface; conductor loss inclusion assumes ``|\gamma_c|\gg|\gamma|``. The qFW was tested for a 100 kHz–100 MHz overhead example and a buried example, but the paper reports best agreement mainly below a few MHz and mode-dependent transitions rather than a universal bound. | Stated — pp. 1874–1878. |
| Earth permeability ``μ_earth`` | Fixed: ``\mu_1=\mu_2=\mu_0``. | Stated — below (1), p. 1874. |
| Arrangement | Single overhead conductor or single bare buried conductor, parallel to the interface; self only. | Stated — abstract, Fig. 1, and §III-A/B. |
| Earth structure | Two homogeneous half-spaces separated by a plane interface. | Stated — Fig. 1 and §II, p. 1874. |
| Conductor and insulation geometry | Infinitely long circular conductor of radius ``r`` at signed region-1 coordinate ``h``; buried application is bare. | Stated — Fig. 1 and §II-A, p. 1874. |
| Constitutive and field assumptions | Linear isotropic lossy media, equal permeability, thin-wire source. Full wave retains longitudinal coupling; qFW changes only the propagation constant supplied to the spectral functions. Internal impedance ``z_i`` is separate. | Stated — (1)–(16), pp. 1874–1875. |
| Conventions | ``e^{j\omega t}`` and ``e^{-\gamma z}``; conductor axis ``z``; ``h>0`` in medium 1 as drawn; per-unit-length ``Z``. | Stated — (1) and surrounding text, p. 1874. |

**Expression.** Full-wave impedance (12), evaluated in the qFW formula with the prescribed substitutions (15) and modal approximation (16), printed p. 1875.

```math
Z=\gamma Z_c=\frac{j\omega\mu_0}{2\pi}
\left[\Lambda_1+S_1-\left(\frac{\gamma}{\gamma_1}\right)^2(S_2+S_4)\right],
\qquad\text{(12)}
```

```math
\Lambda=K_0(\eta_1d)-K_0(\eta_1D),\qquad
d=\sqrt{(h-y)^2+x^2},\quad D=\sqrt{(h+y)^2+x^2},
```

```math
S_1=2\int_0^\infty\frac{e^{-(h+y)u_1}}{u_1+u_2}\cos(x\lambda)\,d\lambda,
\qquad
S_2=2\int_0^\infty\frac{e^{-(h+y)u_1}}{n^2u_1+u_2}\cos(x\lambda)\,d\lambda,
\qquad\text{(4)}
```

```math
S_4=2\int_0^\infty\frac{u_2}{u_1}
\frac{e^{-hu_1}-e^{-2hu_1}}{n^2u_1+u_2}\cos(r\lambda)\,d\lambda,
\qquad\text{(14)}
```

```math
n=\frac{\gamma_2}{\gamma_1},\qquad
u_i=\sqrt{\lambda^2+\gamma_i^2-\gamma^2},\qquad
\eta_1^2=\gamma_1^2-\gamma^2.
```

For qFW, every occurrence in these functions uses

```math
u_i\approx\bar u_i=\sqrt{\lambda^2+\gamma_i^2-\bar\gamma^2},
\qquad
\eta_1\approx\bar\eta=\sqrt{\lambda^2-\bar\gamma^2},
\qquad\text{(15)}
```

where ``\bar\gamma`` is the propagation constant calculated by the cited image approximation and satisfies the noniterative approximate modal equation

```math
\bar M=\frac{2\pi}{j\omega\mu}z_i+
\left(1-\frac{\bar\gamma^2}{\gamma_1^2}\right)\bar\Lambda+
2\left(\bar S_1-\frac{\bar\gamma^2}{\gamma_1^2}\bar S_2\right)=0.
\qquad\text{(16)}
```

**Approximation.** The qFW operation is the source's substitution of an image-approximation propagation constant ``\bar\gamma`` into ``u_i`` and ``\eta_1``. It avoids Newton iteration of the full-wave modal equation (7). It is distinct from qTEM, which sets ``\gamma\approx\gamma_1`` and uses a leading logarithmic term, and from closed-form image evaluation of the Sommerfeld integrals.

**Limitations.** Only a single conductor and single interface are covered. The source does not print the closed image formula that yields ``\bar\gamma`` in this paper, referring to earlier work; the prescription is therefore source-dependent and not independently evaluable from this record. Multiple modal roots and a fast-wave/TL-mode transition can make the choice of full-wave root discontinuous. The qFW tests do not establish a universal frequency range. The paper's derivative (9) is printed with partial derivatives with respect to ``\lambda`` while the surrounding text says differentiation with respect to ``\gamma``; this apparent published inconsistency is preserved as a note, not corrected.

**Reference.** [DeLima2018](@cite).  A. C. S. de Lima et al., “A Noniterative Approximation of a Full-Wave Model of Thin Wire Above and Buried in a Lossy Ground,” *IEEE Transactions on Electromagnetic Compatibility*, 60(6), 1873–1881 (2018), DOI `10.1109/TEMC.2017.2762241`, equations (1)–(16), pp. 1874–1875.

**Transcription source.** Original IEEE page images, printed pp. 1874–1875. The modal dependence, all ``\gamma/\gamma_1`` factors, separate integrals, exponential signs, and qFW barred substitutions were visually checked. The equation is preserved even where the text/OCR and printed symbols are awkward.

## Source transcription

The full-wave parent modal equation is

```math
M=\frac{2\pi}{j\omega\mu}z_i+
\left(1-\frac{\gamma^2}{\gamma_1^2}\right)\Lambda+
\left(S_1-\frac{\gamma^2}{\gamma_1^2}S_2\right)=0,
\qquad\text{(7)}
```

and the source iterates it with ``\gamma^{(n+1)}=\gamma^{(n)}-M/M'`` in (8). The qFW contribution is not a new spectral kernel; it is the noniterative prescription (15)–(16) applied to the parent.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\gamma`` | unchanged | unknown longitudinal modal propagation constant | ``\mathrm m^{-1}``, ``e^{-\gamma z}`` |
| ``\bar\gamma`` | unchanged | image-derived qFW prescribed longitudinal constant | ``\mathrm m^{-1}`` |
| ``\gamma_i`` | unchanged | bulk constant of medium ``i`` | ``\mathrm m^{-1}`` |
| ``u_i,\eta_1`` | unchanged | transverse/vertical spectral roots | ``\mathrm m^{-1}`` |
| ``\Lambda,S_1,S_2,S_4`` | unchanged | Bessel direct/image term and separate Sommerfeld integrals | source normalized |
| ``z_i`` | unchanged | conductor internal impedance | ``\Omega/\mathrm m`` |
| ``Z_c,Z`` | unchanged | characteristic impedance and per-unit-length impedance | ``\Omega`` and ``\Omega/\mathrm m`` |

No notation was renamed.

## Evidence and approximation sources

Equations (1)–(7) give the full-wave parent; (8)–(9) give the iterative solution; (12) gives the impedance; and §III plus (15)–(16) isolate the qFW replacement. Numerical comparisons in §§III–IV are validation evidence only.

## Limitations and discrepancies

-
- The image approximation needed to calculate ``\bar\gamma`` is cited but not reprinted, so this record is equation-verified but dependency-incomplete.
- The derivative-variable conflict in (9) requires an original erratum or author clarification; no repair is supplied.

