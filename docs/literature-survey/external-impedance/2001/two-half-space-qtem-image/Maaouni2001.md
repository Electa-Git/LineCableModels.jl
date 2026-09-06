# Maaouni–Amri two-half-space qTEM impedance approximation

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Self and mutual distances ``X,Y``. Thin bare wire; no finite coating. |
| Calculated quantities | Thin-wire longitudinal impedance near a lossy interface and qTEM logarithmic image approximation |
| Earth structure | Two homogeneous half-spaces. |
| Model and approximation | Equation (12) and the closed ``G`` form are qTEM/asymptotic reductions of the preceding exact spectral representation. |
| Main source | A. Maaouni, A. Amri, and N. Zouhir (2001) |
| Citation key(s) | `:Maaouni2001` |
| Evidence status | Original publication page images checked |

**Description.** Exact spectral parent plus analytical qTEM images for a thin wire in either homogeneous half-space adjoining a plane interface.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | qTEM reduction applied to the exact spectral parent. | §3. |
| Air propagation constant ``γ_air`` | ``k_0`` retained. | Definitions. |
| Earth propagation constant ``γ_earth`` | Ratio ``n`` represents the second medium. | (2)–(7). |
| Earth permittivity and displacement current | Complex material constants retained. | Constitutive definitions. |
| Range of validity | Thin infinite wire parallel to the interface. | Geometry. |
| Earth permeability ``μ_earth`` | Source material ratio retained. | Definitions. |
| Arrangement | Self and mutual distances ``X,Y``. | (2)–(12). |
| Earth structure | Two homogeneous half-spaces. | Figure 1. |
| Conductor and insulation geometry | Thin bare wire; no finite coating. | Model statement. |
| Constitutive and field assumptions | Linear isotropic media and harmonic fields. | Derivation. |
| Conventions | Complex logarithm/root branches follow source convention. | (12), (26). |

**Expression.** The qTEM impedance image integral reduces as

```math
J(X,Y)\simeq\ln\frac{\rho_J^*}{\rho^*},\qquad
\rho_J^*=\sqrt{X^2+(Y+Y_J)^2},\quad
\rho^*=\sqrt{X^2+Y^2},\quad
Y_J=\frac{2}{k_0\sqrt{1-n^2}}.
\qquad\text{(12)}
```

The full closed evaluator for the companion kernel is

```math
G(X,Y)\simeq\frac{n^2}{2(n^4-1)}[Q(bz)+Q(b\bar z)]
-\frac{P(b,z)+P(b,\bar z)-P(-b,z)-P(-b,\bar z)-n^2b[Q(-bz)+Q(-b\bar z)]}{2b(n^4-1)},
\qquad\text{(26)}
```

where ``z=k_0(Y+jX)``, ``b=j/\sqrt{1+n^2}``, ``Q(z)=e^{-z}E_1(-z)``, and ``P`` is defined by (25).

**Approximation.** Equation (12) and the closed ``G`` form are qTEM/asymptotic reductions of the preceding exact spectral representation.

**Limitations.** Two half-spaces and thin-wire geometry only.

**Reference.** [Maaouni2001](@cite), equations (2)–(26).

**Transcription source.** Original page images; conjugates, signs, and denominator ``n^4-1`` verified.

## Source transcription

The source's exact parent and its qTEM image are retained as one formulation with explicit approximation ancestry.

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
