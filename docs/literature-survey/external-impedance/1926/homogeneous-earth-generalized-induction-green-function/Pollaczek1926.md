# Pollaczek generalized induction coefficients for overhead, buried, and mixed conductors

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinitely long straight parallel filamentary conductors. A physical round wire of radius ``\rho`` is treated only through the stated small-radius self averaging; a buried wire is galvanically and dielectrically isolated from earth. No finite insulation impedance/admittance is supplied. |
| Calculated quantities | Complex generalized mutual-induction coefficient for air–air, earth–earth, and both mixed source/observation placements; source-prescribed finite-radius self coefficient and series-impedance assembly |
| Earth structure | Plane homogeneous conducting half-space (medium 1) below homogeneous air (medium 2). |
| Model and approximation | Not an analytical approximation within the explicitly reduced filamentary, two-dimensional physical model. The reductions ``\Gamma=0``, ``\epsilon_1=0``, ``\mu_1=\mu_2=1`` and ``k_2=0`` precede application. The finite-radius self rule is an additional small-radius approximation; separate small- and large-``\|k\eta\|`` asymptotic self expressions (59a)–(59b) are not substituted here. |
| Main source | F. Pollaczek (1926) |
| Citation key(s) | `:Pollaczek1926` |
| Evidence status | Original-publication scan image verified for the displayed equations and definitions |

**Description.** Source-defined complex generalized induction coefficients per unit length between infinitely long parallel filamentary conductors above a homogeneous conducting earth, below its plane surface, or on opposite sides of the interface. The longitudinal electric Green function is converted by the author to mutual induction; the corresponding ``i\omega M`` term is the external series contribution under the source convention.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fixed zero: axial attenuation/variation is neglected and ``\partial/\partial z=0``. | Stated — reduction (b) and equation (2), printed p. 341. |
| Air propagation constant ``γ_air`` | Formal source symbol ``k_2`` is retained while constructing the Green functions, then set to zero in (22) when air conductivity and permittivity are neglected. | Stated — reduction (c), p. 341, and (22), p. 345. |
| Earth propagation constant ``γ_earth`` | Source symbol ``k_1``, renamed ``k`` after (22). Under the selected reductions ``k=e^{-i\pi/4}\sqrt{4\pi\omega\sigma_1}/c`` in source cgs notation; the root has positive real part and negative imaginary part. | Stated — (1a), (10), pp. 340–342, and (22), p. 345. |
| Earth permittivity and displacement current | Neglected: the earth dielectric constant is set aside in reduction (c). Conductivity remains. | Stated — p. 341. |
| Range of validity | Infinite filamentary source; planar homogeneous half-spaces; steady harmonic state; no longitudinal attenuation. For finite-radius self use, radius ``\rho`` must be very small relative to conductor height/depth, and earth-induced asymmetry of conductor current is neglected. | Stated — pp. 339–341 and finite-wire discussion on p. 354. |
| Earth permeability ``μ_earth`` | Relative permeability fixed to one throughout the physical reduction. | Stated — reduction (c), p. 341. |
| Arrangement | Source and observation conductors can both be overhead, both buried, or mixed. Reciprocity is stated in (9). Self is obtained by the source's boundary-circle averaging prescription. | Stated — (9), pp. 342 and 344–345; §8, pp. 354–355. |
| Earth structure | Plane homogeneous conducting half-space (medium 1) below homogeneous air (medium 2). | Stated — Fig. 1 and pp. 340–341. |
| Conductor and insulation geometry | Infinitely long straight parallel filamentary conductors. A physical round wire of radius ``\rho`` is treated only through the stated small-radius self averaging; a buried wire is galvanically and dielectrically isolated from earth. No finite insulation impedance/admittance is supplied. | Stated — reductions (a) and text after (19), pp. 341 and 344. |
| Constitutive and field assumptions | Homogeneous linear isotropic media; time-harmonic Maxwell equations reduced to a two-dimensional axial-electric-field problem; air constitutive terms and earth displacement current omitted in the final specialization. Finite-wire internal field is delegated to the standard solid-conductor skin-effect solution. | Stated — pp. 340–341 and §8, p. 354. |
| Conventions | Common factor ``e^{i\omega t}``; ``z`` along the conductors; ``y`` positive upward with earth ``y<0``; source point ``(\xi,\eta)`` and observation point ``(x,y)``. Original formulas use cgs/practical electromagnetic units. | Stated — coordinate and time convention, p. 340; unit conversion discussion below (58), p. 354. |

**Expression.** The longitudinal electric Green functions for all four source/observation placements, with the author's conversion to generalized induction coefficient, equations (19a), (19c), (20a), (20b), and (58).

Let

```math
\alpha_m(s)=\sqrt{s^2-k_m^2},\qquad
d=\sqrt{(x-\xi)^2+(y-\eta)^2},\qquad
D=\sqrt{(x-\xi)^2+(y+\eta)^2},
```

and ``A=2iJ\omega/c^2`` as in (12). For an overhead source ``\eta\ge0`` and earth observation ``y\le0``,

```math
\mathfrak E_{-+}
=\frac{2A}{i\pi}\int_{-\infty}^{\infty}
\frac{
e^{,is(x-\xi)+y\alpha_1(s)-\eta\alpha_2(s)}
}{\alpha_2(s)+\alpha_1(s)}\,ds.
\tag{19a}
```

For overhead source and observation, ``y\ge0,\eta\ge0``,

```math
\begin{aligned}
\mathfrak E_{++}={}&A\left[H_0^{(1)}(k_2d)-H_0^{(1)}(k_2D)\right]\\
&+\frac{2A}{i\pi}\int_{-\infty}^{\infty}
\frac{e^{,is(x-\xi)-(y+\eta)\alpha_2(s)}}
{\alpha_2(s)+\alpha_1(s)}\,ds.
\end{aligned}
\tag{19c}
```

For a buried source ``\eta\le0`` and air observation ``y\ge0``,

```math
\mathfrak E_{+-}
=\frac{2A}{i\pi}\int_{-\infty}^{\infty}
\frac{
e^{,is(x-\xi)-y\alpha_2(s)+\eta\alpha_1(s)}
}{\alpha_2(s)+\alpha_1(s)}\,ds.
\tag{20a}
```

For buried source and observation, ``y\le0,\eta\le0``,

```math
\begin{aligned}
\mathfrak E_{--}={}&A\left[H_0^{(1)}(k_1d)-H_0^{(1)}(k_1D)\right]\\
&+\frac{2A}{i\pi}\int_{-\infty}^{\infty}
\frac{e^{,is(x-\xi)+(y+\eta)\alpha_1(s)}}
{\alpha_2(s)+\alpha_1(s)}\,ds.
\end{aligned}
\tag{20b}
```

The physical specialization is

```math
k_1=k,\qquad k_2=0,\qquad \xi=0.
\tag{22}
```

For the applicable placement-specific field ``\mathfrak E``, the source defines

```math
\mathfrak E=-J\frac{\partial M}{\partial t}=-i\omega J M,
\qquad
M=\frac{i}{\omega J}\mathfrak E,
\tag{58}
```

and for a physical conductor defines self ``L`` by circumferentially averaging ``M``. The total per-length series impedance with earth return is

```math
i\omega L+r,
\tag{60}
```

where ``r`` is the conductor's per-length alternating-current impedance. Thus ``M`` and ``L`` are complex, frequency-dependent induction coefficients, not purely geometric inductances.

The square roots follow (14a)–(14c): ``\alpha_m(s)\sim s`` as real ``s\to+\infty``, ``\alpha_m(s)\sim-s`` as ``s\to-\infty``, and ``\alpha_m(0)=ik_m`` for the source's ``\operatorname{Im}k_m<0`` branch. ``H_0^{(1)}`` is the first-kind Hankel function.

**Approximation.** Not an analytical approximation within the explicitly reduced filamentary, two-dimensional physical model. The reductions ``\Gamma=0``, ``\epsilon_1=0``, ``\mu_1=\mu_2=1`` and ``k_2=0`` precede application. The finite-radius self rule is an additional small-radius approximation; separate small- and large-``|k\eta|`` asymptotic self expressions (59a)–(59b) are not substituted here.

**Limitations.** The record preserves the source's cgs field normalization and its own conversion relation rather than silently inserting a modern SI prefactor. It supplies no earth displacement current, arbitrary permeability, longitudinal propagation, finite insulation layer, or proximity-distorted finite-conductor boundary solution. The ``k_2\to0`` Hankel differences are limiting expressions and must not be evaluated as two independent singular terms. The exact spectral representation does not by itself establish that later modern formulas attributed to Pollaczek use identical signs or normalizations.

**Reference.** [Pollaczek1926](@cite), equations (19a), (19c), (20a), (20b), (22), and (58)–(60), printed pp. 344–345 and 354–355 (PDF pages 6–7 and 16–17), with assumptions on pp. 340–342.

**Transcription source.** Original 1926 scan. The four signed placement exponents, full real-line measures, Hankel direct-minus-image terms, denominator ordering, root continuation, time convention, and field-to-induction conversion were checked visually on the printed pages. OCR text was used only to locate the sections; it was not used to resolve mathematical tokens.

## Source transcription

The original uses ``\mathfrak E`` for the only nonzero electric-field component after the two-dimensional reduction and indexes the fields by observation medium first and source medium second. Medium 1 is earth and medium 2 is air. Equations (19a), (19c), (20a), and (20b) above preserve that order.

For finite-radius self, the source says to average ``M`` around the conductor boundary. When radius ``\rho`` is small relative to height ``\eta``, this reduces to replacing the singular logarithmic factor

```math
\log\frac{2}{\sqrt{|k|^2\left[x^2+(y-\eta)^2\right]}}
```

by ``\log\!left(2/(\sqrt{|k|^2}\rho)\right)`` and then setting ``x=0,y=\eta`` in the remaining terms. This is a source prescription, not a finite-radius rederivation of the spectral kernel.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\mathfrak E`` | unchanged | Longitudinal electric field induced at ``(x,y)`` | original cgs field units; common time factor suppressed |
| ``J`` | unchanged | Source conductor current amplitude | source current convention |
| ``M,L`` | unchanged | Generalized mutual and self induction coefficient per unit length | source gives H/cm after stated practical conversion |
| ``r`` in (60) | unchanged | Conductor AC impedance per unit length | distinct from distance ``d`` used here |
| ``x,y;\xi,\eta`` | unchanged | Observation and source transverse coordinates | earth ``y,\eta<0``; air ``y,\eta>0`` |
| ``d,D`` | introduced display abbreviations | Direct and reflected-point distances appearing explicitly under the source radicals | one-to-one definitions shown above |
| ``k_1,k_2`` | unchanged | Earth and air medium constants | ``\mathrm{length}^{-1}``; source root convention |
| ``\alpha_m`` | introduced display abbreviation | ``\sqrt{s^2-k_m^2}`` | no coordinate or branch change |
| ``s`` | unchanged | Fourier integration variable | full real line |
| ``A`` | unchanged | ``2iJ\omega/c^2`` | source cgs normalization |
| ``H_0^{(1)}`` | unchanged | First-kind order-zero Hankel function | source branch |
| ``i,\omega,c`` | unchanged | Imaginary unit, angular frequency, speed of light | common factor ``e^{i\omega t}`` |

Only the explicitly defined ``d,D,\alpha_m`` abbreviations are added; their maps are exact and reversible.

## Evidence and approximation sources

- Source geometry, time convention and homogeneous media: printed pp. 339–340.
- Filament, zero-longitudinal-variation, constitutive and permeability reductions: printed p. 341.
- Medium-root definition and branch: (10), (14a)–(14c), printed pp. 342–343.
- Exact Green functions and reciprocity: (9), (19a), (19c), (20a), (20b), printed pp. 342 and 344–345.
- Final air/earth specialization: (22), printed p. 345.
- Generalized mutual coefficient, finite-radius self prescription, and series assembly: §8, (58)–(60), printed pp. 354–355.

## Limitations and discrepancies

- The scan's OCR loses radicals, signs and subscripts in the central equations. Every displayed token above follows the page image.
- The source uses a lossy-medium root with negative imaginary part under ``e^{i\omega t}``. Later ``e^{j\omega t}`` transcriptions often choose a positive-imaginary diffusion constant; no silent conjugation is made here.
- The real part of ``i\omega L`` is identified by the author with the earth-return ohmic resistance, while ``r`` contains the conductor contribution. The record does not relabel the Green function alone as a complete terminal impedance.
- The source contains numerous asymptotic series and self limits beyond this exact representation. Those require separate records and are not claimed complete here.
