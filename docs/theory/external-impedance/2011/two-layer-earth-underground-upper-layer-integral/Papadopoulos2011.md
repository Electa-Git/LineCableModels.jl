# Papadopoulos–Tsiamitros–Papagiannis upper-layer buried-cable earth impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinite parallel single-core cable axes at positive burial depths ``h_i,h_j`` with horizontal separation ``y_{ij}``. Finite radius enters the self substitution. Internal conductor and insulation contributions are separate. |
| Calculated quantities | Per-unit-length mutual earth-return impedance; self earth term by the published radius/depth substitutions |
| Earth structure | Air above a horizontal finite first earth layer of thickness ``d``; second earth layer is the terminal infinite-depth half-space. This is not an arbitrary-layer or different-layer-source kernel. |
| Model and approximation | The source replaces the unknown longitudinal propagation constant by the dielectric value of the upper earth layer and applies identity (5). Equation (6) remains a spectral integral; no analytical truncation is stated. |
| Main source | T. A. Papadopoulos, D. A. Tsiamitros, and G. K. Papagiannis, 2011, DOI `10.1049/iet-gtd.2010.0228` |
| Citation key(s) | `:Papadopoulos2011` |
| Evidence status | Original PDF equations checked; the transformed-factor index and root branch remain unresolved. |

**Description.** Per-unit-length earth-return impedance between infinite parallel single-core cable axes buried in the finite upper layer of a two-layer earth below air. The mutual spectral kernel contains direct and interface-dependent terms; the source obtains the self term by substitution.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Source ``\gamma_x`` enters ``e^{-\gamma_x x}``. The final expression prescribes ``\gamma_x=j\omega\sqrt{\mu_1\varepsilon_1}``, purely imaginary, rather than retaining the lossy bulk constant or setting it to zero. The actual line propagation constant is subsequently calculated from the resulting per-unit-length parameters. | Stated — section 2, printed p. 162, right column before (5); (1). |
| Air propagation constant ``γ_air`` | Appendix indexed definition is ``\gamma_0^2=j\omega\mu_0(\sigma_0+j\omega\varepsilon_0)``. Air is labelled by ``\varepsilon_0,\mu_0`` in Fig. 1; an explicit value for ``\sigma_0`` is not supplied in the inspected definition. The air spectral factor is retained through ``S_{10},D_{10}``, not set to zero. | Stated/equation-implied — Fig. 1, (6b), appendix definition and (19)–(20), printed pp. 162–163, 171. |
| Earth propagation constant ``γ_earth`` | ``\gamma_k^2=j\omega\mu_k(\sigma_k+j\omega\varepsilon_k)``, ``k=1,2``. These bulk constants differ from the prescribed longitudinal ``\gamma_x``. The transformed-factor sentence prints ``a_m=\sqrt{\lambda^2+\gamma_k^2+k_x^2}``; its ``m/k`` index mismatch is retained, not repaired. | Stated — appendix below (10), p. 171, and transformation below (5), p. 162. |
| Earth permittivity and displacement current | Each earth layer retains ``\varepsilon_k`` in its bulk propagation constant. Only the imposed longitudinal constant discards the conductive part. No separate dielectric-loss or complex-permittivity decomposition is specified. | Equation-implied — bulk definition, p. 171; stated longitudinal approximation, p. 162. |
| Range of validity | Quasi-TEM transmission-line mode only; antenna modes ignored. Authors cite satisfactory TL approximation up to 10 MHz for typical earth resistivity, not a universal bound for every geometry/material. Both cable axes must lie in the finite first earth layer. Numerical comparisons in sections 4–5 are tests, not proof of that bound or a truncation-error estimate. | Stated — section 3, p. 163; Fig. 1 and section 2, p. 162. |
| Earth permeability ``μ_earth`` | Distinct scalar ``\mu_1,\mu_2`` retained. Relative permeabilities equal to one in section 4 are numerical-example choices, not restrictions built into (6). | Equation-implied — (6), (19)–(20); stated example restriction — section 4.1, p. 164. |
| Arrangement | Underground mutual interaction with both cables in layer 1; self by replacing ``y_{ij}`` with the outermost cable radius and ``h_j`` with ``h_i``. Both depth orderings are stated to lead to the same final expressions. | Stated — section 2, p. 163, and final appendix paragraph, p. 171. |
| Earth structure | Air above a horizontal finite first earth layer of thickness ``d``; second earth layer is the terminal infinite-depth half-space. This is not an arbitrary-layer or different-layer-source kernel. | Stated — section 2 and Fig. 1, p. 162; proposed future extensions, section 3, p. 163. |
| Conductor and insulation geometry | Infinite parallel single-core cable axes at positive burial depths ``h_i,h_j`` with horizontal separation ``y_{ij}``. Finite radius enters the self substitution. Internal conductor and insulation contributions are separate. | Stated — section 2 and Fig. 1, pp. 162–163. |
| Constitutive and field assumptions | Each layer is represented by spatially uniform scalar ``\mu,\varepsilon,\sigma``; linear superposition and infinite longitudinal geometry underlie the Hertz-vector construction. Quasi-TEM, with an imposed dielectric longitudinal approximation; no new skin/proximity formula is supplied by this earth kernel. | Equation-implied — indexed scalar definitions and (1); stated — sections 2–3, pp. 162–163. |
| Conventions | ``j`` imaginary unit, ``\omega=2\pi f``; longitudinal factor ``e^{-\gamma_x x}``. The constitutive factors are consistent with positive ``j\omega`` time differentiation; no explicit time exponential is supplied in the inspected derivation. Depths are downward from the air interface, whereas appendix ``z`` uses distances ``d-h_i``. Output is per-unit-length earth impedance, not total cable impedance. | Stated/equation-implied — (1), Fig. 1, (6a), appendix (8)–(10) and following definitions. |

**Expression.** Mutual earth-return impedance, (6a)–(6b), printed pp. 162–163. Source notation is retained, including its displayed ``h_1,h_2`` instead of the prose's ``h_i,h_j``.

```math
Z'_{e_{ij}}=\frac{j\omega\mu_1}{2\pi}
\int_0^{+\infty}F(\lambda)\cos(y_{ij}\lambda)\,d\lambda.
\qquad\text{(6a)}
```

```math
F(\lambda)=
\frac{
 S_{10}S_{21}e^{-a_1|h_1-h_2|}
 +S_{10}D_{21}e^{-a_1(d-h_1+d-h_2)}
 -D_{10}S_{21}e^{-a_1(h_1+h_2)}
 -D_{10}D_{21}e^{-a_1(2d-|h_1-h_2|)}
}{a_1(S_{10}S_{21}+D_{10}D_{21}e^{-2a_1d})}.
\qquad\text{(6b)}
```

The author-prescribed transformation is ``u^2-k_x^2=\lambda^2``, with the following printed factor definition and the unprimed versions of appendix (19)–(20). The ``k`` on ``\gamma_k`` in the definition of ``a_m`` is intentionally retained; the appendix separately defines the consistently indexed primed factor below. Equating ``k`` and ``m`` would be an interpretation of the printed mismatch, not a source transcription.

```math
\begin{aligned}
a_m&=\sqrt{\lambda^2+\gamma_k^2+k_x^2} \\
\gamma_k^2&=j\omega\mu_k(\sigma_k+j\omega\varepsilon_k) \\
k&=0,1,2,
\end{aligned}
```

```math
\begin{aligned}
S_{mn}&=(\mu_n a_m+\mu_m a_n) \\
D_{mn}&=(\mu_m a_n-\mu_n a_m),\qquad m,n=0,1,2,
\end{aligned}
```

```math
\begin{aligned}
\gamma_x&=j\omega\sqrt{\mu_1\varepsilon_1} \\
\gamma_x&=jk_x \\
\omega&=2\pi f.
\end{aligned}
```

The equality ``\gamma_x=jk_x`` identifies the exponent of (5) with the preceding prescribed exponent; it is not an independent propagation model. ``Z'`` is in impedance per length; ``h,d,y`` are lengths, ``\lambda,a,\gamma,k_x`` inverse lengths. ``h_1,h_2`` denote the two cable depths in (6b), not the two material-layer thicknesses; the correspondence to ``h_i,h_j`` is contextual and is not a coordinate-sign change. For the self term the source explicitly replaces ``y_{ij}`` by cable ``i``'s outermost radius and ``h_j`` by ``h_i``. No additional free-space or internal term is inserted here.

**Approximation.** The parent is the Hertz-vector longitudinal integral (1), with amplitude (3), before prescribing the unknown ``\gamma_x``. The author's operation is to approximate longitudinal propagation by a purely imaginary dielectric-medium value and apply Bessel/Fourier identity (5). No small-parameter expansion order, discarded integral tail, or analytical truncation is stated. Within that prescribed field model (6) is an integral representation, not a claimed exact solution of the unknown full-wave dispersion relation. A separate ``\gamma_x=0`` reduction is attributed to earlier work, not substituted into this record.

**Limitations.** Both conductors occupy the finite upper earth layer. The source's generalization discussion does not supply a formula for conductors in different layers or a mixed overhead/buried pair. The transformed-factor indices and square-root branch remain unresolved. The equation is preserved as printed; neither appendix inconsistencies nor LCM code are used to repair it.

**Reference.** [Papadopoulos2011](@cite), section 2, (1), (3), (5), (6a)–(6b), printed pp. 162–163; section 3, p. 163; appendix definitions and (19)–(20), p. 171.

**Transcription source.** Original publication; equations and dependent coefficients visually checked on PDF pages 2, 3, and 11 (indices 1, 2, 10), including higher-resolution inspection of the printed ``a_m/\gamma_k`` mismatch. No Markdown conversion was treated as evidence. Image verification confirms the printed witness, not resolution of its inconsistent indices.

## Source transcription

The formula section reproduces (6a)–(6b) in source order. The following auxiliary identity and original primed definitions preserve the operation used to obtain that kernel.

Equation (5), printed p. 162:

```math
\int_{-\infty}^{\infty}
 J_0\!\left(u\sqrt{x^2+y_{ij}^2}\right)e^{-jk_xx}\,dx
=\begin{cases}
0,&u<k_x,\\
2\dfrac{\cos\!\left(y_{ij}\sqrt{u^2-k_x^2}\right)}
{\sqrt{u^2-k_x^2}},&u>k_x.
\end{cases}
\qquad\text{(5)}
```

The source states no separate value at ``u=k_x``. Immediately following (5), it makes ``u^2-k_x^2=\lambda^2``, prints the ``a_m`` definition reproduced above, and states that primed ``S,D,\Delta,A`` become their unprimed counterparts. ``J_0`` is the Bessel function of the first kind, order zero (appendix below (10)).

Appendix unnumbered definition below (10), then (19)–(20), printed p. 171:

```math
\begin{aligned}
a'_k&=\sqrt{u^2+\gamma_k^2} \\
\gamma_k^2&=j\omega\mu_k(\sigma_k+j\omega\varepsilon_k) \\
k&=0,1,2.
\end{aligned}
```

```math
S'_{mn}=(\mu_n a'_m+\mu_m a'_n),\qquad\text{(19)}
```

```math
D'_{mn}=(\mu_m a'_n-\mu_n a'_m).\qquad\text{(20)}
```

Section 3, p. 163, states three reductions without printing new independent kernels: first-earth-layer properties equal to air and ``k_x=k_0=\omega\sqrt{\mu_0\varepsilon_0}`` give the overhead homogeneous-earth expressions attributed to Kikuchi [18]; equal properties of earth layers 1 and 2 give the authors' homogeneous-earth formulas [9]; a further sentence sets ``\gamma_x=0`` and identifies the earlier two-layer impedance [13]. The latter sentence's placement after the equal-layer case is retained as prose, not used to impose simultaneous incompatible layer conditions. The homogeneous-earth original has a [separate record](../../2010/homogeneous-earth-underground-cable-correction/Papadopoulos2010b.md).

The accompanying potential-coefficient/admittance expression is in the [separate admittance record](../../../external-admittance/2011/two-layer-earth-underground-upper-layer-potential-coefficient/Papadopoulos2011.md). Internal and insulation assembly is referenced to [9,15] by the paper; no new full cable matrix is printed here.

## Notation map

No formula section renaming is used. Source primes and the Latin ``a`` are retained. Units below are SI meanings; where not printed next to the equation, they follow the displayed definition rather than constituting an additional source restriction.

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z'_{e_{ij}}`` | unchanged | Mutual per-length earth impedance | ``\Omega/\mathrm{m}``; prime denotes per-unit-length output |
| ``F(\lambda)`` | unchanged | Four-term spectral earth kernel | Source (6b) |
| ``i,j`` | unchanged | Cable source/target labels | Not medium indices |
| ``h_i,h_j``; displayed ``h_1,h_2`` | unchanged | Depths of the two cable axes below the earth surface | Length; source changes subscript notation in (6b) |
| ``d`` | unchanged | Thickness of first earth layer | Positive length |
| ``y_{ij}`` | unchanged | Horizontal axis separation | Length; outer radius replaces it for self |
| ``x,z`` | unchanged | Longitudinal coordinate; appendix vertical coordinate | Length; appendix compares ``z`` with ``d-h_i`` |
| ``u,\lambda`` | unchanged | Original and transformed spectral variables | Inverse length; ``u^2-k_x^2=\lambda^2`` |
| ``a'_k,a_m`` | unchanged | Original/transformed vertical factors | Inverse length; transformed ``m/k`` mismatch and branch unresolved |
| ``\gamma_k`` | unchanged | Bulk propagation constant, air 0 and earth 1,2 | Inverse length; not the imposed longitudinal constant |
| ``\gamma_x,k_x,k_0`` | unchanged | Imposed longitudinal constant; corresponding wavenumber; source's free-space reduction wavenumber | Inverse length; ``e^{-\gamma_x x}`` |
| ``\mu_k,\varepsilon_k,\sigma_k`` | unchanged | Indexed permeability, permittivity, conductivity | ``\mathrm{H/m},\mathrm{F/m},\mathrm{S/m}`` |
| ``S'_{mn},D'_{mn};S_{mn},D_{mn}`` | unchanged | Primed and transformed interface coefficients | Defined by (19)–(20) and the stated transform; ``D`` is not ``\Delta`` |
| ``m,n,k`` | unchanged | Material indices | 0,1,2; mismatch in transformed ``a_m`` retained |
| ``J_0`` | unchanged | First-kind Bessel function, order zero | Dimensionless function of dimensionless argument |
| ``\omega,f,j`` | unchanged | Angular frequency, frequency, imaginary unit | ``\omega=2\pi f``; ``j^2=-1`` |

## Evidence and approximation sources

- Original contribution: two-layer earth-return impedance with the stated dielectric longitudinal prescription; the paper distinguishes earlier zero-longitudinal-propagation impedance [13] from its (6). The admittance contribution is separately recorded.
- Geometry and field parent: Fig. 1 and (1), section 2, p. 162. The finite earth layer is between air and an infinite-depth second earth layer.
- Actual approximation: section 2, right column before (5); the source first discusses unknown lossy propagation, then imposes a purely imaginary approximation motivated by the high-frequency dielectric contribution. It supplies no quantitative error bound for that replacement.
- Transformation: (5) and following sentence; the consistently indexed primed factor in the appendix is a distinct printed witness to the transformed-factor index problem, not permission to repair it.
- Exact coefficients within the printed representation: appendix (19)–(20), with primes removed by the expressly stated transform. Intermediate Hertz coefficients are derivational dependencies, not extra impedance formulas.
- Self and depth-order coverage: p. 163 after (7), and final paragraph of the appendix, p. 171. The latter explicitly says either cable depth ordering gives the same (6) and (7).
- Accuracy evidence is comparison with earlier approximate and homogeneous-earth models in sections 4–5, not an independently established universal error measure. The 10 MHz TL statement is author-attributed to reference [19].

## Limitations and discrepancies

- **Suspected published indexing defect:** immediately after (5), the original prints ``a_m=\sqrt{\lambda^2+\gamma_k^2+k_x^2}``, although the earlier primed definition uses matching ``k``. Both witnesses are retained. No ``k=m`` correction is installed.
- **Not stated:** a branch/radiation prescription for the square roots in the inspected formulas, and an explicit value for air ``\sigma_0`` in the appendix indexed definition. Neither is invented.
- **Published notation change:** the prose/Fig. 1 uses ``h_i,h_j`` while final kernels use ``h_1,h_2``; appendix coefficients (17)–(18) further use an unsubscripted ``h``. The final kernel is not silently rewritten to unify these.
- **Unreconciled parent-coefficient discrepancy:** appendix (17) prints denominator ``S_{10}S_{21}e^{a_1d}+S_{21}D_{10}e^{-a_1d}``, whereas (18) uses ``S_{10}S_{21}e^{a_1d}+D_{21}D_{10}e^{-a_1d}``. The final (6b) is preserved independently; this record does not assert the parent coefficients algebraically establish it without discrepancy.
- **Scope:** different-layer buried conductors, mixed arrangements, and additional earth layers are suggested extensions on p. 163, not extracted coverage from this source.
- **Conversion status:** no corresponding Markdown was found, so there is no conversion-to-PDF comparison to report.

## Numerical interpretation

Numerical evaluation retains the source's upper-earth-layer restriction and physical top-layer thickness. The two printed reflection distances remain distinct for unequal cable depths. Independent evaluations of the magnetic kernel check permeability contrast, reciprocity, and negative-frequency conjugacy.
