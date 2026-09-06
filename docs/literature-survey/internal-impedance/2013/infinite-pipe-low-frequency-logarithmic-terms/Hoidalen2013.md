# Høidalen infinite-pipe low-frequency logarithmic terms

## Identity and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Infinite-wall parent has finite inner radius ``r_{p1}``, no outer radius in (14)–(16). Core offsets ``d_k,d_m`` and angle ``\theta_{km}`` define coupling. The separately labelled mixed-model diagnostic reuses finite-wall outer/transfer terms and pipe DC resistance. |
| Calculated quantities | Low-frequency logarithmic inner-pipe surface term; pipe-mediated core-coupling correction; author's incomplete diagnostic expression for mixed finite/infinite pipe assembly |
| Earth structure | Not applicable. |
| Model and approximation | Høidalen takes the low-frequency behavior of the infinite-wall Bessel parent (2)–(3), identifying its ``K_0`` logarithmic behavior and retaining the additional ``C_1`` term (16). The source does not supply a complete term-by-term expansion, remainder, or a uniform finite-thickness limit. Taking a low-frequency limit after an infinite-wall assumption is not equated here with taking the low-frequency limit at fixed physical wall thickness. |
| Main source | Høidalen's low-frequency analysis of the infinite-wall model attributed to Brown–Rocamora (1976) and Ametani (1980) |
| Citation key(s) | `:Hoidalen2013` |
| Evidence status | Manuscript page images checked for (14)–(16) and their dependencies; the source leaves diagnostic ``X`` unexpanded. |

**Description.** Low-frequency logarithmic expressions for the inner-surface impedance and pipe-mediated core interaction in an infinite-wall cylindrical pipe model. The source compares these with finite-wall results and diagnoses a different, mixed finite/infinite assembly; that diagnostic is not presented as a usable corrected model.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not stated; the source's pipe skin parameter is not an impressed longitudinal constant. No axial exponential or ``Γ=0`` prescription accompanies these formulas. | Equation-implied — (2),(3),(7),(14)–(16); introduction p. 1. |
| Air propagation constant ``γ_air`` | Not applicable to the pipe-surface/core-coupling formulation. | Equation-implied — (14)–(16). |
| Earth propagation constant ``γ_earth`` | Not applicable; no earth term enters the primary expressions or the pipe-only diagnostic (23). | Equation-implied — (14)–(16),(23). |
| Earth permittivity and displacement current | Not applicable. | Equation-implied — pipe skin definition (7) and the primary outputs. |
| Range of validity | Source evaluates ``\omega\to0`` in an already infinite-wall model. It explicitly warns that this model is inadequate for a physically finite pipe at low frequency. No quantitative small-``x_1`` threshold or error remainder is stated. This is not a recommendation to use an infinite pipe at low frequency. | Stated — pp. 1,3; (14)–(16) and discussion around (23). |
| Earth permeability ``μ_earth`` | Not applicable; ``\mu_{rp}`` is relative pipe permeability. | Stated/equation-implied — Table 1 and (7). |
| Arrangement | Core conductors inside a common cylindrical pipe; self/mutual pipe-mediated coupling depends on their radial offsets and angle. No overhead/buried earth-return arrangement is implied. | Stated — Fig. 1 and (1),(3),(5). |
| Earth structure | Not applicable. | Equation-implied — no earth medium in these outputs. |
| Conductor and insulation geometry | Infinite-wall parent has finite inner radius ``r_{p1}``, no outer radius in (14)–(16). Core offsets ``d_k,d_m`` and angle ``\theta_{km}`` define coupling. The separately labelled mixed-model diagnostic reuses finite-wall outer/transfer terms and pipe DC resistance. | Stated/equation-implied — (2),(3),(5),(14)–(16),(23). |
| Constitutive and field assumptions | Uniform scalar pipe conductivity ``\sigma_p`` and permeability ``\mu_0\mu_{rp}``; conductive skin argument omits displacement current. These are pipe terms, not a complete core-to-core proximity treatment. | Equation-implied — (7); section IV.C is separate. |
| Conventions | ``i`` superscript means infinite pipe thickness, not imaginary unit; ``j`` is imaginary. Positive ``j\omega`` factors, natural logarithm and Euler constant ``\gamma``. The source retains ``\omega`` on the right of printed low-frequency ``\lim`` expressions. | Stated — p. 2 below (3), p. 3 below (16); (14)–(16). |

**Expression.** Manuscript equations (14)–(16), retaining all factors and the source's limit notation.

```math
\lim_{\omega\to0} Z_{pi}^{i}
=\frac{j\omega\cdot\mu_0\mu_{rp}}{2\pi}\cdot
(\ln(2/x_1)-\gamma).
\qquad\text{(14)}
```

```math
\lim_{\omega\to0} Z_{p\Sigma,km}^{i}=Z_{p\Sigma,km}+\Delta^{i}.
\qquad\text{(15)}
```

```math
\Delta^{i}=\frac{j\omega\cdot\mu_0}{2\pi}\cdot
\frac{2\mu_{rp}}{1+\mu_{rp}}\cdot C_1\cdot
\frac{(\ln(2/x_1)-\gamma)\cdot x_1^2}
{1+\mu_{rp}+(\ln(2/x_1)-\gamma)\cdot x_1^2}.
\qquad\text{(16)}
```

``\gamma`` is Euler's constant, not a propagation constant. The geometry/skin definitions from (5),(7) are:

```math
C_n=(d_k\cdot d_m/r_{p1}^2)^n\cdot\cos(n\theta_{km}),
\qquad
x_1=r_{p1}\cdot\sqrt{j\omega\cdot\mu_0\cdot\mu_{rp}\cdot\sigma_p}.
```

The source says ``C_1`` in (16) comes from (5), i.e. that definition's ``n=1`` term. It does not replace the parent's entire infinite sum by a new finite summation limit in (15).

The finite-model term referenced by (15) is the low-frequency expression immediately preceding it, (13):

```math
\lim_{\omega\to0} Z_{p\Sigma,km}
=\frac{j\omega\cdot\mu_0}{2\pi}\cdot\frac{2\mu_{rp}}{1+\mu_{rp}}\cdot
\ln\left(
\frac{r_{p1}^2}
{\sqrt{r_{p1}^4+(d_k\cdot d_m)^2-2r_{p1}^2\cdot d_k\cdot d_m\cdot\cos\theta_{km}}}
\right).
\qquad\text{(13)}
```

Equation (15) itself prints ``Z_{p\Sigma,km}`` on the right without a limit symbol; both witnesses are kept rather than rewriting the right-hand side. ``r_{p1},d_k,d_m`` are lengths, ``\theta_{km}`` an angle, ``x_1,C_n`` dimensionless, ``\sigma_p`` conductivity and ``\mu_{rp}`` relative pipe permeability. Outputs are impedance per length. Source branch prescriptions for the square root and complex logarithm are not stated.

**Approximation.** Høidalen takes the low-frequency behavior of the infinite-wall Bessel parent (2)–(3), identifying its ``K_0`` logarithmic behavior and retaining the additional ``C_1`` term (16). The source does not supply a complete term-by-term expansion, remainder, or a uniform finite-thickness limit. Taking a low-frequency limit after an infinite-wall assumption is not equated here with taking the low-frequency limit at fixed physical wall thickness.

**Limitations.** The author states that applying the infinite-wall model to finite pipes can produce erroneous low-frequency inductance. The separate combination of its inner-surface term with finite-wall outer/transfer terms can give negative low-frequency self resistance; the relevant complete author assembly is described below, not diagnosed from an isolated special function. The diagnostic ``X`` in (23) is explicitly unexpanded and cannot be treated as a complete evaluable expression.

**Reference.** [Hoidalen2013](@cite).  Høidalen (2013), DOI `10.1109/TPWRD.2013.2272343`; manuscript (2)–(3),(5),(7), p. 2; (13)–(16),(23), p. 3. The source attributes the infinite-pipe parent to Brown–Rocamora [1] and Ametani [2], not to this later analysis.

**Transcription source.** Original-author low-frequency analysis, checked against manuscript PDF images pp. 2–3. Parent Bessel expressions are secondary witnesses in Høidalen's notation. The manuscript has not been compared with the publisher's final edition, and the parent equations were not independently image-verified through this source.

## Source transcription

Source parent (2)–(3), manuscript p. 2:

```math
Z_{pi}\approx Z_{pi}^{i}
=\frac{j\omega\mu_0}{2\pi}\cdot\mu_{rp}
\frac{K_0(x_1)}{x_1\cdot K_1(x_1)},
\qquad\text{(2)}
```

```math
Z_{p\Sigma,km}\approx Z_{p\Sigma,km}^{i}
=\frac{j\omega\mu_0}{2\pi}\cdot
\sum_{n=1}^{\infty}
\frac{2\mu_{rp}\cdot C_n}
{n(1+\mu_{rp})+\dfrac{x_1\cdot K_{n-1}(x_1)}{K_n(x_1)}}.
\qquad\text{(3)}
```

The full source definitions of ``C_n,x_1`` are repeated in the formula section. ``K_n`` is not explicitly defined by kind in this manuscript; the cited finite-pipe original Kane et al., p. 1647 (HAL PDF page 3), expressly identifies ``K_n`` as the modified Bessel function of the second kind of order ``n``. That definition has been image-checked as reference-chain corroboration, not used to erase any difference between witnesses.

The low-frequency expressions occur in source order (13),(14),(15),(16), all reproduced in the formula section with their equation labels. The finite-pipe analysis has its own [record](../finite-pipe-low-frequency-surface-and-loop-terms/Hoidalen2013.md), which also preserves the finite-wall outer and transfer surface terms (20)–(21).

For the **different mixed-model assembly**, the author prints:

```math
\lim_{\omega\to0} Z_p^{i}
=\lim_{\omega\to0}(Z_{pi}^{i}+Z_{po}-2Z_{pm})
=-R_{p,dc}+\frac{j\omega\mu_0\cdot\mu_{rp}}{2\pi}\cdot X.
\qquad\text{(23)}
```

Here ``Z_{po},Z_{pm}`` are the finite-wall (20),(21), while only ``Z_{pi}^{i}`` uses the infinite-wall model. ``R_{p,dc}`` is the finite pipe's DC resistance per length. The author calls ``X`` an inductance expression obtainable from (14),(22) but does **not print it**. It therefore remains an unresolved dependency of (23), not a reconstructed formula. The reported negative-resistance problem concerns inserting this whole mixed pipe combination into the self cable impedance when pipe DC resistance exceeds core DC resistance. It is not a physical claim about the ``K_0`` function alone.

## Notation map

No renaming or algebraic replacement is used.

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{pi}^{i},Z_{p\Sigma,km}^{i}`` | unchanged | Infinite-wall inner-surface and pipe-mediated coupling terms | Per length; superscript ``i`` is geometric-model label |
| ``Z_{p\Sigma,km}`` | unchanged | Finite-model coupling term referenced by (15) | Source prints low-frequency form in (13), but no limit on (15)'s right-hand occurrence |
| ``\Delta^{i}`` | unchanged | Additional low-frequency logarithmic correction | Per length; contains ``C_1`` |
| ``Z_p^{i}`` | unchanged | Mixed infinite-inner/finite-outer-transfer diagnostic combination | Not a consistently infinite-wall complete pipe model |
| ``Z_{po},Z_{pm},R_{p,dc}`` | unchanged | Finite-wall outer/transfer impedances and DC pipe resistance in (23) | Per length; parent terms preserved in companion finite-pipe record |
| ``X`` | unchanged | Unprinted inductance factor in source diagnostic | Unresolved; no expression fabricated |
| ``\gamma`` | unchanged | Euler constant | Dimensionless; not bulk or longitudinal propagation |
| ``C_n,C_1`` | unchanged | Geometrical series coefficients | Dimensionless; (5) and its stated ``n=1`` use |
| ``r_{p1},d_k,d_m,\theta_{km}`` | unchanged | Pipe inner radius, core-axis offsets, their angular separation | Lengths and angle |
| ``x_1,\sigma_p,\mu_{rp},\mu_0`` | unchanged | Skin argument, pipe conductivity, relative/reference permeability | Dimensionless, ``\mathrm{S/m}``, dimensionless, ``\mathrm{H/m}`` |
| ``K_n,K_{n-1},K_0,K_1`` | unchanged | Source's Bessel functions and orders | Second-kind modified functions corroborated in the cited original; branch unspecified by the manuscript |
| ``k,m,n,\omega,j`` | unchanged | Core labels, positive-integer series order, angular frequency, imaginary unit | Parent ``n=1`` to infinity; positive ``j\omega`` factors |

## Evidence and approximation sources

- Infinite-wall parent and attribution: section IV.A, (2)–(3), p. 2; source expressly distinguishes this model from finite-wall (8)–(9). Numerical comparisons use a pipe-series upper limit of 20 after section V.A's convergence discussion, p. 4; that evaluation setting does not change the source's analytical infinity bound or determine the discarded low-frequency expansion order.
- Mathematical operation: logarithmic low-frequency behavior described before (14), followed by (14)–(16). The retained ``C_1`` correction is explicit; a full expansion/remainder is not.
- Finite versus infinite and mixed assembly are kept distinct. Høidalen's (23) does not justify replacing finite-wall (12) with (14) while preserving other finite terms and calling the result the same physical model.
- Brown and Rocamora's 1976 source was not inspected. Ametani's 1980 pipe section is represented only by the extracted concentric-insulation records. Høidalen's paper is therefore a secondary equation witness for these parents.
- Core-to-core proximity correction (36) is mathematically separate and has its own [record](../solid-core-proximity-low-frequency-subtracted-series/Hoidalen2013.md).

## Limitations and discrepancies

- **Author-stated physical failure:** the infinite-wall model is inappropriate at low frequency for a thin finite pipe. The mixed-model negative resistance discussed after (23) is source-attributed and depends on the full pipe/cable combination and resistance inequality.
- **Unresolved printed dependency:** ``X`` is not expanded in (23). This diagnostic expression is incomplete even though (14)–(16) include their required definitions.
- **Published limit notation:** source right-hand sides retain ``\omega``; (15)'s right side uses finite ``Z_{p\Sigma,km}`` without repeating the limit. No corrected asymptotic notation or limiting order is supplied by the reviewer.
- **No branch supplied:** the inspected manuscript does not specify square-root or logarithm branches. The source expressions are not rationalized, resummed, or regularized.
- **Conversion defects:** Markdown damages/interleaves equation text around the low-frequency formulas. The parent sum, Euler constant, ``x_1^2`` correction factors, and signs in (14)–(16),(23) were read from PDF images.
