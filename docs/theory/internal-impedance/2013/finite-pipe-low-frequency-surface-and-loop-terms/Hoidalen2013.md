# Høidalen finite-pipe low-frequency surface and loop terms

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Finite circular pipe and equal-core geometry using radii ``r_{p1},r_{p2},r_{p3},r_{1k},r_{2k}``, offsets ``d_k``, and separations ``d_{km}``. |
| Calculated quantities | Low-frequency pipe inner-surface impedance, pipe-mediated core interaction, surface-connection combination, and source-provided cable-loop/mode limits |
| Earth structure | Not applicable to pipe/internal and mode-cancelled terms. Not stated for the unspecified ground-return input in the total-loop formula. |
| Model and approximation | Low-frequency limits of the finite-pipe Bessel model combined through the printed cable and modal relations. Equation (35) also assumes a thin wall; the source gives no discarded order or error remainder. |
| Main source | Høidalen's low-frequency analysis of the finite-pipe formulas attributed to Kane et al. (1995), da Silva et al. (2006), and earlier cable assembly |
| Citation key(s) | `:Hoidalen2013` |
| Evidence status | Author manuscript equations checked; printed limit notation documented. |

**Description.** Low-frequency series-impedance terms for insulated solid cores enclosed by a conducting cylindrical pipe of finite wall thickness. The expressions retain the pipe's inner-surface resistance/inductance, pipe-mediated core coupling, the inner/outer-surface connection combination, and the author's symmetrical-core mode limits. Total cable-loop limits include a separately unspecified ground-return term.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not stated as an imposed parameter in these expressions. The introduction describes longitudinal propagation effects and neglect of transversal effects, but gives no ``Γ=0`` prescription or axial exponential. | Stated/equation-implied — manuscript p. 1 introduction; (7)–(24). |
| Air propagation constant ``γ_air`` | Not applicable to the extracted pipe/internal terms. The source does not define a bulk air term in its appended ground-return symbol either. | Equation-implied — (12),(13),(22),(24); separate ``Z_g`` paragraph p. 3. |
| Earth propagation constant ``γ_earth`` | Not applicable to the primary pipe and mode terms. For the total-loop formula (24), ``Z_g(0)`` is an unresolved external input; no earth propagation definition is supplied. | Stated — p. 3 after (18); equation-implied — (24), and cancellation discussion preceding (29),(32). |
| Earth permittivity and displacement current | Not applicable to the primary pipe/internal terms; not stated for the separately appended ``Z_g(0)`` in (24). | Equation-implied — (12),(13),(22),(24); source explicitly does not analyze ``Z_g``. |
| Range of validity | Author's low-frequency limit ``\omega\to0`` at fixed finite pipe geometry. No quantitative small-argument threshold, retained Bessel expansion order, or error bound is given. The further (35) also takes ``r_{p2}\to r_{p1}``; (24),(31),(34),(35) specify unity relative permeability. FEM's 1 Hz–1 MHz sweep is test coverage only. | Stated — (12)–(13),(24),(31),(34)–(35); FEM setup p. 2. |
| Earth permeability ``μ_earth`` | Not applicable to pipe/internal formulas; not stated for ``Z_g(0)``. Pipe permeability ``\mu_0\mu_{rp}`` must not be mistaken for earth permeability. | Equation-implied — pipe definition (7), (12),(13),(22). |
| Arrangement | Cores inside a common conducting pipe. Pipe-mediated self/mutual terms use source indices ``k,m``; total-loop self/mutual limits are given separately. Differential/common-mode limits are restricted to the symmetrical three-core case. | Stated — Fig. 1, (1),(17),(24), sections IV.D–E. |
| Earth structure | Not applicable to pipe/internal and mode-cancelled terms. Not stated for the unspecified ground-return input in the total-loop formula. | Stated — p. 3 ``Z_g`` paragraph and p. 4 mode-cancellation discussion. |
| Conductor and insulation geometry | Circular pipe with inner/outer radii ``r_{p1},r_{p2}``, outer insulation radius ``r_{p3}``; solid core radius ``r_{1k}``, insulation radius ``r_{2k}``, pipe-axis offset ``d_k``, pair-axis separation ``d_{km}``, angle ``\theta_{km}``. Symmetrical mode formulas use equal cores. | Stated — Fig. 1/Table 1 p. 1; (5)–(7),(10),(18), pp. 2–3; mode sections p. 4. |
| Constitutive and field assumptions | Scalar uniform conductor/pipe conductivity and relative permeability; conductive skin argument has no displacement-current term. General pipe formulas retain ``\mu_{rp}``; the specified nonmagnetic total/mode reductions set relative permeabilities to unity. Core-to-core proximity is not fully represented by these pipe terms; its separate correction is not silently included. | Equation-implied — (7),(10),(12),(13),(24); stated — sections IV.C–E and VI. |
| Conventions | Positive ``j\omega`` factors, natural logarithms, per-length impedances. Source prints ``\lim`` with frequency-dependent right-hand sides. ``Z_{pm}`` enters connection impedance as ``-2Z_{pm}``. Common mode has equal core currents with pipe return; differential mode cancels pipe-connection and ground terms. | Stated/equation-implied — (17),(22), sections IV.D–E, pp. 3–4. |

**Expression.** Finite-pipe low-frequency expressions, manuscript (12),(13),(22). The source's ``\lim`` notation is preserved, including the surviving ``\omega`` on the right.

```math
\lim_{\omega\to0} Z_{pi}=R_{p,dc}
+\frac{j\omega\cdot\mu_0\mu_{rp}}{2\pi}\cdot
\left(
\frac{r_{p2}^4}{(r_{p2}^2-r_{p1}^2)^2}\cdot
\ln\left(\frac{r_{p2}}{r_{p1}}\right)
-\frac{3r_{p2}^2-r_{p1}^2}{4(r_{p2}^2-r_{p1}^2)}
\right).
\qquad\text{(12)}
```

```math
\lim_{\omega\to0} Z_{p\Sigma,km}
=\frac{j\omega\cdot\mu_0}{2\pi}\cdot\frac{2\mu_{rp}}{1+\mu_{rp}}\cdot
\ln\left(
\frac{r_{p1}^2}
{\sqrt{r_{p1}^4+(d_k\cdot d_m)^2-2r_{p1}^2\cdot d_k\cdot d_m\cdot\cos\theta_{km}}}
\right).
\qquad\text{(13)}
```

```math
\lim_{\omega\to0} Z_p
=\lim_{\omega\to0}(Z_{pi}+Z_{po}-2Z_{pm})
=\frac{j\omega\mu_0\cdot\mu_{rp}}{2\pi}\cdot\ln(r_{p2}/r_{p1}).
\qquad\text{(22)}
```

``R_{p,dc}`` is the pipe DC resistance per length, explicitly named after (12); the source does not supply an additional area formula for it here. ``Z_{pi},Z_{po},Z_{pm}`` are inner, outer, and mutual pipe-surface impedances; ``Z_{p\Sigma,km}`` is pipe-mediated core coupling, not the separate core-to-core proximity correction. ``d_k,d_m`` are radial offsets from the pipe axis, not burial depths, and ``\theta_{km}`` is the angle between those offsets. All radii/distances are lengths, ``\mu_{rp}`` is relative pipe permeability, and ``\mu_0`` has permeability units. No self singularity is inferred from setting equal cable indices in the complete (13).

The source's solid-round-core low-frequency dependency is (19):

```math
\lim_{\omega\to0} Z_{co}
=R_{c,dc}+\frac{j\omega\mu_0}{2\pi}\cdot\frac{\mu_r}{4}.
\qquad\text{(19)}
```

``R_{c,dc}`` is the core DC resistance per length and ``\mu_r`` its relative permeability in this dependency. This is a reproduced core limit, not a new priority claim for Høidalen.

For the source's unity-permeability total cable-loop reductions:

```math
\lim_{\omega\to0,\,\mu_r=1} Z_{kk}
=R_{c,dc}+\frac{j\omega\mu_0}{2\pi}\cdot
\left[\frac14+\ln(r_{p3}/r_{1k})\right]+Z_g(0),
```

```math
\lim_{\omega\to0,\,\mu_r=1} Z_{km}
=\frac{j\omega\mu_0}{2\pi}\cdot[\ln(r_{p3}/d_{km})]+Z_g(0).
\qquad\text{(24)}
```

``Z_g(0)`` is the source's unspecified low-frequency ground-return term, not evaluated or replaced here. These two formulas are total-loop expressions; they are not relabelled pure pipe-surface impedance.

For the symmetrical three-core arrangement, the source additionally prints:

```math
\lim_{\omega\to0,\,\mu_r=1} Z_1
=R_{c,dc}+\frac{j\omega\mu_0}{2\pi}\cdot
\ln(d_{km}/0.78\cdot r_1),
\qquad\text{(31)}
```

```math
\lim_{\omega\to0,\,\mu_r=1} Z_0
=R_{c,dc}+3R_{p,dc}+\frac{j\omega\mu_0}{2\pi}\cdot
\ln\left(\frac{r_{p1}^3}{0.78\cdot r_1\cdot d_{km}^2}\right)
+j\Delta X_0,
\qquad\text{(34)}
```

```math
\lim_{\omega\to0,\,\mu_r=1,\,r_{p2}\to r_{p1}}\Delta X_0
=\frac{\omega\mu_0}{2\pi}\cdot\frac{r_{p2}-r_{p1}}{r_{p1}}.
\qquad\text{(35)}
```

``r_1`` is the common solid-core radius in these mode formulas. The author states that ``\Delta X_0`` comes from (12); (35) is its further thin-wall expression. Equation (31) prints an inline slash followed by ``\cdot r_1``: no parentheses are inserted to turn it into a different radius ratio. The suspected grouping problem is documented below. No explicit conversion of the printed number ``0.78`` to an exponential is made.

**Approximation.** The author takes low-frequency limits of the finite-pipe Bessel model (8)–(11), with pipe-surface connection terms (20)–(21), and combines them through the printed cable/mode relations. The natural skin arguments ``x_1,x_2`` tend to zero with frequency at fixed radii/materials. The source does not print the full Bessel expansion, discarded order, or error remainder. Equation (35) additionally takes a thin-wall limit. The infinite-wall low-frequency approximation is a different record and is not substituted into this finite-wall model.

**Numerical interpretation.** The finite-wall evaluator uses (12) and the low-frequency expansions of (20)–(21). Their combination recovers (22). At unity relative pipe permeability, (13) gives the same coupling as the finite-wall parent (9)–(11). For magnetic pipes, that parent's limit retains the outer radius through a factor ``[1-r(r_{p1}/r_{p2})^{2n}]/[1-r^2(r_{p1}/r_{p2})^{2n}]``, with ``r=(\mu_{rp}-1)/(\mu_{rp}+1)``. The evaluator retains this factor; the printed (13) omits it and agrees only in the nonmagnetic or infinite-outer-radius limits. The author-printed equations above remain unchanged.

**Limitations.** Circular finite pipe and enclosed cores; mode formulas assume the symmetrical three-core arrangement. The evaluator supplies pipe surface and cavity terms, not a separately fitted modal formula or the unspecified ``Z_g(0)``. General unequal-core proximity is not covered by these pipe terms.

**Reference.** [Hoidalen2013](@cite).  Høidalen (2013), DOI `10.1109/TPWRD.2013.2272343`; actual manuscript pp. 2–4, (8)–(13),(17)–(24),(29),(31)–(35). Earlier finite-pipe expressions are explicitly attributed there to Kane et al. [3] and da Silva et al. [5]; their originality is not reassigned to the later witness.

**Transcription source.** Original-author low-frequency manuscript. Page images 2–4 were inspected, with a high-resolution crop for (12)–(13). Earlier Bessel parent formulas are transcribed below as secondary witnesses in Høidalen's notation, not retroactively verified against their originals. The source's geometry was image-checked on p. 1.

## Source transcription

The low-frequency expressions in the formula section keep their original notation. The following parent set, Høidalen (8)–(11),(20)–(21), supplies the complete finite-pipe surface and coupling model from which he takes those limits. The source identifies the finite-wall model with Kane et al. [3], further elaborated in da Silva et al. [5].

```math
Z_{pi}=\frac{j\omega\mu_0}{2\pi}\cdot\frac{\mu_{rp}}{x_1}
\frac{I_0(x_1)\cdot K_1(x_2)+I_1(x_2)\cdot K_0(x_1)}
{I_1(x_2)\cdot K_1(x_1)-I_1(x_1)\cdot K_1(x_2)}.
\qquad\text{(8)}
```

```math
Z_{p\Sigma,km}=\frac{j\omega\mu_0}{2\pi}\cdot
\sum_{n=1}^{\infty}\frac{2\mu_{rp}\cdot C_n}{x_1}\cdot
\frac{A_{nK2}\cdot I_n(x_1)-A_{nI2}\cdot K_n(x_1)}
{A_{nI1}\cdot A_{nK2}-A_{nI2}\cdot A_{nK1}}.
\qquad\text{(9)}
```

The arguments are (7),(10); the coefficient geometry is (5):

```math
\begin{aligned}
x_1&=r_{p1}\cdot\sqrt{j\omega\cdot\mu_0\cdot\mu_{rp}\cdot\sigma_p} \\
x_2&=r_{p2}\cdot\sqrt{j\omega\cdot\mu_0\cdot\mu_{rp}\cdot\sigma_p},
\end{aligned}
```

```math
C_n=(d_k\cdot d_m/r_{p1}^2)^n\cdot\cos(n\theta_{km}).
\qquad\text{(5)}
```

```math
A_{nI1}=n(\mu_{rp}+1)/x_1\cdot I_n(x_1)-I_{n-1}(x_1),
```

```math
A_{nI2}=n(\mu_{rp}-1)/x_2\cdot I_n(x_2)+I_{n-1}(x_2),
```

```math
A_{nK1}=n(\mu_{rp}+1)/x_1\cdot K_n(x_1)+K_{n-1}(x_1),
```

```math
A_{nK2}=n(\mu_{rp}-1)/x_2\cdot K_n(x_2)-K_{n-1}(x_2).
\qquad\text{(11)}
```

```math
Z_{po}=\frac{j\omega\mu_0}{2\pi}\cdot\frac{\mu_{rp}}{x_2}
\frac{I_0(x_2)\cdot K_1(x_1)+I_1(x_1)\cdot K_0(x_2)}
{I_1(x_2)\cdot K_1(x_1)-I_1(x_1)\cdot K_1(x_2)},
\qquad\text{(20)}
```

```math
Z_{pm}=\frac{j\omega\mu_0}{2\pi}\cdot\frac{\mu_{rp}}{x_1\cdot x_2}
\frac{1}{I_1(x_2)\cdot K_1(x_1)-I_1(x_1)\cdot K_1(x_2)}.
\qquad\text{(21)}
```

Høidalen does not define the precise Bessel family in this manuscript. His cited original Kane et al., printed p. 1647 (HAL PDF page 3), explicitly defines ``I_n,K_n`` as modified Bessel functions of first and second kind, order ``n``; that prose was image-checked. It corroborates the parent notation, not a proof that the later witness has no algebraic transcription differences. No branch is stated in Høidalen's square-root definitions. ``\sigma_p`` is pipe conductivity.

For the source-provided cable assembly, (17) is:

```math
Z_{kk}=Z_{co}+Z_{ins1,k}+Z_{ins2,kk}+Z_{pi}+Z_{p\Sigma,kk}
+Z_{po}-2Z_{pm}+Z_{ins3}+Z_g,
```

```math
Z_{km}=Z_{ins2,km}+Z_{pi}+Z_{p\Sigma,km}
+Z_{po}-2Z_{pm}+Z_{ins3}+Z_g.
\qquad\text{(17)}
```

Its required insulation terms, (4),(6),(18), are kept as assembly dependencies rather than attributed as new Høidalen insulation formulas:

```math
Z_{ins2,km}=\frac{j\omega\mu_0}{2\pi}\cdot Q_{km},\qquad\text{(4)}
```

```math
Q_{kk}=\ln\left(\frac{r_{p1}}{r_{2k}}\cdot[1-(d_k/r_{p1})^2]\right),
```

```math
Q_{km}=\ln\sqrt{
\frac{r_{p1}^4+(d_k\cdot d_m)^2-2\cdot r_{p1}^2\cdot d_k\cdot d_m\cdot\cos\theta_{km}}
{r_{p1}^2\cdot(d_k^2+d_m^2-2\cdot d_k\cdot d_m\cdot\cos\theta_{km})}}.
\qquad\text{(6)}
```

```math
\begin{aligned}
Z_{ins1,k}&=\frac{j\omega\mu_0}{2\pi}\cdot\ln(r_{2k}/r_{1k}) \\
Z_{ins3}&=\frac{j\omega\mu_0}{2\pi}\cdot\mu_{ri3}\cdot\ln(r_{p3}/r_{p2}).
\end{aligned}\qquad\text{(18)}
```

``\mu_{ri3}`` is relative permeability of the outer insulation in (18). It is not restored to the first insulation term, which prints only ``\mu_0``. Full-frequency ``Z_{co}`` and ``Z_g`` are named external dependencies; only the core's low-frequency (19) is supplied here. Thus (17) does not constitute a new completely specified full-frequency cable formula in this record.

The source's symmetrical-mode assembly, p. 4, is:

```math
Z_1=Z_{co}+Z_{ins1}+(Z_{p\Sigma,kk}-Z_{p\Sigma,km})
+\frac{j\omega\mu_0}{2\pi}\cdot(Q_{kk}-Q_{km}),
\qquad\text{(29)}
```

```math
Z_0=Z_{co}+Z_{ins1}+3Z_{pi}+Z_{p\Sigma,kk}+2Z_{p\Sigma,km}
+\frac{j\omega\mu_0}{2\pi}\cdot(Q_{kk}+2Q_{km}).
\qquad\text{(32)}
```

The source states that connection and ground-return terms cancel in these modes; it does not set the physical earth impedance identically zero in the whole cable model. Section V.A, across the bottom of p. 4 and top of p. 5, gives a geometry-specific series observation: in the symmetrical arrangement, all terms except ``n=3,6,9,\ldots`` contribute to ``Z_{p\Sigma,kk}-Z_{p\Sigma,km}``, whereas only those multiples of three contribute to ``Z_{p\Sigma,kk}+2Z_{p\Sigma,km}``. The complete sentence spans the page break; it does not say that only multiples of three contribute to the differential term. The full parent summations are retained above. The separate core-to-core correction and its mode factors are in the [proximity record](../solid-core-proximity-low-frequency-subtracted-series/Hoidalen2013.md).

## Notation map

No notation renaming is used. The plain source ``\lim`` and the ambiguous inline grouping in (31) are retained.

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{pi},Z_{po},Z_{pm}`` | unchanged | Pipe inner, outer, mutual surface impedances | Per length; mutual term enters as ``-2Z_{pm}`` |
| ``Z_p`` | unchanged | ``Z_{pi}+Z_{po}-2Z_{pm}`` | Pipe surface-connection combination, not total cable impedance |
| ``Z_{p\Sigma,km}`` | unchanged | Pipe-mediated core self/mutual coupling | Per length; separate from core-to-core proximity |
| ``Z_{kk},Z_{km},Z_1,Z_0`` | unchanged | Total loop self/mutual and differential/common-mode impedance | Per length; formula-specific restrictions retained |
| ``Z_{co},R_{c,dc},R_{p,dc}`` | unchanged | Core internal impedance, core DC resistance, pipe DC resistance | Per length; DC resistances are named inputs |
| ``Z_{ins1,k},Z_{ins2,km},Z_{ins3}`` | unchanged | Core insulation, internal pipe region, outer pipe insulation terms | Source-provided assembly dependencies |
| ``Z_g,Z_g(0)`` | unchanged | Ground-return impedance and printed zero-frequency value | Unspecified by this source; not inferred as zero |
| ``r_{p1},r_{p2},r_{p3}`` | unchanged | Pipe inner, pipe outer, outer-insulation radii | Length |
| ``r_{1k},r_{2k},r_1`` | unchanged | Core/insulation radii; common core radius in mode limits | Length |
| ``d_k,d_m,d_{km},\theta_{km}`` | unchanged | Pipe-axis offsets, pair-axis separation, angular separation | Lengths and angle; not earth coordinates |
| ``\mu_0,\mu_{rp},\mu_r,\mu_{ri3}`` | unchanged | Reference permeability and relative pipe/core/outer-insulation values | ``\mathrm{H/m}`` and dimensionless ratios |
| ``\sigma_p,x_1,x_2`` | unchanged | Pipe conductivity and skin arguments | ``\mathrm{S/m}``; dimensionless arguments |
| ``C_n,A_{nI1},A_{nI2},A_{nK1},A_{nK2}`` | unchanged | Geometry and boundary coefficients | Dimensionless; (5),(11) |
| ``I_n,K_n`` | unchanged | Modified Bessel functions, first/second kind | Definition corroborated in the cited Kane original, not stated in Høidalen |
| ``Q_{kk},Q_{km}`` | unchanged | Logarithmic insulation-region geometry factors | Dimensionless |
| ``\Delta X_0`` | unchanged | Residual common-mode pipe reactance term | Per-length reactance; (35) is further thin-wall formula |
| ``k,m,n,\omega,j`` | unchanged | Core labels, positive-integer series order, angular frequency, imaginary unit | ``n=1`` to infinity in parent (9); ``j^2=-1`` |

## Evidence and approximation sources

- Finite-wall parent: (8)–(11), pp. 2–3; Høidalen expressly attributes it to Kane [3] and da Silva [5]. The new recorded object here is his low-frequency analysis, not a priority claim for those Bessel formulas.
- Primary low-frequency pipe formulas: (12),(13),(22); (19) supplies the separately identified core dependency. Parent finite-wall geometry and conductive skin arguments are explicit; the source does not give a full asymptotic remainder calculation.
- Total and modal assembly: (17),(24),(29),(31),(32),(34),(35). Definitions for all displayed insulation/geometry terms are preserved above. An external ``Z_g`` is intentionally not filled using a remembered Pollaczek formula.
- Tests: sections V–VI compare finite/infinite pipe models and FEM, and distinguish common/differential mode behavior. Section V.A, p. 4, uses an upper summation limit of 20 for the subsequent pipe-series calculations after discussing convergence of (3) and (9). This numerical test setting does not replace the printed infinite series or define a low-frequency expansion order. Agreement in the selected tests does not define a universal error bound.
- The source warns that combining an infinite-wall inner term with finite-wall connection terms can produce negative low-frequency self resistance. That different model is retained in the [separate infinite-wall approximation record](../infinite-pipe-low-frequency-logarithmic-terms/Hoidalen2013.md), not mixed into the finite-wall result.

## Limitations and discrepancies

- **Published limit notation:** several equations print ``\lim_{\omega\to0}`` but keep an order-``\omega`` term on the right. They are retained as the source's low-frequency expressions, not replaced by equalities with ``\omega=0`` or reviewer-generated asymptotic notation.
- **Suspected published grouping defect:** (31) prints ``\ln(d_{km}/0.78\cdot r_1)``; the source calls it the GMD solution, while that literal slash/product grouping does not state a dimensionless radius ratio. The original grouping is preserved, not rewritten as ``d_{km}/(0.78r_1)``.
- **Output distinctions:** total-loop (24) has unspecified ``Z_g(0)``; it is not the same quantity as (12),(13),(22). Source mode cancellation does not establish an earth-property assumption.
- **No isolated-term singularity inference:** (12)'s factors with ``r_{p2}^2-r_{p1}^2`` are part of a complete expression at finite thickness. No physical singularity is diagnosed from a factor alone, and no stabilized or reconstructed thin-wall limit is substituted.
- **Conversion defects:** Markdown damages many equation rows and mixes columns. Source page images settle the square root in (13), powers/radius ordering in (12), all parent coefficients, and the inline grouping of (31); the conversion is not used to repair the source.
