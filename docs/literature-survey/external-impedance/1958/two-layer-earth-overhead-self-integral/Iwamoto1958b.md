# Two-layer overhead self earth-return integral — Iwamoto witness

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Overhead line-current self geometry specified by height; conductor skin, insulation, and free-space geometric terms are excluded. |
| Calculated quantities | Two-layer self earth-return impedance, nonconducting-lower-region specialization, normalization, and graphical evaluator. |
| Earth structure | An upper layer of resistivity ``\rho_1`` and thickness ``D`` above a lower half-space of resistivity ``\rho_2``; ``\rho_2=\infty`` gives a nonconducting lower region. |
| Model and approximation | Conduction-only, fixed-permeability integral attributed to Sunde. The chart uses the approximate ``r_e\sim0.1\sqrt{f/R_e}``; it does not replace printed equation (21). |
| Main source | K. Iwamoto (August 1958), citing Sunde (1949). |
| Citation key(s) | Primary: `:Iwamoto1958b`; auxiliary evaluator: `:Iwamoto1958a`; attributed sources: `:Sunde1949`, `:Carson1926` |
| Evidence status | All 12 August pages checked. The root branch, unit scaling, and the January operator's step argument remain unresolved. |

**Description.** Self earth-return impedance of an overhead wire above a finite conducting earth layer and a lower earth half-space, including a nonconducting lower-region specialization and normalized resistance/reactance representation.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not stated as an independent imposed parameter in the appendix. There is no independent longitudinal variable in the final kernel. The main paper's solved line-propagation quantity is not an impressed-current prescription; it does not establish ``Γ=0`` here. | Not stated — appendix, pp. 1048–1049. Equation-implied — final (付1) contains only the displayed transverse spectral variable. |
| Air propagation constant ``γ_air`` | No air bulk propagation constant or air wavenumber is defined in the appendix kernel. Do not restore one or assert a separately printed zero value. Air permeability is ``\mu_0``. | Not stated — appendix. Stated — Fig. 付第1, p. 1048. |
| Earth propagation constant ``γ_earth`` | No named bulk gamma is supplied. The source prints normalized radicals ``\sqrt{\lambda^2+j}`` and ``\sqrt{\lambda^2+j\rho_1/\rho_2}`` with scales ``D',h'`` below. These are retained without manufacturing another gamma definition. The square-root branch is not specified. | Stated — definitions after (付1), p. 1048. Unresolved — explicit root branch. |
| Earth permittivity and displacement current | Earth permittivity is neglected in the Carson-based appendix construction; the radicals and scales retain resistivity and permeability only. No dielectric-loss or conductivity/permittivity double-counting rule is relevant to this conduction-only formula. | Stated — opening appendix paragraph, p. 1048, “誘電率を無視した”. Equation-implied — (付1) definitions omit a permittivity term. |
| Range of validity | Horizontal two-layer model with the stated omitted permittivity and fixed permeability. No universal frequency cutoff or analytical error bound is printed. The finite-depth chart calculations enumerate ``n'`` values, not a proof of uniform validity for all parameters. | Stated — appendix geometry/model and p. 1049 numerical procedure. Not stated — universal bound or truncation error. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0=4\pi\times10^{-7}`` in air and both earth regions; ``\mu_{12},\mu_{10}`` are dimensionless reflection coefficients, **not** layer permeabilities. | Stated — Fig. 付第1 and definitions after (付1), p. 1048. |
| Arrangement | Single overhead wire at height ``h`` above the top interface; self earth-return contribution. No source/target mutual separation or buried/mixed formula is printed in this appendix. | Stated — Fig. 付第1 and appendix heading, p. 1048. Equation-implied — ``e^{-2h'\lambda}`` and absence of a mutual-placement variable. |
| Earth structure | Upper earth layer of resistivity ``\rho_1``, thickness ``D``; lower half-space of resistivity ``\rho_2``. Finite conducting depth is represented by ``\rho_2=\infty``, hence ``\mu_{12}=\mu_{10}``. This is a nonconducting lower region, not a perfectly conducting one. | Stated — Fig. 付第1 and prose after (付1), p. 1048. |
| Conductor and insulation geometry | Overhead line-current geometry specified by height. No wire radius or insulation thickness enters the appendix correction, and no conductor-internal or free-space geometric term is added. No inferred self-radius regularization is used. | Equation-implied — Fig. 付第1 and (付1). Not stated — finite-radius correction for this formula. |
| Constitutive and field assumptions | Each layer uses one scalar resistivity; horizontal uniformity and a linear current-to-impedance relation are equation-implied. The author invokes extension of Carson's calculation, not a newly specified full-wave boundary problem. Nonlinear, anisotropic and dispersive constitutive models are not supplied. | Stated — appendix Carson attribution. Equation-implied — uniform scalar layer labels and (付1). |
| Conventions | Vertical coordinate is ``y``, surface ``y=0``, lower interface ``y=-D``. Preserve positive ``j`` in the radicals; a full time exponential is not stated. The dot on ``\dot Z_e`` denotes the source's complex impedance, not a time derivative. MKS dimensions give Ω/m; ``R_e`` is Ω/m, while chart quantity ``R'_e=1000R_e`` is Ω/km. | Stated — Fig. 付第1, (付1) MKS label, §7 (21), Fig. 付第3 inset. Equation-implied — complex chart output and per-length normalization. |

**Expression.** The self earth-return impedance, appendix (付1), p. 1048. In particular, the prefactor has **no additional ``j``**:

```math
\dot Z_e=\frac{\omega\mu_0}{\pi}
\int_0^\infty \{\sqrt{\lambda^2+j}-\lambda\}
\frac{1+\mu_{12}e^{-2D'\sqrt{\lambda^2+j}}}
     {1-\mu_{12}\mu_{10}e^{-2D'\sqrt{\lambda^2+j}}}
e^{-2h'\lambda}\,d\lambda .
\tag{付1}
```

```math
\mu_{12}=
\frac{\sqrt{\lambda^2+j}-\sqrt{\lambda^2+j\rho_1/\rho_2}}
     {\sqrt{\lambda^2+j}+\sqrt{\lambda^2+j\rho_1/\rho_2}},
\qquad
\mu_{10}=\frac{\sqrt{\lambda^2+j}-\lambda}
              {\sqrt{\lambda^2+j}+\lambda},
```

```math
D'=\sqrt{\omega\mu_0/\rho_1}\,D,\qquad
h'=\sqrt{\omega\mu_0/\rho_1}\,h,\qquad
\mu_0=4\pi\times10^{-7}.
```

The source-prescribed finite-depth specialization is ``\rho_2=\infty`` and ``\mu_{12}=\mu_{10}``. The normalized finite-depth expression and all its helpers, §7 and appendix:

```math
n=D/h,\tag{20}
```

```math
r_e=4\sqrt{\omega\times10^{-7}/R_e}
    =0.10026\sqrt{f/R_e},\tag{21}
```

```math
R_e=\frac{\rho_e}{\pi h^2},\qquad R'_e=1000R_e,\qquad n'=nr_e,
```

```math
F(\lambda)=\{\sqrt{\lambda^2+j}-\lambda\}
\frac{1+\mu_{10}e^{-n'\sqrt{\lambda^2+j}}}
     {1-\mu_{10}^{2}e^{-n'\sqrt{\lambda^2+j}}},
```

```math
\frac{\dot Z_e(r_e)}{R_e}
=\frac{r_e^2}{4}\int_0^\infty F(\lambda)e^{-r_e\lambda}\,d\lambda .
\tag{付2}
```

Here ``\rho_e`` is the conducting earth's resistivity in the finite-depth chart model (the upper-layer ``\rho_1`` by geometry/context), not an independently fitted second-layer quantity. ``f`` is frequency and ``\omega`` angular frequency. The printed decimal and equality in (21) are retained. A 450-dpi crop confirms that both denominators are unprimed. **Suspected published unit-switch:** under the conventional angular/cyclic-frequency relation and the same numerical resistance scale, the two printed prefactors do not agree (their ratio is approximately the square root of 1000). The separately printed per-kilometre scale may explain the mismatch, but this is an inference, not a source-authorized correction; the exact normalization convention remains unresolved. ``n`` here is a depth/height ratio, not the conductor count denoted by the same letter in §6.

The separately printed evaluation relation is retained, with only the crossed 相乗 operator glyph displayed as ``\star`` (notation-only alias; definition and its printed discrepancy are transcribed below):

```math
\frac{\dot Z_{el}(u)}{R_e}
=\frac{r_e}{4}F_l(u)\star E_l(u),\tag{付3}
```

```math
\dot Z_{el}(u)=\dot Z_e(e^u),\qquad
F_l(u)=F(e^{-u}),\qquad
E_l(u)=1-e^{-e^u}.
```

The lower-case subscript ``l`` is preserved. The operator is not ordinary multiplication. Its separately inspected January definition (18), §4, p. 25, supplies the full domain and endpoint term:

```math
\varphi(u)=\varphi_0(u)H(u-\infty),\qquad
\psi(u)=\psi_0(u)H(u-\infty),
```

```math
\varphi(u)\star\psi(u)=\varphi_0(-\infty)\psi_0(\infty)
+\int_{-\infty}^{\infty}\varphi'_0(u-\zeta)\psi_0(\zeta)\,d\zeta.
\tag{18}
```

Here ``\varphi,\psi`` are generic operands, subscript 0 denotes the printed helper functions, ``H`` is the unit step, and the prime is differentiation. The accompanying prose specifies a discontinuity at **negative** infinity despite the literal ``H(u-\infty)`` helpers. This inconsistency remains unresolved; the endpoint term is not dropped or altered by inference. The fuller source explanation and separate witness locator are retained below.

**Approximation.** (付1) is an integral representation within the author's conduction-only, fixed-permeability layer model, not a newly introduced analytical approximation. The finite-depth choice is an explicitly stated material specialization. (付2)–(付3) provide normalization and a graphical/numerical evaluation route, not a separate physical contribution merely because notation changes. The chart inset uses the expressly approximate ``r_e\sim0.1\sqrt{f/R_e}``; it must not replace the printed (21) in a source transcription. Section 7's effective penetration-depth fits are measurement outputs, not new impedance kernels.

**Limitations.** Secondary verification of the credited Sunde kernel is limited to Iwamoto's printed witness. The inspected Sunde source is a corrected later edition. No explicit root branch, imposed longitudinal prescription, general mutual geometry or full external self assembly is supplied here. Equation (21)'s unprimed frequency forms have an unresolved unit-scale mismatch. The January operator's prose places a discontinuity at negative infinity while its printed step arguments use a minus-infinity subtraction; retain that inconsistency rather than changing a sign.

**Reference.** [Iwamoto1958b](@cite); [Iwamoto1958a](@cite). Iwamoto, August 1958, appendix (付1)–(付3), pp. 1048–1049, with §7 (20)–(21), p. 1045. Iwamoto, January 1958, §4 (18), p. 25, for the separately inspected operator definition. [Sunde1949](@cite), p. 119 as credited by Iwamoto, **1949 original not inspected**; the inspected witness is the 1968 corrected edition. [Carson1926](@cite) is the expressly credited parent, not the source of Iwamoto's layer coefficients.

**Transcription source.** All displayed August equations were checked against the original August publication page images, including magnified coefficient and exponent crops. They remain a **secondary witness for the earlier-attributed kernel**; original-page inspection of Iwamoto does not establish priority over Sunde. January (18) is separately checked against that earlier Iwamoto publication, not reconstructed from the August transform.

## Source transcription

Notation is unchanged except for the explicitly mapped 相乗 glyph. Source labels retain the original Japanese appendix prefix. The complete source equation set is repeated here in source order.

### August §7, p. 1045/PDF page 8

```math
n=D/h.\tag{20}
```

```math
r_e=4\sqrt{\omega\times10^{-7}/R_e}
    =0.10026\sqrt{f/R_e}.\tag{21}
```

The definition ``R_e=\rho_e/(\pi h^2)`` is printed in the main paper's conductor-resistance definitions and repeated in Fig. 付第3. No reinterpretation as the real part of ``\dot Z_e`` is intended.

### August appendix, p. 1048/PDF page 11

```math
\dot Z_e=\frac{\omega\mu_0}{\pi}
\int_0^\infty\{\sqrt{\lambda^2+j}-\lambda\}
\frac{1+\mu_{12}e^{-2D'\sqrt{\lambda^2+j}}}
     {1-\mu_{12}\mu_{10}e^{-2D'\sqrt{\lambda^2+j}}}
e^{-2h'\lambda}\,d\lambda
\quad\text{(MKS)}.\tag{付1}
```

```math
\mu_{12}=
\frac{\sqrt{\lambda^2+j}-\sqrt{\lambda^2+j\rho_1/\rho_2}}
     {\sqrt{\lambda^2+j}+\sqrt{\lambda^2+j\rho_1/\rho_2}},
\quad
\mu_{10}=\frac{\sqrt{\lambda^2+j}-\lambda}
              {\sqrt{\lambda^2+j}+\lambda},
```

```math
D'=\sqrt{\omega\mu_0/\rho_1}\,D,\quad
h'=\sqrt{\omega\mu_0/\rho_1}\,h,\quad
\mu_0=4\pi\times10^{-7}.
```

The intervening prose credits the same form to Sunde [7] and explicitly gives ``\rho_2=\infty``, ``\mu_{12}=\mu_{10}`` for finite conducting depth.

```math
\frac{\dot Z_e(r_e)}{R_e}
=\frac{r_e^2}{4}\int_0^\infty F(\lambda)e^{-r_e\lambda}\,d\lambda.
\tag{付2}
```

```math
\frac{\dot Z_{el}(u)}{R_e}
=\frac{r_e}{4}F_l(u)\star E_l(u).\tag{付3}
```

```math
\dot Z_{el}(u)=\dot Z_e(e^u),\qquad F_l(u)=F(e^{-u}).
```

### August appendix continuation, p. 1049/PDF page 12

```math
E_l(u)=1-e^{-e^u},
```

```math
F(\lambda)=\{\sqrt{\lambda^2+j}-\lambda\}
\frac{1+\mu_{10}e^{-n'\sqrt{\lambda^2+j}}}
     {1-\mu_{10}^{2}e^{-n'\sqrt{\lambda^2+j}}},
\qquad n'=nr_e.
```

The numerical procedure expressly enumerates ``n'=0.1,0.15,0.2,0.3,0.5,0.8,1,2,3,5,7``, then uses ``n=n'/r_e`` to construct the charts. This is a tested/chart-construction list, not a range-of-validity theorem.

Fig. 付第2's inset repeats the finite-depth normalized impedance with **unsubscripted** reflection coefficient ``\mu``:

```math
P_e+jQ_e=\frac{r_e^2}{4}
\int_0^\infty\{\sqrt{\lambda^2+j}-\lambda\}
\frac{1+\mu e^{-nr_e\sqrt{\lambda^2+j}}}
     {1-\mu^2 e^{-nr_e\sqrt{\lambda^2+j}}}
e^{-r_e\lambda}\,d\lambda,
```

```math
\mu=\frac{\sqrt{\lambda^2+j}-\lambda}{\sqrt{\lambda^2+j}+\lambda},
\qquad r_e\sim0.1\sqrt{f/R_e}.
```

Fig. 付第3's unit box prints:

```math
R'_e=1000R_e,\qquad R_e=\rho_e/(\pi h^2).
```

The inset labels ``\rho_e`` in Ωm and ``D`` in m. The chart paragraph uses a per-kilometre conductor-resistance scale; the primed scale must not silently replace the per-metre quantity in (21). ``P_e,Q_e`` are the normalized impedance's resistance/reactance parts, not potential coefficients.

### Separately attributed January operator dependency

January §4, p. 25/PDF page 4/index 3, defines the same crossed 相乗 symbol as a Duhamel operation, extended to functions with a discontinuity at ``u=-\infty``. The **literal printed** helper lines are:

```math
\varphi(u)=\varphi_0(u)H(u-\infty),\qquad
\psi(u)=\psi_0(u)H(u-\infty).
```

```math
\varphi(u)\star\psi(u)
=\varphi_0(-\infty)\psi_0(\infty)
+\int_{-\infty}^{\infty}
\varphi'_0(u-\zeta)\psi_0(\zeta)\,d\zeta .
\tag{18}
```

The equation number in the publication is **(18)**; this subsection's January locator distinguishes it from August's unrelated (18). ``H`` is the unit step, as described after January (12), p. 24/PDF page 3/index 2; ``\varphi,\psi`` are generic functions and the prime is the source's differentiation notation. The following prose says that the endpoint term vanishes when either or both functions have no discontinuity at negative infinity, and the subscript 0 can then be removed.

**Printed discrepancy:** the prose says **negative** infinity, but the helpers visibly print ``H(u-\infty)``. The endpoint factor visibly uses ``\psi_0(\infty)``, not ``\psi_0(-\infty)``. Neither is silently repaired. This record does not select altered endpoints or assert a resolved ordinary-convolution implementation of (付3).

## Notation map

Only the crossed 相乗 glyph is renamed, one-to-one to ``\star``. All other source notation is retained. The repeated use of ``\mu`` is context-specific, not consolidated into one permeability symbol.

| Source → display | Meaning and units | Convention / locator |
| --- | --- | --- |
| ``\dot Z_e,\dot Z_e(r_e)`` → same | Complex self earth-return impedance, Ω/m | Dot is complex-quantity notation; (付1)–(付2) |
| ``\omega,f,j`` → same | Angular frequency, frequency (s⁻¹/Hz), imaginary unit | Positive ``j`` retained in radicals; explicit time exponential not stated |
| ``\mu_0`` → same | Fixed permeability, H/m in MKS | Air and both earth layers, Fig. 付第1 |
| ``\rho_1,\rho_2`` → same | Upper/lower earth resistivity, Ωm | Layer 1 finite, layer 2 half-space |
| ``\rho_e`` → same | Conducting-layer resistivity in finite-depth charts, Ωm | Upper-layer identification is geometry/context-implied |
| ``\mu_{12},\mu_{10}`` → same | Dimensionless reflection coefficients | Definitions after (付1); not layer permeabilities |
| Chart ``\mu`` → same | Dimensionless coefficient equal by its definition to ``\mu_{10}`` | Fig. 付第2 inset only |
| ``y,h,D`` → same | Vertical coordinate, wire height, upper-layer thickness, m | Surface 0, lower interface −D; do not change vertical coordinate to the later paper's z |
| ``\lambda`` → same | Dimensionless integration variable | Domain 0 to infinity, measure ``d\lambda`` |
| ``h',D'`` → same | Dimensionless scaled height/thickness | Prime denotes scaling, not differentiation |
| ``R_e,R'_e`` → same | “Conductor resistance” scale, Ω/m and Ω/km respectively | Not simply the real part of the earth impedance |
| ``n,n',r_e`` → same | Dimensionless depth ratio and normalized parameters | Here ``n=D/h``; §6's wire count is a different use |
| ``F,F_l,E_l`` → same | Dimensionless spectral and transformed helper functions | Exact arguments ``e^{-u}`` and nested exponential preserved |
| ``u`` → same | Logarithmic transformed coordinate, dimensionless | ``\dot Z_{el}(u)=\dot Z_e(e^u)`` |
| ``\dot Z_{el}`` → same | Transformed-argument complex impedance, Ω/m | Subscript is lower-case l, not numeral 1 |
| ``P_e,Q_e`` → same | Dimensionless normalized resistance/reactance | Fig. 付第2; not potential/admittance |
| Crossed 相乗 glyph → ``\star`` | Source-defined Duhamel operation | January (18), including endpoint term and printed helper inconsistency |
| ``\varphi,\psi,\varphi_0,\psi_0,H,\zeta`` → same | Generic transform functions, their source helper functions, unit step, and integration coordinate | January operator definition; dimensions depend on operands; ``\zeta`` dimensionless here |
| Prime on ``\varphi'_0`` → same | Differentiation in the operator definition | Not the scaling prime on ``D',h',n',R'_e`` |

## Evidence and approximation sources

1. **Credited parent and priority boundary.** August p. 1048 explicitly credits Carson's infinite-depth construction and states agreement in form with Sunde p. 119 (1949). This record preserves a useful secondary kernel witness. It does not count Iwamoto's new notation or charting as a new physical contribution. The [1973 Nakagawa record](../../1973/three-layer-earth-overhead-integral/Nakagawa1973.md) independently documents the later attribution to Iwamoto and its restricted comparison; that citation is not proof that the 1958 paper first introduced the kernel.
2. **Sunde witness inspected.** The book's copyright page identifies an unabridged, corrected 1968 republication of the 1949 work; its preface reports mostly typographical corrections. Page 119 contains a two-layer integral for mutual inductance (4.55), not Iwamoto's normalized self-impedance notation. The record does not infer algebraic conversion or identity with the uninspected 1949 edition.
3. **Finite-depth specialization.** August p. 1048 sets lower resistivity to infinity and identifies the two reflection coefficients. This is a source-stated material formula, not a thin-layer expansion. The exponent in ``F`` is the product ``-n'\sqrt{\lambda^2+j}``; magnified imagery confirms that the root is not a denominator.
4. **Evaluation route and precision.** The source normalizes the integral and uses its earlier 相乗 method to generate finite-depth charts. It supplies the listed ``n'`` sample values. Section 9, p. 1048, distinguishes roughly three-digit tabulated ordinates from about two useful interpolation digits, with poorer dotted-line portions. That statement concerns the paper's graphical calculations, not a demonstrated global relative-error bound for (付1).
5. **Other August sections.** Sections 1–5 concern resistivity inference and waveforms; §§7–8 use the finite-depth impedance charts for penetration-depth criteria and a measurement case. Their fitted depth outputs and waveform convolutions are not new line-parameter kernels.
6. **Earlier dependency.** The January paper's (18) defines the operation used in August (付3). The [January record](../homogeneous-earth-overhead-logarithmic-integral-evaluation/Iwamoto1958a.md) preserves its distinct logarithmic-integral evaluator, main/appendix/chart disagreements, and normalization. The standard cylindrical-conductor and Carson constituents retain their earlier attributions.

## Limitations and discrepancies

1. **Edition used:** the formula was checked in the corrected 1968 Dover republication of Sunde's 1949 work.
2. **Material and propagation specification:** no named earth bulk gamma, explicit air bulk gamma, imposed longitudinal dependence, or root-branch rule is printed with the appendix. Omission is not recast as a source-stated zero or a newly supplied square-root branch.
3. **Output boundary:** ``\dot Z_e`` is a self earth-return term. No conductor skin impedance, free-space logarithmic term, wire-radius substitution, general mutual cosine factor, or admittance conversion has been appended.
4. **Printed operator inconsistency:** January (18)'s discontinuity prose and ``H(u-\infty)`` helpers disagree in their infinity sign. Their literal tokens, including the positive-infinity endpoint factor, are retained. The main impedance integral should be evaluated directly instead of using this operator.
5. **Suspected published normalization discrepancy:** a 450-dpi image of August (21) confirms unprimed resistance denominators in both equalities. With the conventional angular/cyclic-frequency relation and identical numerical resistance units, their coefficients disagree by approximately the square root of 1000. The separate per-kilometre chart resistance provides a possible explanation, not a verified correction. Neither a prime nor a replacement constant has been inserted.
6. **Notation and units:** reflection coefficients are not permeabilities; ``n`` is not the §6 wire count; ``R'_e`` is not the per-metre normalization; lower-case ``l`` is not numeral 1; the prefactor lacks an additional ``j``. The chart's approximate decimal does not replace (21).
