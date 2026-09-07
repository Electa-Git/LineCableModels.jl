# Papadopoulos–Papagiannis–Labridis two-layer overhead earth impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Two infinite, electrically thin overhead conductors; the self term uses the source-prescribed radius substitution. Conductor internal and insulation terms are excluded. |
| Calculated quantities | Mutual earth-return impedance correction, self correction, homogeneous reduction, and auxiliary 2010 matrix assembly. |
| Earth structure | Air above a finite earth layer of thickness ``d`` and a lower earth half-space. |
| Model and approximation | Integral representation under the thin-conductor quasi-TEM model. The 2010 derivation replaces the unknown longitudinal constant by the air value and applies a Bessel transform; no quadrature is specified. |
| Main source | T. A. Papadopoulos, G. K. Papagiannis, and D. A. Labridis (2009); the 2010 paper supplies the auxiliary derivation. |
| Citation key(s) | Primary: `:Papadopoulos2009`; auxiliary derivation: `:Papadopoulos2010a` |
| Evidence status | Both publications checked against page images. Prime/index conflicts, root branches, and the appendix identity remain unresolved. |

**Description.** Per-unit-length self and mutual earth-return impedance corrections for parallel overhead thin conductors above a finite upper earth layer and a lower earth half-space. A separately identified later witness supplies the perfectly conducting ground contribution and its assembly with the correction.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | 2009: quasi-TEM is stated, but no explicit imposed longitudinal constant or exponential is printed; the definition ``a_m=\sqrt{u^2+\gamma_m^2-\gamma_0^2}`` is not an explicit author statement about ``Γ``. 2010 parent (4a) has unknown ``\gamma_x`` and ``e^{-\gamma_xx}``; final (5) prescribes ``\gamma_x=\gamma_0=jk_0=j\omega\sqrt{\mu_0\varepsilon_0}``. It is neither zero nor the lossy earth bulk constant. | Stated — 2009 I–II.A, p. 1064; unnumbered appendix definition p. 1067; 2010 p. 1162 between (4) and (5). Unresolved — original explicit longitudinal prescription. |
| Air propagation constant ``γ_air`` | 2009 air has ``\varepsilon_0,\mu_0``; the indexed definition is ``\gamma_0^2=j\omega\mu_0(\sigma_0+j\omega\varepsilon_0)`` and no separate value of ``\sigma_0`` is printed there. The ``-\gamma_0^2`` term is retained. 2010 explicitly uses the free-space ``\gamma_0=j\omega\sqrt{\mu_0\varepsilon_0}``, not zero. | Stated — 2009 II.A and appendix; 2010 p. 1162 longitudinal prescription. Equation-implied — nonconducting air in that 2010 prescription. |
| Earth propagation constant ``γ_earth`` | ``\gamma_m^2=j\omega\mu_m(\sigma_m+j\omega\varepsilon_m)``, earth indices 1,2. Retains conductivity and permittivity separately. 2009 ``a_m`` and 2010 ``\alpha_k`` are distinct source notations; their formulas are retained below. Square-root branches are not explicitly given. | Stated — 2009 appendix p. 1067; 2010 appendix p. 1169, below (A.3) and (A.12). |
| Earth permittivity and displacement current | Arbitrary scalar layer permittivities retained in the bulk constants. Axial displacement effects are retained in this impedance model; the radial kernel belongs to the separately recorded potential coefficient. No independent loss tangent or complex-permittivity/conductivity decomposition is specified. | Stated/equation-implied — 2009 II–III, (1)–(2), appendix; 2010 sections 3.1–3.3. |
| Range of validity | Electrically thin perfect parallel conductors and quasi-TEM TL mode; no spectral expansion order or universal numerical error bound. 2010 states the TL criterion ``f\ll c/h`` and discusses conditional use beyond it. Its 50 Hz–10 MHz study is a tested range, not a universal validity interval. | Stated — 2009 p. 1064; 2010 pp. 1161–1163,1168. |
| Earth permeability ``μ_earth`` | Layer-specific scalar ``\mu_1,\mu_2`` retained. Unity relative permeability is a numerical-example choice only. | Stated — 2009 II.B, IV, pp. 1065–1066; 2010 sections 3.3,4.2, pp. 1163–1164. |
| Arrangement | Both conductors in air. Mutual correction for horizontal separation ``y_{ij}``, heights ``h_i,h_j``; self by replacing ``y_{ij}`` with the conductor's outer radius and ``h_j`` with ``h_i``. No buried or mixed-placement formula is provided here. | Stated — 2009 Fig. 1, II.A, pp. 1064–1065; 2010 Fig. 2 and prose after (6), p. 1162. |
| Earth structure | Two horizontal earth layers: upper thickness ``d``, lower infinite extent, air above. Identical-layer prescription gives homogeneous earth. An arbitrary-number-of-layers extension is only proposed in the 2010 prose, not supplied as a recursion. | Stated — 2009 II.A,II.C; 2010 section 3.1, p. 1162. |
| Conductor and insulation geometry | 2009: two uniform, electrically thin, perfect conductors of infinite length. 2010 derivation takes infinite source conductor ``i`` and target ``j`` per unit length. Outer radius regularizes the source-prescribed self geometry; no insulation contribution or internal skin/proximity expression is part of this earth correction. | Stated — 2009 II.A; 2010 section 3.1 and the separate ``Z'_w`` in (2). |
| Constitutive and field assumptions | Uniform scalar electromagnetic properties within each layer; linear Hertz-vector field superposition is equation-implied. Final kernels retain the TL mode under a quasi-TEM restriction, not a solved full-wave dispersion model. Grounding/bonding constraints beyond the reference-conductor TL description are not stated for the kernel. | Stated/equation-implied — 2009 I–II.A; 2010 section 2 and appendix A. |
| Conventions | ``j`` imaginary and ``\omega=2\pi f`` in 2009. No explicit time exponential is printed; positive ``j\omega`` is the equation's harmonic convention. Heights are upward from the upper earth surface. 2010 appendix uses upward ``z`` from the lower interface: earth 1 at ``0\le z<d``, earth 2 at ``z<0``, air at ``z\ge d``. ``Z'_g`` is a per-length correction, not total series impedance. | Stated — 2009 p. 1065/Fig. 1; 2010 (1)–(2), Fig. 18, appendix subheadings, pp. 1161,1168–1169. |

**Expression.** The 2009 earth-return correction, (1a)–(1b), p. 1065:

```math
Z'_{g_{ij}}=\frac{j\omega\mu_0}{\pi}
\int_0^\infty F(u)\cdot e^{-u(h_i+h_j)}\cos(y_{ij}u)\cdot du.
\qquad\text{(1a)}
```

```math
F(u)=\mu_1\frac{s_{12}+d_{12}e^{-2\alpha_1d}}
{s_{01}s_{12}+d_{01}d_{12}e^{-2\alpha_1d}}.
\qquad\text{(1b)}
```

Required definitions, 2009 appendix p. 1067:

```math
s_{mn}=(a_m\mu_n+a_n\mu_m),\qquad\text{(A.1)}
```

```math
d_{mn}=(a_m\mu_n-a_n\mu_m),\qquad\text{(A.2)}
```

```math
\begin{aligned}
a_m&=\sqrt{u^2+\gamma_m^2-\gamma_0^2} \\
\gamma_m^2&=j\omega\mu_m(\sigma_m+j\omega\varepsilon_m),
\qquad m,n=0,1,2 \\
\omega&=2\pi f.
\end{aligned}
```

The main equation really uses Greek ``\alpha_1`` while the appendix defines Latin ``a_m``. Their identification is **unresolved in the 2009 witness**; the later definition below is separately attributed, not silently inserted. ``F`` has length units inferred from these coefficients, giving ``Z'_g`` in ``\mathrm{\Omega/m}``. Layer 0 is air, layer 1 finite earth and layer 2 the lower half-space. All lengths are in metres under the source's SI examples.

For self, apply the source's substitutions ``y_{ij}\mapsto`` conductor ``i`` outer radius and ``h_j\mapsto h_i`` to (1). No extra absolute value or internal impedance is added. For the homogeneous-earth formula, section II.C prescribes replacing ``a_2`` with ``a_1`` and ``\gamma_2`` with ``\gamma_1``; the context is identical earth media. This statement does not independently resolve the Greek/Latin mismatch.

**Approximation.** Integral representation within the stated thin-conductor quasi-TEM model; no analytical truncation order or discarded spectral tail is specified. The 2010 derivation identifies an initially unknown longitudinal constant and replaces it by the air value while retaining only the TL mode, then applies a Bessel-transform identity. This later evidence clarifies the model family but is not backdated as a statement printed in 2009. Numerical integration is cited to an earlier Papagiannis et al. paper; no quadrature is implemented here.

**Limitations.** This original output is the earth correction alone. The 2009 ``\alpha_1/a_m`` ambiguity and missing root branches remain explicit. The 2010 external assembly is preserved separately, including its questionable appendix logarithmic identity, without using it to repair the original. No arbitrary-layer, buried, mixed, internal-conductor or insulation coverage is claimed.

**Reference.** [Papadopoulos2009](@cite), (1a)–(1b), p. 1065; appendix (A.1)–(A.2) and definitions, p. 1067. Separately inspected later witness: [Papadopoulos2010a](@cite), (2), (4a), (5a)–(5c), appendix A, pp. 1161–1162,1168–1170.

**Transcription source.** Original 2009 publication, PDF pages 1–4 inspected, with the principal equations checked on pages 2 and 4. The additional 2010 equation set below is verified against that later publication's own page images. These witnesses are not merged into an inferred corrected edition.

## Source transcription

The formula section reproduces the 2009 main equations and their original appendix notation. The following is a **separate 2010 source-ordered witness**, not another independent attribution for the same principal kernel.

2010 (2), p. 1161, supplies the matrix decomposition:

```math
\mathbf Z'(\omega)=\mathbf Z'_w+\mathbf Z'_e
=\mathbf Z'_w+\mathbf Z'_{pg}+\mathbf Z'_g.
\qquad\text{(2010:2)}
```

``\mathbf Z'_w`` is the diagonal internal conductor term, referred by the authors to Ametani's skin-effect expressions; it is not evaluated by this record. ``pg`` denotes perfectly conducting ground and ``g`` its imperfect-earth correction.

2010 (4a), p. 1162, before the longitudinal approximation:

```math
Z'_{e_{ij}}(\gamma_x)=\frac{j\omega\mu_0}{4\pi}
\int_0^\infty
\left[\frac{u}{a_0}
\left(e^{-a_0|h_j-h_i|}+e^{-a_0(h_j+h_i)}\cdot T'_1\right)\right]
\times\left[\int_{-\infty}^{\infty}
J_0\!\left(u\sqrt{x^2+y_{ij}^2}\right)e^{-\gamma_xx}\,dx\right]du.
\qquad\text{(2010:4a)}
```

The two integrations remain nested, with their original measures and domains. The main text uses unprimed ``a_0`` whereas the appendix's original spectral factor is ``a'_0``; this notation mismatch is not erased. Here ``J_0`` is the first-kind order-zero Bessel function, defined in the appendix p. 1169. No convergence condition for arbitrary complex ``\gamma_x`` is supplied with (4a); the record does not label this a solved dispersion relation.

2010 (5a)–(5c), after the explicitly prescribed ``\gamma_x=\gamma_0=jk_0=j\omega\sqrt{\mu_0\varepsilon_0}``:

```math
Z'_{e_{ij}}=Z'_{pg_{ij}}+Z'_{g_{ij}}
=\frac{j\omega\mu_0}{2\pi}\ln\frac{D_{ij}}{d_{ij}}
+\frac{j\omega\mu_0}{\pi}(P+jQ),\qquad\text{(2010:5a)}
```

```math
P+jQ=\int_0^\infty F_{strat}(\lambda)\cdot
e^{-\lambda(h_i+h_j)}\cos(y_{ij}\lambda)\cdot d\lambda,\qquad\text{(2010:5b)}
```

```math
F_{strat}(\lambda)=\mu_1
\frac{s_{12}+d_{12}e^{-2\alpha_1d}}
{s_{01}s_{12}+d_{01}d_{12}e^{-2\alpha_1d}}.\qquad\text{(2010:5c)}
```

``P+jQ`` here is a dimensionless integral correction, **not** the potential-coefficient matrix of the admittance record. No logarithmic term is retroactively appended to the 2009 correction-only output.

The original spectral definitions below 2010 (A.3), p. 1169, and the parent reflection coefficients (A.6a), (A.7a)–(A.7b), (A.7d)–(A.7e) are:

```math
\begin{aligned}
\gamma_k^2&=j\omega\mu_k(\sigma_k+j\omega\varepsilon_k) \\
a'_k&=\sqrt{u^2+\gamma_k^2} \\
k&=0,1,2,
\end{aligned}
```

```math
T'_1=\frac{\Delta'_1}{\Delta'}=
\frac{d'_{01}s'_{12}+s'_{01}d'_{12}e^{-2a'_1d}}
{s'_{01}s'_{12}+d'_{01}d'_{12}e^{-2a'_1d}},\qquad\text{(2010:A.6a)}
```

```math
\Delta'=s'_{01}s'_{12}+d'_{01}d'_{12}e^{-2a'_1d},\qquad\text{(2010:A.7a)}
```

```math
\Delta'_1=d'_{01}s'_{12}+s'_{01}d'_{12}e^{-2a'_1d},\qquad\text{(2010:A.7b)}
```

```math
\begin{aligned}
s'_{mn}&=(a'_m\mu_n+a'_n\mu_m) \\
d'_{mn}&=(a'_m\mu_n-a'_n\mu_m),\qquad m,n=0,1,2.
\end{aligned}\qquad\text{(2010:A.7d--e)}
```

2010 (A.12), p. 1169, is the source-provided longitudinal-transform identity, retained as printed:

```math
\int_{-\infty}^{\infty}
J_0\!\left(u\sqrt{x^2+y_{ij}^2}\right)e^{-jk_0x}\,dx
=\begin{cases}
0,&u<k_0,\\
2\dfrac{\cos\!\left(y_{ij}\sqrt{u^2-k_0^2}\right)}
{\sqrt{u^2-k_0^2}},&u>k_0.
\end{cases}\qquad\text{(2010:A.12)}
```

No value at ``u=k_0`` is printed. The immediately following text prescribes:

```math
\begin{aligned}
u^2-k_0^2&=\lambda^2 \\
a'_k\longmapsto\alpha_k&=\sqrt{\lambda^2+\gamma_k^2+k_0^2} \\
k&=0,1,2,\qquad T'_1\longmapsto T_1.
\end{aligned}
```

Thus the final ``s_{mn},d_{mn}`` are the corresponding coefficients after the author's substitution of ``\alpha`` for the original ``a'``; their dependency rendering is:

```math
\begin{aligned}
s_{mn}&=(\alpha_m\mu_n+\alpha_n\mu_m) \\
d_{mn}&=(\alpha_m\mu_n-\alpha_n\mu_m).
\end{aligned}
```

This is the **2010 transform prescription**, not a replacement definition for the printed 2009 symbols. Appendix (A.13), p. 1170, prints the following logarithmic identity with a **plus** and primed factors even after the transform:

```math
\begin{aligned}
\int_0^\infty
\left(\frac{e^{-a'_0|h_j-h_i|}}{a'_0}
+\frac{e^{-a'_0(h_i+h_j)}}{a'_0}\right)
\cos(y_{ij}\lambda) \\
d\lambda&=\ln\frac{D_{ij}}{d_{ij}}.
\end{aligned}\qquad\text{(2010:A.13)}
```

Immediately below it:

```math
\begin{aligned}
D_{ij}&=\sqrt{y_{ij}^2+(h_i+h_j)^2} \\
d_{ij}&=\sqrt{y_{ij}^2+(h_i-h_j)^2}.
\end{aligned}
```

The plus in (A.13) and its post-transform ``a'_0`` are retained as a suspected published defect/notation conflict, not accepted as a repaired proof of (5a). No singularity claim is inferred from an isolated term. The source's self prescription applies to (5) as a whole. Its homogeneous reduction sets the electromagnetic properties of the two earth layers equal, explicitly ``\gamma_2=\gamma_1,a_2=a_1`` (p. 1162); the ``a`` notation is printed there too.

## Notation map

No normalization is performed inside either witness. The table distinguishes symbols that look similar but are not automatically interchangeable. Units are equation-implied SI unless stated otherwise.

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z'_{g_{ij}},Z'_{pg_{ij}},Z'_{e_{ij}},\mathbf Z'_w,\mathbf Z'`` | unchanged | Earth correction; perfect-ground external term; full external term; internal and full series matrices | ``\mathrm{\Omega/m}``; full assembly only explicitly printed in 2010 |
| ``F(u),F_{strat}(\lambda)`` | unchanged | 2009 and 2010 impedance kernels | Length; their exponential weight stays outside the kernel |
| ``P,Q`` in 2010 (5b) | unchanged | Real/imaginary names of the integral correction | Dimensionless; not potential coefficients |
| ``u`` (2009), ``\lambda`` (2010 final) | unchanged | Final transverse integration arguments | ``\mathrm{m^{-1}}``; only comparison correspondence, no renaming inside equations |
| ``u`` (2010 parent) | unchanged | Original Hankel-transform argument | ``\mathrm{m^{-1}}``; not identified with final 2009 ``u`` |
| ``a_m,\alpha_1,a'_1`` (2009) | unchanged | Defined Latin spectral factor, main Greek exponent factor, appendix primed occurrence | ``\mathrm{m^{-1}}``; identification gaps retained |
| ``a'_k,\alpha_k,a_0`` (2010) | unchanged | Original, transformed, and main-parent spectral notation | ``\mathrm{m^{-1}}``; no explicit root-branch prescription |
| ``s_{mn},d_{mn},s'_{mn},d'_{mn}`` | unchanged | Magnetic interface combinations in each witness | Defined products of permeability and spectral factor |
| ``T'_1,\Delta',\Delta'_1`` | unchanged | Parent reflection quotient and denominator/numerator | 2010 appendix definitions; not derivatives |
| ``\gamma_m,\gamma_k`` | unchanged | Bulk-medium constants | ``\mathrm{m^{-1}}``; indices 0 air,1 finite earth,2 lower earth |
| ``\gamma_x,k_0`` | unchanged | 2010 longitudinal constant and free-space wavenumber | ``\mathrm{m^{-1}}``; parent unknown, final prescribed; ``e^{-\gamma_xx}`` |
| ``\mu_m,\varepsilon_m,\sigma_m`` | unchanged | Layer permeability, permittivity, conductivity | ``\mathrm{H/m},\mathrm{F/m},\mathrm{S/m}``; scalar |
| ``h_i,h_j,y_{ij},d,D_{ij},d_{ij},x,z`` | unchanged | Heights, separation, layer depth, image/direct distances, longitudinal/upward appendix coordinates | Metres; ``d`` is not ``d_{ij}`` |
| ``J_0,j,\omega,f,c,h`` | unchanged | First-kind Bessel function, imaginary unit, angular frequency, frequency, light speed, height in TL criterion | Bessel order zero; ``\mathrm{s^{-1}},\mathrm{Hz},\mathrm{m/s},\mathrm m`` as applicable |
| ``\alpha_{other\ model},\alpha_{proposed}`` | unchanged | Modal attenuation constants in the later source's comparison measure | ``\mathrm{neper/m}`` in its plots; not the spectral ``\alpha_k`` |

## Evidence and approximation sources

2009 I–II.A supplies the quasi-TEM restriction and refers the Hertz-vector methodology to Sunde [2]. Its appendix contains no explicit longitudinal exponential or root branch. The 2010 section 3.1 gives the missing **later-author** longitudinal explanation: the parent (4) has unknown ``\gamma_x``; an additional surface-attached/fast-wave mode is discussed; the retained TL mode is represented by the prescribed air propagation constant. Appendix (A.12) then removes the inner longitudinal integral. There is no source-stated small-parameter order, series truncation, quadrature tolerance or error theorem to attach to this reduction.

The 2009 (1b) and 2010 (5c) main quotients match term by term, while their spectral-variable/definition witnesses differ as recorded. No claim of fully resolved equivalence across the inconsistent tokens is made. The 2010 logarithmic/external assembly and parent integral are useful additional published evidence, not counted as a second newly invented earth-correction kernel. The [potential-coefficient record](../../../external-admittance/2009/two-layer-earth-overhead-potential-correction/Papadopoulos2009.md) retains the additional displacement kernel and distinguishes potential from admittance.

Author-stated reductions: 2009 II.B identifies ``\mu_i=\mu_0,\varepsilon_0=0,\varepsilon_i\ne0`` with Sunde's two-layer impedance, and ``\mu_i\ne\mu_0,\varepsilon_i\ne\varepsilon_0`` with its reference [4]. These are attributed identifications, not independently derived kernels or proof of the cited source's priority. The 2009 [4] reference names Nakagawa's **1981 admittance-correction** article; the 2010 comparison instead cites Nakagawa–Ametani–Iwamoto's **1973 stratified-earth impedance** article as [14]. The 2009 homogeneous paragraph also calls the result Kikuchi's model but points to [2], which is Sunde; its Kikuchi reference is [3].

Accuracy evidence is model comparison, not mathematical certification. 2009 IV–V tests the 150 kV horizontal three-phase arrangement (10 m height, 6.5 m adjacent spacing), layer resistivity ratios 10,5,0.1,0.2 with lower resistivity 10,100,1000 ``\mathrm{\Omega\,m}``, relative permittivities 1–20, upper depths 5–20 m and unit relative permeability. The 2010 study states 50 Hz–10 MHz and includes ground-mode and transient comparisons. Its Eq. (8) defines percentage difference as ``|\alpha_{other\ model}-\alpha_{proposed}|/|\alpha_{proposed}|\times100``; these modal ``\alpha`` values are attenuation constants, not the spectral ``\alpha_k``. No quoted result is treated as a universal error bound or as an independent full-wave reference.

## Limitations and discrepancies

1. **Printed 2009 notation gap:** Greek ``\alpha_1`` in (1b) is not explicitly connected to the appendix's Latin ``a_m``. The related 2009 appendix (A.5), used by the admittance denominator, also prints ``a'_1`` without defining that prime. The later paper supplies a distinct transform definition, not an erratum.
2. **Printed 2010 notation/identity issues:** (4a) has unprimed ``a_0`` while the original appendix factor is primed; (A.13) retains primed factors and a plus between the two exponentials. The image was checked; neither feature is an OCR repair opportunity. The exact scope of the claimed identity needs author clarification.
3. **Branch and parent-domain gaps:** neither inspected publication gives an explicit square-root branch for the final factors. The 2010 parent integral does not state convergence conditions for arbitrary complex ``\gamma_x``; the final prescribed model is kept separate.
4. **Output restriction:** the original correction is not total external impedance; the latter is explicitly supplied only by the separately transcribed 2010 (5a). Neither expression includes the conductor's internal skin/proximity term.
6. **Earlier-source comparison limits:** Sunde and the Nakagawa sources remain separately attributed dependencies/candidates. The independently inspected [Ametani–Nagaoka–Koide snow record](../../../external-admittance/2000/snow-layer-overhead-potential-coefficient/Ametani2001.md) now retains the 2000 Japanese original and 2001 published translation; that paper develops a potential coefficient but uses Carson for series impedance after neglecting snow effects on that quantity. No new snow-modified impedance kernel is inferred from its admittance attribution. Inspecting these sources does not establish first-ever priority, a universal range, or an arbitrary-layer recurrence. No separate published erratum was inspected.

## Numerical interpretation

Numerical evaluation uses the final two-layer magnetic kernel with the prescribed air propagation reference. The corresponding transverse constants are used consistently in the interface coefficients and finite-layer exponential. The printed self-radius prescription is retained. Independent source-kernel comparisons include unequal permeabilities and negative-frequency conjugacy.
