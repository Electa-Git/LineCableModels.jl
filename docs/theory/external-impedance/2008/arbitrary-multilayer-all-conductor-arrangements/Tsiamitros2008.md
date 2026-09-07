# Tsiamitros–Papagiannis–Dokopoulos arbitrary-multilayer earth-return impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinite straight single-core conductors parallel to the ``x`` axis, represented externally by filamentary horizontal dipoles. Insulation thickness is assumed negligible relative to conductor diameter for this earth term; insulation impedance is a separate addend (8). Self replaces ``y_{ij}`` by the outermost cable radius and sets ``m=l,h_2=h_1``. |
| Calculated quantities | Per-unit-length self and mutual impedance for overhead, buried same-layer/cross-layer, and mixed conductor placements in an arbitrary number of horizontal earth layers |
| Earth structure | ``n`` horizontal earth layers below air; finite layers ``1,\ldots,n-1`` and a semi-infinite ``n``th layer. |
| Model and approximation | Not an analytical low- or high-frequency approximation within the quasi-TEM, filamentary, horizontal-layer model. Equation (21) is the direct field-equation result. The recursions are exact algebraic interface assemblies in that model. The paper separately prescribes numerical quadrature; no quadrature rule is substituted into the equation here. |
| Main source | D. A. Tsiamitros, G. K. Papagiannis, and P. S. Dokopoulos (2008) |
| Citation key(s) | `:Tsiamitros2008` |
| Evidence status | Original publication page images checked for the main kernel, both recursions, terminal cases, and self prescription. |

**Description.** Per-unit-length earth-return impedance between infinitely long parallel conductors placed above, within the same layer of, or within different layers of a horizontally stratified earth. The spectral kernel combines upward and downward interface recursions and includes layer conductivity, permittivity, and permeability.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | No independent ``\Gamma`` is retained. The paper assumes quasi-TEM propagation and infinite uniform conductors, integrates the horizontal dipole along the conductor, and omits end-gradient terms. | Stated — p. 2393 and derivation (19)–(21), pp. 2394–2395. |
| Air propagation constant ``γ_air`` | Region 0: ``\bar\gamma_0^2=j\omega\mu_0(\sigma_0+j\omega\varepsilon_0)`` with ``\sigma_0=0`` and free-space ``\mu_0,\varepsilon_0``; ``\bar\alpha_0=\sqrt{u^2+\bar\gamma_0^2}``. | Stated — p. 2392 and definition below (13), p. 2394. |
| Earth propagation constant ``γ_earth`` | Each earth layer ``i`` uses ``\bar\gamma_i^2=j\omega\mu_i(\sigma_i+j\omega\varepsilon_i)`` and ``\bar\alpha_i=\sqrt{u^2+\bar\gamma_i^2}``. | Stated — definition below (13), p. 2394. |
| Earth permittivity and displacement current | Retained independently in every layer through ``j\omega\varepsilon_i``. | Stated — pp. 2392 and 2394. |
| Range of validity | The authors state that (21) is neither a low- nor a high-frequency approximation. Physical validity remains limited by quasi-TEM propagation, infinite uniform parallel conductors, filamentary external source, and horizontally stratified linear media. No universal numerical frequency bound is stated in Part I. | Stated — pp. 2393 and 2395. |
| Earth permeability ``μ_earth`` | Arbitrary scalar layer values ``\mu_i`` are retained. | Stated — p. 2392 and (15)–(21), pp. 2394–2395. |
| Arrangement | Overhead/overhead, buried/buried in the same or different layers, and mixed overhead/buried arrangements; self and mutual. The displayed general construction orders the layer indices as in Fig. 1 and uses reciprocity for the reversed excitation. | Stated — abstract, Figs. 1–2, (14), (20)–(21), and pp. 2392–2395. |
| Earth structure | ``n`` horizontal earth layers below air; finite layers ``1,\ldots,n-1`` and a semi-infinite ``n``th layer. | Stated — Fig. 1 and §II-A, pp. 2392–2393. |
| Conductor and insulation geometry | Infinite straight single-core conductors parallel to the ``x`` axis, represented externally by filamentary horizontal dipoles. Insulation thickness is assumed negligible relative to conductor diameter for this earth term; insulation impedance is a separate addend (8). Self replaces ``y_{ij}`` by the outermost cable radius and sets ``m=l,h_2=h_1``. | Stated — p. 2393 and self prescription below (21), p. 2395. |
| Constitutive and field assumptions | Linear, bilateral, homogeneous isotropic material within each layer; quasi-TEM; end effects absent. Skin and proximity are not contained in (21) and require separate internal impedance terms. | Stated — pp. 2393–2394. |
| Conventions | Time-harmonic ``j\omega`` convention as printed; conductors run along ``x``; ``y_{ij}`` is horizontal transverse separation; layer depths/thickness coordinates follow Figs. 1–2. Output is per unit length. | Stated — Figs. 1–2 and §§II–III, pp. 2392–2395. |

**Expression.** General mutual earth-return impedance for conductor ``i`` in layer ``m`` and conductor ``j`` in layer ``l`` under the source ordering, equation (21), printed p. 2395.

```math
\begin{aligned}
\bar Z_{ij}={}&\frac{j\omega\mu_m}{2\pi}\int_0^\infty
\frac{\cos(uy_{ij})}{\bar\alpha_m}
\Bigg\{
2^{m-l}
\frac{
(\mu_1\mu_2\cdots\mu_{m-1})
(\bar\alpha_1\bar\alpha_2\cdots\bar\alpha_{m-1}\bar\alpha_m)
}{
(\mu_1\mu_2\cdots\mu_{l-1})
(\bar\alpha_1\bar\alpha_2\cdots\bar\alpha_{l-1}\bar\alpha_l)
}\\
&\qquad\times
\left(e^{-\bar\alpha_l d_l}e^{-\bar\alpha_{l+1}d_{l+1}}\cdots
e^{-\bar\alpha_m d_m}\right)
\frac{\bar F_1\bar F_2}{\overline{DTD}_0}
\Bigg\}\,du,
\end{aligned}
\qquad\text{(21)}
```

```math
\bar F_1=
\overline{DTD}_m e^{\bar\alpha_m(d_m-h_1)}
+\overline{DTN}_m e^{-\bar\alpha_m(d_m-h_1)},
```

```math
\bar F_2=
\overline{TDD}_{l-1}e^{\bar\alpha_lh_2}
+\overline{TDN}_{l-1}e^{-\bar\alpha_lh_2}.
```

The spectral and constitutive definitions are

```math
\begin{aligned}
\bar\alpha_i&=\sqrt{u^2+\bar\gamma_i^2} \\
\bar\gamma_i^2&=j\omega\mu_i(\sigma_i+j\omega\varepsilon_i) \\
i&=0,1,\ldots,n.
\end{aligned}
```

The downward-to-top interface recursion, source (15)–(18), is

```math
\overline{DTN}_i=
(\mu_{i+1}\bar\alpha_i-\mu_i\bar\alpha_{i+1})\overline{DTD}_{i+1}
+(\mu_{i+1}\bar\alpha_i+\mu_i\bar\alpha_{i+1})
\overline{DTN}_{i+1}e^{-2\bar\alpha_{i+1}d_{i+1}},
\qquad\text{(15)}
```

```math
\overline{DTD}_i=
(\mu_{i+1}\bar\alpha_i+\mu_i\bar\alpha_{i+1})\overline{DTD}_{i+1}
+(\mu_{i+1}\bar\alpha_i-\mu_i\bar\alpha_{i+1})
\overline{DTN}_{i+1}e^{-2\bar\alpha_{i+1}d_{i+1}},
\qquad\text{(16)}
```

```math
\begin{aligned}
\overline{DTN}_n&=0 \\
\overline{DTD}_n&=1.
\end{aligned}\qquad\text{(17--18)}
```

The top-to-down recursion, source (A.15)–(A.18), is

```math
\begin{aligned}
\overline{TDD}_{-1}&=1 \\
\overline{TDN}_{-1}&=0,
\end{aligned}\qquad\text{(A.15--A.16)}
```

```math
\overline{TDD}_{l-1}=
(\mu_{l-1}\bar\alpha_l+\mu_l\bar\alpha_{l-1})\overline{TDD}_{l-2}
+(\mu_{l-1}\bar\alpha_l-\mu_l\bar\alpha_{l-1})
\overline{TDN}_{l-2}e^{-2\bar\alpha_{l-1}d_{l-1}},
\qquad\text{(A.17)}
```

```math
\overline{TDN}_{l-1}=
(\mu_{l-1}\bar\alpha_l-\mu_l\bar\alpha_{l-1})\overline{TDD}_{l-2}
+(\mu_{l-1}\bar\alpha_l+\mu_l\bar\alpha_{l-1})
\overline{TDN}_{l-2}e^{-2\bar\alpha_{l-1}d_{l-1}}.
\qquad\text{(A.18)}
```

For self impedance, the source sets ``m=l``, replaces ``y_{ij}`` by the cable outermost radius ``r_{ii}``, and replaces ``h_2`` by ``h_1``.

**Approximation.** Not an analytical low- or high-frequency approximation within the quasi-TEM, filamentary, horizontal-layer model. Equation (21) is the direct field-equation result. The recursions are exact algebraic interface assemblies in that model. The paper separately prescribes numerical quadrature; no quadrature rule is substituted into the equation here.

**Limitations.** The kernel does not contain conductor skin/proximity, a finite-thickness insulation field, conductor pitching, lateral earth variation, end effects, or higher propagation modes. The source does not print a square-root branch statement next to ``\bar\alpha_i``; decay selection is therefore unresolved rather than invented. The notation uses ``d_i`` in layer exponentials and local conductor coordinates ``h_1,h_2`` as drawn; those coordinates must not be replaced by unsigned global depths without reproducing the source transformations.

**Reference.** [Tsiamitros2008](@cite), equations (15)–(18), (21), and (A.15)–(A.18), printed pp. 2394–2395 and 2399 (PDF pages 3–4 and 8), with specializations on pp. 2395–2398.

**Transcription source.** Original IEEE Part I page images. The product bounds, ``2^{m-l}``, exponential chain, permeability/vertical-root factors, ``\bar F_1,\bar F_2`` signs, recursion coefficients and terminal values were visually checked. PDF text extraction omits most equation bodies and was navigation only.

## Source transcription

The displayed equations retain every source overbar and the literal DT/TD letter order. The letters are source mnemonics: DT means construction from lower layers toward the top; TD means construction from the top toward lower layers; final D/N indicates the role normally played as denominator/numerator in the derivation.

The paper verifies (21) by reduction, without defining those reductions as new contributions:

- overhead conductors above homogeneous earth: (28)–(29), identified with Carson;
- underground cables in homogeneous earth: (37)–(38), identified with Pollaczek;
- mixed overhead/underground in homogeneous earth: (45), identified with Pollaczek;
- overhead above two-layer earth: (50), identified with Sunde/Wedepohl;
- buried and mixed cases in two-layer earth: (56) and (63), identified with the authors' earlier papers.

These limiting witnesses support the arrangement scope of (21) and are not counted as separate 2008 original formulas.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\bar Z_{ij}`` | unchanged | External earth-return impedance between conductor pair | ``\Omega/\mathrm m`` |
| ``u`` | unchanged | Transverse spectral variable | ``\mathrm m^{-1}`` |
| ``\bar\gamma_i`` | unchanged | Bulk medium constant in layer ``i`` | ``\mathrm m^{-1}`` |
| ``\bar\alpha_i`` | unchanged | Vertical spectral root | ``\mathrm m^{-1}`` |
| ``\mu_i,\varepsilon_i,\sigma_i`` | unchanged | Permeability, permittivity and conductivity of region ``i`` | SI; region 0 is air |
| ``n,m,l`` | unchanged | Final earth-layer index and conductor layer indices | source construction uses the ordering shown in Fig. 1 |
| ``d_i`` | unchanged | Source layer depth/thickness coordinate in exponential propagation factors | metres; geometry follows Figs. 1–2 |
| ``h_1,h_2`` | unchanged | Local vertical coordinates of conductors ``i,j`` in their layers | metres |
| ``y_{ij}`` | unchanged | Horizontal separation normal to the conductor axes | metres; self replaced by outer radius |
| ``\overline{DTD},\overline{DTN}`` | unchanged | Lower-to-upper interface-recursion terms | source-defined |
| ``\overline{TDD},\overline{TDN}`` | unchanged | Upper-to-lower interface-recursion terms | source-defined |
| ``\bar F_1,\bar F_2`` | unchanged | Source and observation layer factors | source-defined |

No notation was renamed.

## Evidence and approximation sources

1. Air and layer material definitions, final semi-infinite layer and general arrangement: abstract and §II-A, pp. 2392–2393.
2. Quasi-TEM, filamentary external field, insulation-thickness reduction and infinite-length assumption: p. 2393.
3. Spectral/root definitions: below (13), p. 2394.
4. Down-to-top recursion and mixed construction: (14)–(18), p. 2394.
5. General same-/cross-layer impedance, self substitution and non-approximation statement: (20)–(21), p. 2395.
6. Limiting arrangements: §§IV-A–F, pp. 2395–2398.
7. Top-to-down recursion: (A.15)–(A.18), p. 2399.

## Limitations and discrepancies

- The source DOI is printed as `10.1109/TPWRS.2008.923816` despite publication in *IEEE Transactions on Power Delivery*. The printed DOI is retained as bibliographic identity evidence, not altered.
- The source's claim that (21) can be combined with nonlinear time-domain phenomena does not make its materials nonlinear; (21) itself is derived from linear layer field equations.
- Square-root decay branches and treatment of a conductor in the terminal infinite layer require source-consistent limiting evaluation; this record does not insert a stabilized reformulation.

## Numerical interpretation

The implementation uses the printed upward and downward interface recursions for physical layer intervals. Square roots follow the decaying or outgoing branch, including negative-frequency conjugacy. Independent magnetic boundary-condition solves check unequal permeabilities and source and observation points in different layers. This selection supplies impedance only.
