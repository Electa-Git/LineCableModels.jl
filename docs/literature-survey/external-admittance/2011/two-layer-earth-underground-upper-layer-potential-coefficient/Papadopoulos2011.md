# Papadopoulos–Tsiamitros–Papagiannis upper-layer buried-cable earth potential coefficient

## Identity and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite parallel single-core cable axes; outer radius enters self substitution. The earth coefficient is separate from dielectric-insulation and conductor contributions; the source references their assembly to [9,15]. |
| Calculated quantities | Mutual earth-return potential coefficient and source-printed admittance relation; self earth term by radius/depth substitutions |
| Earth structure | Air half-space, finite upper earth layer of thickness ``d``, infinite-depth lower earth half-space. General different-layer/mixed/arbitrary-layer arrangements are not supplied by this final kernel. |
| Model and approximation | The source replaces the unknown longitudinal propagation constant by the dielectric value of the upper earth layer and applies transform (5). ``G`` is retained; no spectral truncation or tail approximation is specified. |
| Main source | T. A. Papadopoulos, D. A. Tsiamitros, and G. K. Papagiannis (2011), DOI `10.1049/iet-gtd.2010.0228` |
| Citation key(s) | `:Papadopoulos2011` |
| Evidence status | Original PDF equations checked; the transformed-factor index, root branch, and scalar/matrix interpretation remain unresolved. |

**Description.** Per-unit-length earth-related potential coefficient between parallel single-core cable axes in the finite upper layer of a two-layer earth, with air above. The source adds four radial-displacement-current kernel contributions to its impedance kernel and prints a corresponding admittance relation. Self interaction follows the source's geometric substitutions.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Source ``\gamma_x`` with ``e^{-\gamma_x x}``; final model sets ``\gamma_x=j\omega\sqrt{\mu_1\varepsilon_1}``, not the lossy bulk ``\gamma_1`` and not zero. Actual line propagation is calculated later from the resulting per-unit-length parameters. | Stated — section 2 before (5), p. 162; parent (2). |
| Air propagation constant ``γ_air`` | Indexed definition gives ``\gamma_0^2=j\omega\mu_0(\sigma_0+j\omega\varepsilon_0)``. Air is labelled ``\mu_0,\varepsilon_0``; no explicit ``\sigma_0`` value is supplied in the inspected appendix definition. Air enters ``G_3,G_4`` and interface coefficients; its propagation is not set to zero. | Stated/equation-implied — Fig. 1, (7f)–(7g), appendix below (10), pp. 162–163,171. |
| Earth propagation constant ``γ_earth`` | ``\gamma_k^2=j\omega\mu_k(\sigma_k+j\omega\varepsilon_k)``, earth ``k=1,2``. These lossy bulk terms are retained separately from imposed ``\gamma_x``. Printed transformation ``a_m=\sqrt{\lambda^2+\gamma_k^2+k_x^2}`` has unresolved ``m/k`` indexing. | Stated — p. 162 below (5), p. 171 below (10). |
| Earth permittivity and displacement current | ``\varepsilon_1,\varepsilon_2`` retained; ``G=2a_1(G_1+G_2+G_3+G_4)`` represents radial displacement current. Source says ignoring ``G`` results in propagation constant ``\gamma_1``. No separate complex-permittivity/loss-conductivity decomposition is specified. | Stated/equation-implied — (7) and section 3, p. 163. |
| Range of validity | Quasi-TEM TL mode, excluding antenna modes. The stated up-to-10-MHz TL approximation is conditional on typical earth resistivity and is cited to [19], not a universal admittance-kernel bound. The source/target must both be in the finite first earth layer. No expansion order or universal error estimate is supplied. | Stated — section 3 and Fig. 1, pp. 162–163. |
| Earth permeability ``μ_earth`` | Layer-dependent scalar ``\mu_1,\mu_2`` retained; unity relative permeability is only the test choice in section 4. | Equation-implied — (7d)–(7g), (19)–(22); stated example choice p. 164. |
| Arrangement | Mutual buried cable-axis interaction; self by replacing ``y_{ij}`` with outermost cable radius and ``h_j`` with ``h_i``. Either depth ordering is explicitly included by the appendix's statement about (6)–(7). | Stated — p. 163 after (7), p. 171 final paragraph. |
| Earth structure | Air half-space, finite upper earth layer of thickness ``d``, infinite-depth lower earth half-space. General different-layer/mixed/arbitrary-layer arrangements are not supplied by this final kernel. | Stated — section 2/Fig. 1 and proposed extensions in section 3. |
| Conductor and insulation geometry | Infinite parallel single-core cable axes; outer radius enters self substitution. The earth coefficient is separate from dielectric-insulation and conductor contributions; the source references their assembly to [9,15]. | Stated — sections 2–3, pp. 162–163. |
| Constitutive and field assumptions | Spatially uniform scalar electromagnetic properties within each layer; linear Hertz-vector superposition with the prescribed quasi-TEM longitudinal approximation. This earth expression supplies no new internal skin/proximity or insulation constitutive model. | Equation-implied — (2), appendix scalar definitions; stated — pp. 162–163. |
| Conventions | ``j`` imaginary, ``\omega=2\pi f``, longitudinal ``e^{-\gamma_x x}``. Positive ``j\omega`` constitutive factors imply the corresponding harmonic differentiation convention, but no explicit time exponential is stated here. Positive depths measured downward from air interface; appendix ``z`` is compared with ``d-h_i``. ``P`` is a potential coefficient; ``Y'`` is per-length admittance as printed. | Stated/equation-implied — Fig. 1, (2), (7a)–(7b), appendix (8)–(10) and definitions. |

**Expression.** Source-printed admittance relation and mutual earth potential coefficient, (7a)–(7g), with the complete shared kernel (6b). The inverse in (7a) is retained in its printed scalar/subscripted form; this record does not replace it with a matrix inverse.

```math
Y'_{e_{ij}}=j\omega P_{e_{ij}}^{-1}.
\qquad\text{(7a)}
```

```math
P_{e_{ij}}=\frac{j\omega}{2\pi(\sigma_1+j\omega\varepsilon_1)}
\int_0^{+\infty}[F(\lambda)+G(\lambda)]\cos(y_{ij}\lambda)\,d\lambda.
\qquad\text{(7b)}
```

```math
G(\lambda)=2a_1[G_1(\lambda)+G_2(\lambda)+G_3(\lambda)+G_4(\lambda)].
\qquad\text{(7c)}
```

```math
G_1(\lambda)=
\frac{\mu_1\mu_2(\gamma_1^2-\gamma_2^2)
 (S_{10}A_{10}e^{-a_1(2d-h_1-h_2)}-D_{10}A_{10}e^{-a_1(2d+h_1-h_2)})}
 {(A_{10}A_{12}-\Delta_{10}\Delta_{12}e^{-2a_1d})
  (S_{10}S_{21}+D_{10}D_{21}e^{-2a_1d})}.
\qquad\text{(7d)}
```

```math
G_2(\lambda)=
\frac{\mu_1\mu_2(\gamma_1^2-\gamma_2^2)
 (S_{10}\Delta_{10}e^{-a_1(2d+h_2-h_1)}-D_{10}\Delta_{10}e^{-a_1(2d+h_1+h_2)})}
 {(A_{10}A_{12}-\Delta_{10}\Delta_{12}e^{-2a_1d})
  (S_{10}S_{21}+D_{10}D_{21}e^{-2a_1d})}.
\qquad\text{(7e)}
```

```math
G_3(\lambda)=
\frac{\mu_1\mu_0(\gamma_1^2-\gamma_0^2)
 (S_{21}\Delta_{12}e^{-a_1(2d+h_1-h_2)}+D_{21}\Delta_{12}e^{-a_1(4d-h_1-h_2)})}
 {(A_{10}A_{12}-\Delta_{10}\Delta_{12}e^{-2a_1d})
  (S_{10}S_{21}+D_{10}D_{21}e^{-2a_1d})}.
\qquad\text{(7f)}
```

```math
G_4(\lambda)=
\frac{\mu_1\mu_0(\gamma_1^2-\gamma_0^2)
 (S_{21}A_{12}e^{-a_1(h_1+h_2)}+D_{21}A_{12}e^{-a_1(2d+h_2-h_1)})}
 {(A_{10}A_{12}-\Delta_{10}\Delta_{12}e^{-2a_1d})
  (S_{10}S_{21}+D_{10}D_{21}e^{-2a_1d})}.
\qquad\text{(7g)}
```

The shared ``F`` is repeated without algebraic regrouping so the record is independently readable:

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

The transformation following (5) prescribes unprimed versions of appendix (19)–(22):

```math
S_{mn}=(\mu_n a_m+\mu_m a_n),\qquad
D_{mn}=(\mu_m a_n-\mu_n a_m),
```

```math
A_{mn}=(a_n\gamma_m^2\mu_n+a_m\gamma_n^2\mu_m),\qquad
\Delta_{mn}=(a_n\gamma_m^2\mu_n-a_m\gamma_n^2\mu_m),
\qquad m,n=0,1,2.
```

```math
a_m=\sqrt{\lambda^2+\gamma_k^2+k_x^2},\qquad
\gamma_k^2=j\omega\mu_k(\sigma_k+j\omega\varepsilon_k),
\qquad k=0,1,2,
```

```math
\gamma_x=j\omega\sqrt{\mu_1\varepsilon_1},\qquad
\gamma_x=jk_x,\qquad \omega=2\pi f.
```

The ``a_m/\gamma_k`` mismatch is printed in the source and deliberately not changed. ``\gamma_x=jk_x`` identifies (5)'s exponent with the imposed exponent in the immediately preceding prose. ``\Delta`` and ``D`` are distinct coefficient families and cannot be interchanged. ``h_1,h_2`` are the source's final-kernel labels for the two cable depths; prose uses ``h_i,h_j``. ``P`` has potential-coefficient meaning, with units inferred from (7a) rather than separately printed beside it; ``Y'`` is admittance per length. For self, the source prescribes ``y_{ij}`` replaced by cable ``i``'s outermost radius and ``h_j`` by ``h_i``. No insulation coefficient or full cable-admittance matrix is added to (7).

**Approximation.** Parent expression (2), with Hertz amplitude (4), contains an initially unknown longitudinal propagation constant. The author approximates it by the purely imaginary dielectric value of the upper earth layer, then evaluates the longitudinal Bessel transform by (5). No spectral-series truncation, expansion parameter, retained order, or tail approximation is specified. ``G`` is retained in this formulation; the source's separate observation about omitting it is not applied here. Integral form does not make the prescribed quasi-TEM model full-wave.

**Limitations.** Both cable axes must lie in the upper finite layer. The potential coefficient is not itself an admittance or the total cable shunt matrix. The printed scalar inverse in (7a), transformed-factor index mismatch, missing root branch, and differing appendix/final coefficient indices remain visible rather than repaired.

**Reference.** [Papadopoulos2011](@cite), (6b), (7a)–(7g), printed p. 163; propagation prescription and transform (5), p. 162; appendix unnumbered definitions and (19)–(22), p. 171.

**Transcription source.** Original publication, visually checked against PDF pages 2–3 and 11 (indices 1–2 and 10). ``D`` versus ``\Delta``, each exponent, numerator sign, layer index, and both denominator factors of all four ``G`` terms were inspected in the page image. No converted Markdown or LCM implementation supplies missing mathematics.

## Source transcription

The formula section preserves (7a)–(7g) in their source equation order; (6b) is repeated afterward solely to provide the required shared dependency. No ``G`` terms or separate denominator factors are combined. The following original primed coefficient definitions, appendix (19)–(22), establish the symbol distinctions before the author's unpriming transformation:

```math
S'_{mn}=(\mu_n a'_m+\mu_m a'_n),\qquad\text{(19)}
```

```math
D'_{mn}=(\mu_m a'_n-\mu_n a'_m),\qquad\text{(20)}
```

```math
A'_{mn}=(a'_n\gamma_m^2\mu_n+a'_m\gamma_n^2\mu_m),\qquad\text{(21)}
```

```math
\Delta'_{mn}=(a'_n\gamma_m^2\mu_n-a'_m\gamma_n^2\mu_m),\qquad m,n=0,1,2.
\qquad\text{(22)}
```

Immediately below appendix (10), p. 171:

```math
a'_k=\sqrt{u^2+\gamma_k^2},\qquad
\gamma_k^2=j\omega\mu_k(\sigma_k+j\omega\varepsilon_k),\qquad k=0,1,2.
```

The transformed factor quoted in the formula section is the **separate printed witness** below (5), p. 162. The original transform identity is:

```math
\int_{-\infty}^{\infty}
 J_0\!\left(u\sqrt{x^2+y_{ij}^2}\right)e^{-jk_xx}\,dx
=\begin{cases}
0,&u<k_x,\\
2\dfrac{\cos\!\left(y_{ij}\sqrt{u^2-k_x^2}\right)}
{\sqrt{u^2-k_x^2}},&u>k_x,
\end{cases}
\qquad\text{(5)}
```

followed by ``u^2-k_x^2=\lambda^2``. No value at ``u=k_x`` is stated. ``J_0`` is the first-kind Bessel function of order zero (appendix below (10)).

Source-prescribed reductions, section 3, p. 163: making the first earth layer match air and setting ``k_x=k_0=\omega\sqrt{\varepsilon_0\mu_0}`` is identified with Kikuchi's overhead homogeneous-earth case [18]; making the two earth layers identical is identified with the authors' homogeneous-earth case [9]. That original has a [separate record](../../2010/homogeneous-earth-underground-cable-correction/Papadopoulos2010b.md). The source states that ignoring ``G`` results in propagation constant ``\gamma_1``. These are source assertions, not independently reconstructed limiting formulas.

For full cable matrices, section 2 explicitly refers to [9,15] and their procedure for adding ground terms to the appropriate conductor/insulation contributions. It does not print a new matrix-inversion assembly here. The [associated impedance record](../../../external-impedance/2011/two-layer-earth-underground-upper-layer-integral/Papadopoulos2011.md) identifies the shared kernel's output separately.

## Notation map

No renaming is performed. Original ``a`` and capital Greek ``\Delta`` are retained. SI units below identify meanings from the definitions where not explicitly printed in the source.

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``P_{e_{ij}}`` | unchanged | Earth-return potential coefficient | ``\mathrm{m/F}`` equation-implied by (7a), not separately printed |
| ``Y'_{e_{ij}}`` | unchanged | Source-defined mutual per-length earth admittance | ``\mathrm{S/m}``; scalar/subscripted inverse notation retained |
| ``F,G,G_1,G_2,G_3,G_4`` | unchanged | Shared impedance kernel and radial-displacement kernels | Source (6b), (7c)–(7g); decomposition unchanged |
| ``S_{mn},D_{mn},A_{mn},\Delta_{mn}`` and primed counterparts | unchanged | Interface coefficient families before/after the transform | Source (19)–(22); ``D`` and ``\Delta`` distinct |
| ``a'_k,a_m`` | unchanged | Original/transformed vertical spectral factors | Inverse length; transformed index mismatch unresolved |
| ``\gamma_k`` | unchanged | Bulk medium propagation constant | Inverse length; air 0, first earth 1, second earth 2 |
| ``\gamma_x,k_x,k_0`` | unchanged | Imposed longitudinal constant; associated wavenumber; free-space limit wavenumber | Inverse length; ``e^{-\gamma_x x}`` |
| ``\mu_k,\varepsilon_k,\sigma_k`` | unchanged | Material permeability, permittivity, conductivity | ``\mathrm{H/m},\mathrm{F/m},\mathrm{S/m}``; scalar by source model |
| ``i,j`` | unchanged | Cable labels | Distinct from medium indices |
| ``h_i,h_j`` and displayed ``h_1,h_2`` | unchanged | Burial depths of the two cable axes | Length, downward from surface; source subscript change retained |
| ``d,y_{ij}`` | unchanged | Upper earth-layer thickness; horizontal cable separation | Length; self uses outermost radius instead of separation |
| ``x,z`` | unchanged | Longitudinal coordinate; appendix vertical coordinate | Length; appendix comparisons use ``d-h_i`` |
| ``u,\lambda`` | unchanged | Original/transformed integration variables | Inverse length; transform ``u^2-k_x^2=\lambda^2`` |
| ``J_0`` | unchanged | First-kind Bessel function of order zero | Dimensionless |
| ``m,n,k`` | unchanged | Medium indices | 0,1,2; printed transformed ``m/k`` discrepancy retained |
| ``\omega,f,j`` | unchanged | Angular frequency, frequency, imaginary unit | ``\omega=2\pi f``; ``j^2=-1`` |

## Evidence and approximation sources

- Original contribution: section 1 explicitly identifies two-layer earth shunt corrections as the new target, contrasting earlier two-layer impedance work and the authors' earlier homogeneous-earth Z/Y model. This record does not assign priority from a later book.
- Geometry, longitudinal factor and approximation: Fig. 1, parent (2)/(4), section 2's right column and transform (5), p. 162.
- Complete final potential/admittance expression: (7a)–(7g) and shared (6b), p. 163. The four ``G`` formulas are contributions to one kernel, not separate new formulation records.
- All required interface definitions: appendix (19)–(22) and the stated primed-to-unprimed transform. The source's separate ``D`` and ``\Delta`` symbols were image-checked, not inferred from extracted text that confuses them.
- Self substitution and the actual quantity being computed: p. 163 after (7). Potential/admittance conversion (7a) is author-printed; replacing its scalar notation with a matrix operation would need separate evidence from the referenced assembly sources.
- Quasi-TEM applicability and treatment of radial displacement: section 3. Published numerical comparisons cover selected stratified/homogeneous soils and cable cases in sections 4–5; they do not establish a universal accuracy bound.
- Appendix final paragraph states that both depth orderings produce (6)–(7); no additional formula is manufactured by swapping depths or simplifying the printed kernel.

## Limitations and discrepancies

- **Suspected published indexing defect:** the transformed definition prints ``a_m`` with ``\gamma_k``; the appendix primed definition uses matching ``k``. Both are preserved without choosing a repair.
- **Not stated:** square-root branch/radiation prescription and explicit air ``\sigma_0`` value in the inspected indexed definition. No zero, sign, or branch is fabricated.
- **Unresolved scalar/matrix interpretation:** (7a) prints ``P_{e_{ij}}^{-1}`` with element subscripts. The paper references full cable-matrix assembly elsewhere. This record does not identify elementwise inverse with a general matrix inverse.
- **Unreconciled appendix/final dependency:** final (7d)–(7g) uses ``A_{10}A_{12}-\Delta_{10}\Delta_{12}e^{-2a_1d}``; appendix (15)–(16) instead prints ``A'_{10}A'_{12}e^{a'_1d}-\Delta'_{10}\Delta'_{21}e^{-a'_1d}``. The differing ``12/21`` subscripts are not normalized away or used to change a sign. Original-image transcription does not establish that these parent and final expressions agree algebraically.
- **Other parent notation conflict:** appendix (17) and (18) have different second products in the printed ``K'`` denominator, documented in the companion impedance record. Unsubscripted appendix ``h`` is not silently substituted for a chosen cable depth.
- **Scope restriction:** suggested extension to cables in different earth layers, mixed arrangements, or additional layers is future coverage, not a printed result to be inferred from the paper's broad prose.
- **Conversion availability:** no Markdown counterpart was found; the PDF is the inspected witness.
