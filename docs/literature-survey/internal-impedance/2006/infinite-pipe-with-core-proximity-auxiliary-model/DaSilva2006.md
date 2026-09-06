# da Silva–Fernández–Rivas infinite-pipe model with core proximity

## Identity and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Round solid cores at arbitrary offsets inside an infinite-wall circular pipe; no outer pipe radius, insulation radius, or bonding network is included. |
| Calculated quantities | Self ``Z(i,i)`` and mutual ``Z(j,i)`` loop impedances for proposed Method 2; core skin and pairwise proximity, infinite-wall pipe surface and pipe-mediated coupling contributions |
| Earth structure | Not applicable; no earth boundary or layers in the extracted model. |
| Model and approximation | Method 2 combines core skin, core-to-core proximity, and infinite-wall pipe coupling. The positive-order series remain infinite, and the source gives no truncation error bound. |
| Main source | D. da Silva, G. Fernández, and R. A. Rivas (2006), Method 2. |
| Citation key(s) | `:DaSilva2006` |
| Evidence status | English equation pages checked; the Spanish pages corroborate them. The summation bound and root branch remain unresolved. |

**Description.** Per-unit-length self and mutual core–pipe loop impedance for solid round conductors inside an infinitely thick conducting cylindrical pipe. The model combines core skin and pairwise core-proximity terms with pipe surface and pipe-mediated coupling terms. The source's hollow-region inductive contribution is included in the complete loop, not classified as earth return.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not stated as an imposed parameter or axial exponential. The formulas use cross-sectional geometry and conductor/pipe penetration constants; they do not prescribe ``Γ=0`` explicitly. | Equation-implied — (1)–(6), p. 2; no impressed axial propagation definition in nomenclature, p. 1. |
| Air propagation constant ``γ_air`` | Not applicable to the extracted internal core–pipe model; no outside-air wave kernel or bulk propagation constant is supplied. | Equation-implied — (1)–(6); loop return definition in nomenclature, p. 1. |
| Earth propagation constant ``γ_earth`` | Not applicable; ``m_g`` is the pipe parameter, not an earth constant. | Stated — nomenclature, p. 1; (6), p. 2. |
| Earth permittivity and displacement current | Not applicable. The absence of permittivity in (5),(6) concerns the conductor/pipe penetration model, not an earth approximation. | Equation-implied — (5),(6), p. 2. |
| Range of validity | Method 2 assumes infinite pipe thickness and includes core-to-core proximity. The author relates agreement with finite-pipe results to penetration depth becoming smaller than physical pipe thickness; this is not a universal frequency bound or a printed small-parameter expansion. Nonmagnetic pipe errors near 50–60 Hz and the separate magnetic-pipe examples are discussed. The printed outer bound ``P-1`` in the self proximity sum conflicts with counting all other cores. | Stated — abstract, pp. 1–2 section III, sections IV–V pp. 3–5; equation-implied counting issue — (1) and nomenclature definition of ``P``. |
| Earth permeability ``μ_earth`` | Not applicable. ``\mu_{rg}`` is relative pipe permeability; ``\mu_0`` is free-space permeability. | Stated — nomenclature, p. 1. |
| Arrangement | Self core–pipe loop and mutual interaction of two core–pipe loops; no external overhead/underground classification. Published examples have one, two or three solid cores, including cradle and triangular placement. | Stated — nomenclature and section IV, pp. 1–5. |
| Earth structure | Not applicable; no earth boundary or layers in the extracted model. | Equation-implied — (1)–(6). |
| Conductor and insulation geometry | Solid round core ``i`` has radius ``a_i`` and pipe-axis offset ``b_i``; ``b_{ik}`` is intercore axis distance; pipe inner radius ``c_1``. ``P`` is the number of cores. ``\phi_{j,i}`` is the angle between offset vectors about the pipe axis. Method 2 contains no outer pipe radius ``c_2``. No independent insulation-layer radius or bonding network is part of these equations. | Stated — nomenclature, p. 1; Figs. 1,2,7,8; equation-implied — (1)–(6). |
| Constitutive and field assumptions | Scalar core/pipe conductivities and relative permeabilities; the conductive penetration definitions omit a ``j\omega\varepsilon`` displacement term. Spatially uniform scalar material parameters are equation-implied by the Bessel arguments; the source does not give a complete statement of linearity, isotropy, or axial field conditions. Magnetic-vector-potential parent construction is stated. Core proximity is included by the printed pairwise terms, not an explicitly complete interaction expansion. | Stated — section III, p. 2; equation-implied — (1),(2),(4)–(6). |
| Conventions | ``Z(i,i)`` is the self impedance with pipe return; ``Z(j,i)`` is the source's mutual loop label. Positive ``j\omega`` factors; explicit time exponential, axial sign convention and square-root branch not stated. Output per length in ``\Omega/\mathrm m``. | Stated — nomenclature and section IV matrices/tables; equation-implied — (1)–(6). |

**Expression.** Method 2 self impedance, equation (1), English p. 2. The source's additive order and the two distinct summations are preserved, including the printed outer upper bound ``P-1``:

```math
\begin{aligned}
Z(i,i)={}&\frac{j\omega\mu_0}{2\pi}
\ln\left(\frac{c_1^2-b_i^2}{c_1a_i}\right)
+\frac{m_i}{2\pi a_i\sigma_i}\frac{I_0(m_ia_i)}{I_1(m_ia_i)}\\
&+\sum_{\substack{k=1\\k\ne i}}^{P-1}
\frac{j\omega\mu_0}{\pi}\sum_{n=1}^{\infty}
\left(\frac{a_k}{b_{ik}}\right)^{2n}D_n(k)\\
&+\frac{m_g}{2\pi c_1\sigma_g}\frac{K_0(m_gc_1)}{K_1(m_gc_1)}
+\frac{1}{\sigma_g}\sum_{n=1}^{\infty}
\left(\frac{b_i}{c_1}\right)^{2n}C_n.
\end{aligned}\tag{1}
```

Mutual impedance, equation (2):

```math
\begin{aligned}
Z(j,i)={}&\frac{j\omega\mu_0}{4\pi}
\ln\left(
\frac{c_1^4+(b_ib_j)^2-2(b_ib_j)c_1^2\cos(\phi_{j,i})}
{c_1^2(b_i^2+b_j^2-2b_ib_j\cos(\phi_{j,i}))}
\right)\\
&+\frac{j\omega\mu_0}{\pi}\sum_{n=1}^{\infty}
\left(\frac{a_j}{b_{ij}}\right)^{2n}D_n(j)\\
&+\frac{m_g}{2\pi c_1\sigma_g}\frac{K_0(m_gc_1)}{K_1(m_gc_1)}\\
&+\frac{1}{\sigma_g}\sum_{n=1}^{\infty}
\left(\frac{b_jb_i}{c_1^2}\right)^n C_n\cos(n\phi_{j,i}).
\end{aligned}\tag{2}
```

Complete coefficient and penetration definitions, equations (3)–(6):

```math
C_n=\frac{m_g^2K_n(m_gc_1)}
{\pi[n\mu_{rg}K_n(m_gc_1)-m_gc_1K'_n(m_gc_1)]},\tag{3}
```

```math
D_n(k)=\frac{(a_k)^{-1}I_n(m_ka_k)}
{\dfrac{n}{a_k}I_n(m_ka_k)+\dfrac{m_k}{\mu_{rk}}I'_n(m_ka_k)},\tag{4}
```

```math
m_i=\sqrt{j\omega\mu_0\mu_{ri}\sigma_i},\tag{5}
```

```math
m_g=\sqrt{j\omega\mu_0\mu_{rg}\sigma_g}.\tag{6}
```

``m_k`` and ``m_j`` use the same core-labelled definition (5) for the corresponding core. The source's ``D_n(j)`` in (2) means the ``k=j`` instance of (4); it is not symmetrized between the two cores. ``I_n,K_n`` are the modified Bessel functions of the first and second kind, order ``n``. Nomenclature explicitly identifies ``I'_n,K'_n`` as their first derivatives; their arguments are retained as printed, without replacement by recurrence identities. ``m_i,m_g`` are called reciprocals of complex penetration depth and have inverse-length dimensions. ``C_n`` has inverse-area dimensions and ``D_n`` is dimensionless, as implied by their equations.

``\mu_0`` is free-space permeability; ``\mu_{ri},\mu_{rg}`` are relative material permeabilities, ``\sigma_i,\sigma_g`` conductivities, and ``\omega`` angular frequency in radians per second. The radii/axis distances are lengths and ``\phi`` is an angle. No square-root branch, finite series cutoff, numerical convergence rule, dielectric permittivity, or external reference transformation is supplied with these expressions.

**Approximation.** A proposed auxiliary combination: infinite-thickness pipe representation with core skin and core-to-core proximity included. Section III states that its magnetic-vector-potential formulation draws on Tegopoulos–Kriezis [7], Brown–Rocamora [2], and Kane's thesis [5]. The paper does not print a stepwise limit ``c_2\to\infty`` of a chosen parent equation, an expansion parameter/order, or discarded interaction series. The physical infinite-wall substitution is explicit; no reviewer-derived asymptotic calculation is added. The positive-order series remain infinite.

**Limitations.** The self proximity sum has an indexing inconsistency under the stated core count. Core-labelled mutual proximity remains as printed and is not forced symmetric for unequal cores. Infinite-wall behavior is not generally a low-frequency finite-pipe approximation. The source does not specify the square-root branch.

**Reference.** [DaSilva2006](@cite). Da Silva, Fernández, and Rivas (2006), DOI `10.1109/TDCLA.2006.311519`; English p. 1 nomenclature and attribution, p. 2 section III.A equations (1)–(6), and pp. 3–6 comparison domain. The Spanish version repeats the equation set on its printed p. 2/PDF page 8.

**Transcription source.** Original conference-publication English page images control this transcription. All equations (1)–(6), coefficients, powers and limits were visually checked, with a high-resolution crop of the self sum and mutual/kernel expressions. Spanish PDF page 8 independently repeats the same equation set and upper bound.

## Source transcription

The formula section displays the complete source-notation set in original equation order (1)–(6); no mathematical renaming is used. English section III.A is on the left of p. 2, followed by section III.B's distinct finite-pipe model. The numbered ``(1)`` printed just after ``C_n`` is the equation label, not an argument ``C_n(1)``. The Spanish set is likewise (1)–(6), at its printed p. 2/PDF page 8.

The required definitions from English nomenclature p. 1 are explicitly distinguished: ``b_i`` measures conductor-axis to pipe-axis distance, while ``b_{ik}`` is core-to-core distance; ``P`` is the total number of conductors inside the pipe. Both distinctions matter for the two different powers and the core-count issue in (1). ``\phi_{j,i}`` is an angle around the pipe center, not a spatial crossing angle between nonparallel cables.

The source does not print a separate finite outer-surface or through-wall transfer impedance for Method 2. Nor does it derive admittance from the complete loop impedance. Its field-region logarithms overlap the geometric constituents of the [Kane record](../../1995/finite-shield-core-loop-skin-and-proximity-series/Kane1995.md), but the infinite-pipe combination and the printed sum are retained as a distinct published model. The [companion Method 3 record](../finite-pipe-without-core-proximity-auxiliary-model/DaSilva2006.md) uses different pipe coefficients and omits core-to-core proximity; the two records are not blended.

## Notation map

Notation is unchanged. The letter ``j`` is a core label in ``Z(j,i),a_j,D_n(j)`` and the imaginary unit in ``j\omega``; the distinct uses are not normalized into new symbols.

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z(i,i),Z(j,i)`` | unchanged | Self/mutual core–pipe loop impedances | ``\Omega/\mathrm m``; pipe return/reference |
| ``a_i,a_j,a_k`` | unchanged | Solid-core radii | Metres |
| ``b_i,b_j`` | unchanged | Core offsets from pipe axis | Metres |
| ``b_{ik},b_{ij}`` | unchanged | Core-axis separation | Metres; not an offset from pipe axis |
| ``c_1`` | unchanged | Inner radius of infinite-thickness pipe | Metres; no ``c_2`` in this model |
| ``\phi_{j,i}`` | unchanged | Angle between core-offset vectors about pipe axis | Radians; cosine and integer multiple as printed |
| ``P,i,j,k,n`` | unchanged | Total core count, core labels and positive Bessel/sum order | Self sum upper bound is literally ``P-1`` with ``k\ne i`` |
| ``\mu_0,\mu_{ri},\mu_{rk},\mu_{rj},\mu_{rg}`` | unchanged | Free-space and relative core/pipe permeabilities | ``\mathrm H/\mathrm m``; relative quantities dimensionless |
| ``\sigma_i,\sigma_j,\sigma_k,\sigma_g`` | unchanged | Core and pipe conductivities | ``\mathrm S/\mathrm m`` |
| ``m_i,m_j,m_k,m_g`` | unchanged | Reciprocals of complex material penetration depths | ``\mathrm m^{-1}``; no root branch stated |
| ``I_n,K_n,I'_n,K'_n`` | unchanged | Modified first-/second-kind Bessel functions and first derivatives | Nomenclature definitions, p. 1 |
| ``C_n`` | unchanged | Infinite-pipe coupling coefficient (3) | ``\mathrm m^{-2}``, equation-implied; not capacitance |
| ``D_n(k),D_n(j)`` | unchanged | Core-proximity coefficient (4) and its core-``j`` instance | Dimensionless |
| ``\omega,j`` in factors | unchanged | Angular frequency and imaginary unit | ``\mathrm{rad/s}``, dimensionless; time exponential not stated |

## Evidence and approximation sources

- Title/byline, conference header and all 12 pages of the bilingual file were image-inspected. The English and Spanish versions each have independently restarted printed pagination 1–6. Their displayed mathematical sets (1)–(14) were compared; no separate formula is created solely for translation. Case-study and appendix layout differs between versions, so all quoted numerical locators here identify the English pages.
- The paper explicitly proposes auxiliary Methods 2 and 3 instead of limiting the comparison to unchanged earlier formulae. The 2006 publication supplies the equations used in this record.
- The infinite-pipe approximation is an assumption of Method 2, while core-to-core terms are retained. Agreement with finite-pipe Method 4 as penetration depth becomes smaller than wall thickness is discussed in sections IV.A–B, pp. 3–4. No numerical bound on a dimensionless wall-depth ratio is supplied.
- Author-reported discrepancy example: Table VIII, English p. 6, gives Method 2 mutual-resistance relative error against Kane Method 4 as 95.01% at 1 Hz and 0.26% at 100 kHz. These are not analytic error bounds or a comparison to an exact benchmark. The two-core test has ``a_1=a_2=0.0135`` m, ``b_1=b_2=0.0815`` m, ``\phi_{1,2}=1.096`` rad, ``c_1=0.127`` m, physical ``c_2=0.133`` m, nonmagnetic cores/pipe, core conductivity ``5.80\times10^7`` and pipe conductivity ``3.57\times10^7`` S/m specified at 20°C, with core/pipe temperatures 80/60°C and temperature coefficients 0.0043/0.0042 K⁻¹ (English p. 3). The paper lists these temperature data but prints no conductivity-temperature correction equation; none is supplied here.
- The single-core example states that methods using the same pipe model coincide because there is no second core for intercore proximity. The three-core magnetic-pipe cases at 60 Hz show geometry-dependent self entries; these example matrices are not a new general reference-conductor transformation.

## Limitations and discrepancies

- **Demonstrable source-level counting issue:** nomenclature says ``P`` is the total core count, but (1) uses ``k=1,\ k\ne i`` through ``P-1``. For ``P=2,i=1`` this range has no included index, although there is another core. Both English and Spanish images print this bound. It is not corrected to the ``P`` bound printed in Kane's original (9),(11), and no hidden relabelling convention is invented.
- **Source mathematics versus results:** the above counting issue prevents treating the printed self expression as an unambiguous specification of all reported multi-core results. The tables are not used to infer corrected mathematics or missing code.
- **Mutual asymmetry retained:** the core-proximity part of (2) uses ``a_j,m_j,\mu_{rj},\sigma_j`` via ``D_n(j)``. No averaging with the reverse pair is supplied. Høidalen's later discussion of the parent asymmetry is a distinct witness, not permission to symmetrize this paper.
- **Unstated branch/constitutive extension:** the root choice and complete axial dependence are not stated. Conductor displacement current, anisotropic materials or dielectric losses are not restored to (5),(6).
- **No finite sum inferred:** the displayed positive-order sums are infinite; neither the case-study plots nor a later author's chosen cutoff establishes this paper's truncation order.
- **Conversion status:** no original-paper Markdown was located. Text extraction damages layout and could join the equation number to ``C_n``; the PDF images control the transcription.
