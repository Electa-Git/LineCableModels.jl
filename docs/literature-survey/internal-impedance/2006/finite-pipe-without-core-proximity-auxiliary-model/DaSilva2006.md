# da Silva–Fernández–Rivas finite-pipe model without core-to-core proximity

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Round solid cores at arbitrary offsets inside a finite-wall circular pipe; insulation parameters and bonding transformations are excluded. |
| Calculated quantities | Self ``Z(i,i)`` and mutual ``Z(j,i)`` impedances for proposed Method 3; finite-wall pipe skin and pipe-mediated coupling, and self core skin; no core-to-core eddy-current increment |
| Earth structure | Not applicable. |
| Model and approximation | Method 3 retains full finite-wall Bessel coefficients and omits core-to-core proximity. No expansion parameter or error bound is specified. |
| Main source | D. da Silva, G. Fernández, and R. A. Rivas (2006), Method 3. |
| Citation key(s) | `:DaSilva2006` |
| Evidence status | English equation pages checked; the Spanish pages corroborate the equations. |

**Description.** Per-unit-length self and mutual impedance of solid-core return loops through a finite-thickness cylindrical pipe. This auxiliary model retains the pipe's nonuniform surface response and core skin effect but omits the additional core-to-core proximity contribution. “Without core-to-core proximity” does not mean that the pipe-mediated coupling series is removed.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not stated as an imposed constant or axial exponential. The equations are cross-sectional core–pipe parameter expressions; no explicit ``Γ=0`` prescription is supplied. | Equation-implied — (5)–(14), p. 2; no longitudinal parameter in nomenclature. |
| Air propagation constant ``γ_air`` | Not applicable to the extracted core–pipe loop; no exterior-air propagation kernel is part of Method 3. | Equation-implied — (7)–(14); loop definitions on p. 1. |
| Earth propagation constant ``γ_earth`` | Not applicable. ``m_g`` is explicitly the pipe reciprocal complex penetration depth. | Stated — nomenclature p. 1 and definition (6). |
| Earth permittivity and displacement current | Not applicable; no earth medium in these expressions. The conductive material definitions omit a permittivity contribution. | Equation-implied — (5),(6). |
| Range of validity | Finite pipe thickness retained; proximity among the solid cores is deliberately omitted. The source reports high-frequency error from this omission, with mutual resistance most affected. Examples are one-/two-core nonmagnetic pipe and three-core magnetic-pipe arrangements, not universal frequency/geometry limits. No controlled proximity expansion order or error bound is provided. | Stated — abstract, section III.B and conclusion pp. 1,2,5; section IV test cases. |
| Earth permeability ``μ_earth`` | Not applicable; ``\mu_{rg}`` denotes the pipe material, while ``\mu_0`` is free-space permeability. | Stated — nomenclature p. 1. |
| Arrangement | Self and mutual core–pipe return loops. Not a ground-return arrangement. Single-, two- and three-core cross sections are studied. | Stated — nomenclature and section IV, pp. 1–5. |
| Earth structure | Not applicable. | Equation-implied — (7)–(14). |
| Conductor and insulation geometry | Round solid core radii ``a_i``, offsets ``b_i`` from pipe axis; pipe inner and outer radii ``c_1,c_2``. Mutual angle ``\phi_{j,i}`` is measured at pipe center. No separate insulation-layer parameters or arbitrary bonding transformation are supplied. | Stated — nomenclature and Figs. 1,2,7,8; equation-implied — (7)–(14). |
| Constitutive and field assumptions | Scalar material conductivities and relative permeabilities; Bessel penetration parameters (5),(6) retain conductive diffusion without a displacement-current term. Spatially uniform scalar parameters are equation-implied; complete linearity, isotropy, and axial field assumptions are not expressly enumerated. Magnetic-vector-potential derivation is stated. Core-to-core eddy-current terms omitted, pipe-mediated ``H_n`` series retained. | Stated — section III p. 2; equation-implied — (5)–(14), compared with Method 2 (1),(2). |
| Conventions | Pipe return/reference; self ``Z(i,i)`` and mutual ``Z(j,i)`` are loop quantities, not isolated metal surfaces. Positive ``j\omega`` convention in geometric terms; explicit time/axial exponent and root branch not stated. ``\Omega/\mathrm m`` output. | Stated — nomenclature and section IV matrices/tables; equation-implied — (7),(8). |

**Expression.** Method 3 self and mutual expressions, (7),(8), English p. 2. The source order of the geometric, core, zeroth-order pipe and positive-order pipe terms is preserved.

```math
\begin{aligned}
Z(i,i)={}&\frac{j\omega\mu_0}{2\pi}
\ln\left(\frac{c_1^2-b_i^2}{c_1a_i}\right)
+\frac{m_i}{2\pi a_i\sigma_i}\frac{I_0(m_ia_i)}{I_1(m_ia_i)}\\
&+\frac{m_g}{2\pi c_1\sigma_g}
[A_0I_0(m_gc_1)+B_0K_0(m_gc_1)]\\
&+\frac{1}{\sigma_g}\sum_{n=1}^{\infty}
\left(\frac{b_i}{c_1}\right)^{2n}H_n.
\end{aligned}\qquad\text{(7)}
```

```math
\begin{aligned}
Z(j,i)={}&\frac{j\omega\mu_0}{4\pi}
\ln\left(
\frac{c_1^4+(b_ib_j)^2-2(b_ib_j)c_1^2\cos(\phi_{j,i})}
{c_1^2(b_i^2+b_j^2-2b_ib_j\cos(\phi_{j,i}))}
\right)\\
&+\frac{m_g}{2\pi c_1\sigma_g}
[A_0I_0(m_gc_1)+B_0K_0(m_gc_1)]\\
&+\frac{1}{\sigma_g}\sum_{n=1}^{\infty}
\left(\frac{b_ib_j}{c_1^2}\right)^n H_n\cos(n\phi_{j,i}).
\end{aligned}\qquad\text{(8)}
```

The distinct zeroth-order coefficients are (9),(10); ``B_0`` has the printed positive numerator:

```math
A_0=\frac{K_1(m_gc_2)}
{I_1(m_gc_1)K_1(m_gc_2)-I_1(m_gc_2)K_1(m_gc_1)},\qquad\text{(9)}
```

```math
B_0=\frac{I_1(m_gc_2)}
{I_1(m_gc_1)K_1(m_gc_2)-I_1(m_gc_2)K_1(m_gc_1)}.\qquad\text{(10)}
```

Positive-order pipe coefficients and the complete denominator, (11)–(14):

```math
H_n=A_nI_n(m_gc_1)+B_nK_n(m_gc_1),\qquad\text{(11)}
```

```math
A_n=\frac{m_g}{\pi c_1\Delta_n}
\left[\frac{n\mu_{rg}}{m_gc_2}K_n(m_gc_2)+K'_n(m_gc_2)\right],\qquad\text{(12)}
```

```math
B_n=-\frac{m_g}{\pi c_1\Delta_n}
\left[\frac{n\mu_{rg}}{m_gc_2}I_n(m_gc_2)+I'_n(m_gc_2)\right],\qquad\text{(13)}
```

```math
\begin{aligned}
\Delta_n={}&
\left[\frac{n\mu_{rg}}{m_gc_1}I_n(m_gc_1)-I'_n(m_gc_1)\right]
\cdot\left[\frac{n\mu_{rg}}{m_gc_2}K_n(m_gc_2)+K'_n(m_gc_2)\right]\\
&-\left[\frac{n\mu_{rg}}{m_gc_2}I_n(m_gc_2)+I'_n(m_gc_2)\right]
\cdot\left[\frac{n\mu_{rg}}{m_gc_1}K_n(m_gc_1)-K'_n(m_gc_1)\right].
\end{aligned}\qquad\text{(14)}
```

Their required penetration definitions appear immediately before Method 3 in the source:

```math
m_i=\sqrt{j\omega\mu_0\mu_{ri}\sigma_i},\qquad\text{(5)}
```

```math
m_g=\sqrt{j\omega\mu_0\mu_{rg}\sigma_g}.\qquad\text{(6)}
```

Nomenclature p. 1 explicitly defines ``I_n,K_n`` as modified first-/second-kind Bessel functions, and their primes as first derivatives. ``m_i,m_g`` are reciprocal complex penetration depths. ``\mu_0`` is free-space permeability, ``\mu_{ri},\mu_{rg}`` relative core/pipe permeabilities, and ``\sigma_i,\sigma_g`` their conductivities. All radii/axis offsets are lengths, ``\phi_{j,i}`` is an angle and ``\omega`` angular frequency. ``A_0,B_0,\Delta_n`` are dimensionless, while positive-order ``A_n,B_n,H_n`` have inverse-area dimensions as implied by (11)–(14). The separate zeroth-order coefficients must not be replaced by an unprinted ``n=0`` substitution of the positive-order formulas.

No finite-thickness skin-depth expansion, source-defined root branch, numerical sum cutoff, external earth impedance, or dielectric-loss model accompanies this set. The complete output remains the loop impedance printed in (7),(8), not a newly isolated pipe transfer impedance.

**Approximation.** A proposed auxiliary finite-pipe model with core-to-core proximity omitted. The paper states the physical omission and gives the resulting self/mutual equations; it does not specify a small geometric parameter, interaction order or discarded higher-order remainder. The pipe boundary dependence on both ``c_1`` and ``c_2`` is retained through full Bessel coefficients. This is not an infinite-pipe formula or a claim of exact complete multi-core electromagnetic response.

**Limitations.** Missing core-to-core proximity affects multi-core high-frequency comparisons. Pipe-mediated coupling remains present. The source does not specify the square-root branch. Separate outer-surface, transfer, and terminal-bonding formulae are not supplied by this set.

**Reference.** [DaSilva2006](@cite). Da Silva, Fernández, and Rivas (2006), DOI `10.1109/TDCLA.2006.311519`; English nomenclature p. 1, section III.B p. 2 equations (7)–(14) with (5)–(6), and case studies, conclusion, and appendix on pp. 3–6. The Spanish version repeats the equation set on its printed p. 2/PDF page 8.

**Transcription source.** Original conference-publication page images in the English part of the source PDF. Equations (5)–(14), the boundary-product signs, distinct zeroth/positive orders and geometry definitions were visually inspected and compared with Spanish PDF page 8. This does not verify earlier Kane/Tegopoulos equations by attribution alone.

## Source transcription

The complete set above uses original notation. In source order, the shared material definitions (5),(6) precede (7),(8), then (9),(10), (11), (12),(13), and (14). English (7) occupies the bottom left of p. 2, with (8)–(14) at the top right; reading columns in the wrong order would separate the self/mutual set from its definitions. The Spanish version groups Method 3 in the right column of its p. 2.

The absence of the ``D_n`` core-proximity terms in (7),(8) is deliberate, as section III states. The positive-order ``H_n`` terms are explicitly present in both formulas. Therefore replacing “no proximity among conductors” by “uniform pipe current” would contradict the printed expressions.

The source's method comparison is a two-assumption comparison, not a new terminal-network operation: Method 1 is described as infinite pipe without core proximity, Method 2 infinite pipe with it, Method 3 finite pipe without it, and Method 4 finite pipe with it. The [separate Method 2 record](../infinite-pipe-with-core-proximity-auxiliary-model/DaSilva2006.md) and [Kane original record](../../1995/finite-shield-core-loop-skin-and-proximity-series/Kane1995.md) retain those distinct witnesses. No coefficients are copied between them to make a reconstructed fifth model.

## Notation map

No mathematical renaming is performed. In particular, ``A_0,B_0`` are distinct source-provided zeroth-order definitions, and ``\Delta_n`` retains the order subscript that differs from Kane's original notation.

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z(i,i),Z(j,i)`` | unchanged | Self/mutual core–pipe loop impedances | ``\Omega/\mathrm m``; pipe return/reference |
| ``a_i`` | unchanged | Solid-core radius in self skin/geometric terms | Metres |
| ``b_i,b_j`` | unchanged | Core-axis distances to pipe axis | Metres |
| ``c_1,c_2`` | unchanged | Inner/outer pipe radii | Metres; finite wall |
| ``\phi_{j,i}`` | unchanged | Core-offset angular separation about pipe center | Radians |
| ``\mu_0,\mu_{ri},\mu_{rg}`` | unchanged | Free-space and relative core/pipe permeabilities | ``\mathrm H/\mathrm m`` and dimensionless relatives |
| ``\sigma_i,\sigma_g`` | unchanged | Core/pipe conductivities | ``\mathrm S/\mathrm m`` |
| ``m_i,m_g`` | unchanged | Reciprocals of complex penetration depths | ``\mathrm m^{-1}``; source branch not stated |
| ``I_n,K_n,I'_n,K'_n`` | unchanged | Modified Bessel functions of first/second kind and first derivatives | Source nomenclature, p. 1 |
| ``A_0,B_0`` | unchanged | Separate zeroth-order pipe coefficients | Dimensionless |
| ``A_n,B_n,H_n`` | unchanged | Positive-order pipe coefficients and combination | ``\mathrm m^{-2}``, equation-implied |
| ``\Delta_n`` | unchanged | Order-dependent finite-pipe boundary denominator | Dimensionless |
| ``n`` | unchanged | Positive-order summation index | One through infinity, no stated cutoff |
| ``\omega,j`` | unchanged | Angular frequency and imaginary unit | ``\mathrm{rad/s}``, dimensionless; explicit time exponential not stated |

## Evidence and approximation sources

- The abstract and section III identify this as an auxiliary finite-pipe model that omits core proximity. These assumptions differ from Method 2.
- The stated parent development uses magnetic vector potential and cites Tegopoulos–Kriezis, Brown–Rocamora, and Kane's thesis. Earlier shell papers remain separate candidates.
- Section IV.A states that, for its single-core case, Methods 3 and 4 coincide because there is no intercore proximity contribution. This is an author-stated geometry-specific comparison, not a proof of exactness for arbitrary multi-core arrangements.
- In the two-core nonmagnetic-pipe example, Table VIII (English p. 6) reports Method 3 mutual-resistance discrepancies from Kane Method 4 of 26.43% at 1 kHz and 27.25% at 100 kHz. The tested geometry is ``a_1=a_2=0.0135`` m, ``b_1=b_2=0.0815`` m, ``\phi_{1,2}=1.096`` rad, ``c_1=0.127`` m and ``c_2=0.133`` m, with all relative permeabilities one. Core/pipe conductivities are specified at 20°C as ``5.80\times10^7``/``3.57\times10^7`` S/m, temperatures 80/60°C and coefficients 0.0043/0.0042 K⁻¹ (English p. 3). A temperature-correction equation is not printed. These numbers are reported comparisons to another model, not guaranteed error bounds or independent measurement validation.
- The conclusion p. 5 states that omission of core proximity causes high-frequency error, with mutual resistance most affected. The magnetic-pipe three-core examples at 60 Hz do not prescribe a universal permeability-dependent cutoff.

## Limitations and discrepancies

- **Source-defined omission:** no core-to-core ``D_n`` term belongs to Method 3. Adding the Kane/Method 2 increment would change the recorded physical model.
- **Pipe coupling retained:** the ``H_n`` sums must remain, with both radii in (9)–(14). Neglect of core proximity is not evidence that these pipe coefficients vanish.
- **Boundary and zeroth-order fidelity:** the positive numerator of ``B_0``, the negative prefactor of positive-order ``B_n``, and the four signed factors of ``\Delta_n`` were checked in both language versions. They are not replaced by scaled Bessel functions, recursion identities or a later normalized form.
- **Scope of matching language witnesses:** the English and Spanish equation sets agree on the inspected (5)–(14). This is not a claim of word-for-word equality of their surrounding prose or all case-study pagination.
- **Root and field definitions:** the source does not supply the root branch, explicit axial dependence, or a complete list of linearity and isotropy restrictions. Parent-source assumptions are not assigned to the 2006 formulation.
- **No original-source repair:** this paper's nomenclature can corroborate what later authors mean by primed Bessel functions, but it cannot retroactively change a missing definition or a printed defect in Kane's original.
