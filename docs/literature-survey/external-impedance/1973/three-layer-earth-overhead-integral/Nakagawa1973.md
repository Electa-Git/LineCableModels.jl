# Nakagawa–Ametani–Iwamoto three-layer overhead earth-return impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinite parallel line currents above stratified earth; conductor skin, insulation, and a general finite-radius self term are excluded. |
| Calculated quantities | Mutual external and earth-return impedance, own-image self correction, and two-layer and homogeneous reductions. |
| Earth structure | Air above two finite earth layers and a lower earth half-space; ``d_2`` is cumulative depth. No arbitrary-layer recursion is supplied. |
| Model and approximation | Integral representation with the source's longitudinal prescription ``\gamma_0=jk`` and Bessel-transform reduction. The two- and one-layer formulae are geometric limits; no approximation error bound is supplied. |
| Main source | M. Nakagawa, A. Ametani, and K. Iwamoto (1973). |
| Citation key(s) | Primary: `:Nakagawa1973`; later witness: `:Ametani1975` |
| Evidence status | All 1973 and 1975 pages checked. Appendix grouping, derivation defects, the root branch, and the finite-radius self prescription remain unresolved. |

**Description.** Per-length mutual external impedance and earth-return correction for parallel overhead conductors above three horizontal earth layers, with explicit two-layer and homogeneous limits and an own-image self earth-correction formula.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | The current propagation constant is assigned pure imaginary ``jk`` and the Hertz vector is proportional to ``\exp(-\gamma x)``; the final reduction explicitly takes ``\gamma_0=jk``. This is a prescribed longitudinal approximation, not ``Γ=0`` or an unknown solved by a dispersion relation. The source reuses bulk/longitudinal gamma notation; the role distinction is recorded without changing its equations. | Stated — §2.1 and §2.3, p. 1522, prose before (5) and (12); appendix §8.2. |
| Air propagation constant ``γ_air`` | ``\gamma_0=jk`` is retained. The indexed appendix relation is ``\gamma_i^2=-k_i^2=j\omega\mu_i(\sigma_i+j\omega\varepsilon_i)``. Air is index 0; the main expression writes unindexed ``k``. No separate explicit value of air conductivity is printed with that prescription. The pure-imaginary relation is not a zero-air-propagation approximation. | Stated — p. 1521 symbols, Fig. 2, p. 1522 prescription, p. 1527 definitions. Equation-implied — lossless bulk propagation under the prescribed pure-imaginary constant and real scalar constitutive inputs; not elevated to a printed ``\sigma_0=0`` statement. |
| Earth propagation constant ``γ_earth`` | ``\gamma_i^2=-k_i^2=j\omega\mu_i(\sigma_i+j\omega\varepsilon_i)``, earth indices 1,2,3. Final ``a_i=\sqrt{v^2+k^2-k_i^2}``, ``b_i=a_i/\mu_i``. The square-root branch is not explicitly stated. | Stated — definitions after (14), p. 1523, and after (19), p. 1527. |
| Earth permittivity and displacement current | Layer-specific scalar permittivity and conductivity are separately retained in the original bulk relation. Some numerical graphs expressly neglect displacement current; later captions use ``\varepsilon_i=\varepsilon_0`` for corresponding comparisons. These tested formulas do not remove permittivity from the main kernel. No dielectric loss-tangent or complex-permittivity/conductivity bookkeeping rule is specified. | Stated — introduction, §3.1–3.2, Figs. 6–9,11; equation-implied — p. 1527 bulk relation. |
| Range of validity | Infinite-line Hertz-vector construction with the prescribed plane-wave longitudinal approximation. No analytical truncation order or uniform error bound is stated. The authors associate visible displacement-current differences at ``w>1`` with frequencies above 1 MHz in their high-resistivity/low-height examples; this is not a universal cutoff. They explicitly question the plane-wave assumption for some low-``w`` results. | Stated — §2.1, p. 1522; §3.1, p. 1523, final paragraphs; §3.2, p. 1524; abstract. |
| Earth permeability ``μ_earth`` | Independent layer values ``\mu_1,\mu_2,\mu_3`` retained. Unity relative permeability is an example choice, not a kernel restriction. Ferromagnetic middle-layer examples use ``\mu_2=100\mu_0`` and ``500\mu_0``; no dispersive or hysteretic permeability law is supplied. | Stated — introduction, (14), §3.3 and Fig. 9, pp. 1521,1523,1525–1526. |
| Arrangement | Both source and target conductors in air, parallel to ``x``; lateral separation ``y`` and heights ``h_1,h_2``. Mutual ``Z_{12}``; self earth correction uses the reference conductor's own image (``y=0,h_1=h_2=h`` is geometry-implied). No buried or mixed-placement formula is printed. | Stated — Figs. 1–2, §2.3, and §3/§3.1 own-image/self explanation, pp. 1521–1523. |
| Earth structure | Air 0; finite upper earth 1, finite middle earth 2, lowest earth half-space 3. Interfaces at ``z=0,-d_1,-d_2``. ``d_2`` is cumulative upper-plus-middle depth, not middle-layer thickness. Two layers by ``d_2\to\infty``; homogeneous limit by subsequently ``d_1\to\infty``. No arbitrary-layer recursion is supplied. | Stated — Fig. 2 and §2.2, p. 1522; p. 1523 following (14)–(15); appendix interface locations p. 1527. |
| Conductor and insulation geometry | Infinite parallel line-current construction from horizontal dipoles. Height is above the top interface; the output excludes a separately evaluated conductor skin/proximity or insulation contribution. Numerical wire radii do not establish a general finite-radius self substitution in the logarithmic mutual term. | Stated/equation-implied — §2.1–2.3; Figs. 11–12 example radii, p. 1526. Unresolved — general full-external self geometric prescription. |
| Constitutive and field assumptions | Uniform scalar properties within each horizontal medium and linear field superposition are equation-implied. Source invokes continuity of tangential fields and prints ``H_y=0`` in the appendix. The printed derivative/component mismatches in those conditions are preserved, not silently repaired. No nonlinear earth or conductor model is specified. | Stated — §2.2 and appendix §8.1, p. 1527; equation-implied — scalar bulk constants and dipole integration. |
| Conventions | ``x`` longitudinal, ``z`` upward from ground surface, ``y`` lateral. The source's positive ``j\omega`` harmonic factors and negative longitudinal exponential are retained; a full time exponential is not explicitly printed. ``Z_e`` is the earth correction and ``Z_{12}`` includes the geometrical external term. SI constitutive quantities imply impedance in Ω/m; graph resistivities are Ωm and lengths m. | Stated — (3)–(4), Figs. 1–2 and captions; equation-implied — per-length longitudinal field/current relation and dimensions. Not stated — explicit time exponential and root branch. |

**Expression.** Original mutual external impedance and its separately printed earth correction, (12)–(13), p. 1522:

```math
Z_{12}=j\omega\frac{\mu_0}{2\pi}\ln\left(\frac{D_2}{D_1}\right)+Z_e.
\qquad\text{(12)}
```

```math
Z_e=\omega\frac{\mu_0}{2\pi}(P+jQ)
=j\omega\frac{\mu_0}{2\pi}\,2\int_0^\infty
B_2\exp\{-(h_1+h_2)v\}\cos(yv)\,dv.
\qquad\text{(13)}
```

```math
D_1=\sqrt{y^2+(h_1-h_2)^2},\qquad
D_2=\sqrt{y^2+(h_1+h_2)^2}.
```

The three-layer coefficient and complete final spectral definitions, (14) and following lines, p. 1523:

```math
B_2=\frac{c_1+c_2}{(v+\mu_0b_1)c_1+(v-\mu_0b_1)c_2},
\qquad\text{(14)}
```

```math
c_1=(b_1+b_2)(b_2+b_3)
 +(b_1-b_2)(b_2-b_3)\exp\{2a_2(d_1-d_2)\},
```

```math
c_2=\bigl[(b_1-b_2)(b_2+b_3)
 +(b_1+b_2)(b_2-b_3)\exp\{2a_2(d_1-d_2)\}\bigr]
 \exp(-2a_1d_1),
```

```math
a_i=\sqrt{v^2+k^2-k_i^2},\qquad b_i=a_i/\mu_i,
\qquad i=1,2,3.
```

Required original propagation definitions, appendix p. 1527 and the p. 1522 prescription:

```math
\gamma_i^2=-k_i^2=j\omega\mu_i(\sigma_i+j\omega\varepsilon_i),
\qquad \gamma_0=jk.
```

No branch is added. In the original notation ``k`` is the unindexed air wavenumber of the imposed relation, ``k_i`` belongs to medium ``i``, and ``v`` is the final inverse-length spectral variable. ``\sigma_i`` is conductivity (S/m), ``\varepsilon_i`` permittivity (F/m), ``\mu_i`` permeability (H/m); ``\omega`` is angular frequency. ``P,Q`` are the dimensionless real/reactive correction quantities in the **printed** ``\omega\mu_0(P+jQ)/(2\pi)`` decomposition; the first prefactor must not acquire an extra ``j``. ``B_2`` has length units, so its integral is dimensionless. ``h_1,h_2,D_1,D_2,y,d_1,d_2`` are lengths. Depth ``d_2`` includes both finite earth layers.

For the two-layer formula the source takes ``d_2\to\infty`` and replaces ``B_2`` in (13) by (15):

```math
B_3=\frac{b_1+b_2+(b_1-b_2)\exp(-2a_1d_1)}
{(v+\mu_0b_1)(b_1+b_2)
 +(v-\mu_0b_1)(b_1-b_2)\exp(-2a_1d_1)}.
\qquad\text{(15)}
```

The subsequent homogeneous ``d_1\to\infty`` replacement is printed, unnumbered below (15), as

```math
(v+\mu_0b_1)^{-1}.
```

These are related formulas, not independent priority assignments to Nakagawa et al. for earlier two-/one-layer results. For the self **earth correction**, §3's own-image definition and §3.1's self designation give ``y=0,h_1=h_2=h,D_2=2h`` by geometry. No ``D_1``-to-radius prescription is inserted into the full mutual (12); that separate finite-radius self quantity remains unresolved.

**Approximation.** Integral representation within an explicitly restricted longitudinal/plane-wave model, not a solved full-wave dispersion relation. Section 2.1 assigns pure-imaginary longitudinal ``jk`` following Wise, and §2.3/appendix §8.2 apply ``\gamma_0=jk`` and a Bessel-transform reduction to obtain (12)–(14). No expansion parameter, retained order or error bound for this prescription is given here. The two-/one-layer formulas are source-stated geometric limits, not fitted equivalent-earth approximations. The appendix's infinite series is a transform-identity step, not a new truncated parameter formulation.

**Limitations.** Three earth layers maximum in the inspected expression set. The root branch and complete finite-radius self geometric term are not established. Ambiguous appendix ``B_1`` slash grouping and published field/transform inconsistencies remain visible below; the clear final ``B_2`` fraction is not a license to repair its parent. Later 1975 notation and a defective constitutive token are retained separately. No arbitrary-layer recursion, admittance, buried/mixed kernel or conductor-internal term is inferred.

**Reference.** [Nakagawa1973](@cite), §§2.1–2.3, equations (12)–(15), pp. 1521–1523, and appendix §§8.1–8.2, pp. 1527–1528; later witness [Ametani1975](@cite), equations (1)–(6), pp. 500–501, DOI `10.1541/ieejpes1972.95.500`.

**Transcription source.** Original 1973 journal page images are the mathematical authority; all eight pages were inspected. All seven pages of the 1975 publication were inspected as evidence of a later restatement. Neither an LCM implementation nor a Markdown conversion was used to reconstruct equations.

## Source transcription

Notation is unchanged throughout; the formula section gives original (12)–(15) and the required final definitions. The following additional equations follow their source order, with explicit locators. They document the declared parent, field conventions and reduction rather than separate new formulas.

### Original 1973 field and integral parents

Section 2.1, p. 1521, (1)–(2):

```math
\mathbf E=-\gamma^2\boldsymbol\Pi+\mathop{\mathrm{grad}}\mathop{\mathrm{div}}\boldsymbol\Pi,
\qquad\text{(1)}
```

```math
\mathbf H=\frac{\gamma^2}{j\omega\mu}\mathop{\mathrm{curl}}\boldsymbol\Pi.
\qquad\text{(2)}
```

Continuation p. 1522, (3)–(4):

```math
E_x=-\gamma^2\Pi_x+
\frac{\partial}{\partial x}\left(
\frac{\partial\Pi_x}{\partial x}+\frac{\partial\Pi_y}{\partial y}
+\frac{\partial\Pi_z}{\partial z}\right)
=-IZ_{12}-\frac{\partial V}{\partial x},
\qquad\text{(3)}
```

```math
Z_{12}=\frac{\gamma^2\Pi_x}{I}.
\qquad\text{(4)}
```

``\boldsymbol\Pi`` is the Hertz vector, ``I`` the longitudinal source current, ``V`` the potential in the printed field decomposition. These field equations do not supply a separate potential-coefficient or admittance output.

The two-conductor integral parent (11), p. 1522, retains its **two separate integral terms** and nested integration:

```math
\begin{aligned}
Z_{12}={}&j\omega\frac{\mu_0}{4\pi}
\int_{-\infty}^{\infty}
\left\{\frac{\exp(-\gamma_0R_1)}{R_1}
-\frac{\exp(-\gamma_0R_2)}{R_2}\right\}\exp(-\gamma_0x)\,dx\\
&+j\omega\frac{\mu_0}{2\pi}\int_0^\infty
B_1\exp\{-\alpha_0(h_1+h_2)\}
\left\{\int_{-\infty}^{\infty}J_0(r\lambda)
\exp(-\gamma_0x)\,dx\right\}\lambda\,d\lambda.
\end{aligned}
\qquad\text{(11)}
```

The source defines ``r=\sqrt{x^2+y^2}``, with ``J_0`` the Bessel function of the first kind, order zero, on p. 1522. Appendix §8.2, p. 1527, makes the distances explicit:

```math
R_1=\sqrt{x^2+D_1^2},\qquad R_2=\sqrt{x^2+D_2^2},
```

with ``D_1,D_2`` as in the formula section. The final ``\gamma_0=jk`` prescription leads to (12)–(15), transcribed above.

### Original appendix conditions and spectral coefficient

Appendix §8.1, p. 1527, (18) prints the following conditions. In particular, the second row's derivative variables and the fourth row's last component are **not changed**:

```math
\begin{aligned}
\gamma_i^2\Pi_{ix}&=\gamma_{i+1}^2\Pi_{(i+1)x},\\
\frac{1}{\mu_i}\gamma_i^2\frac{\partial\Pi_{ix}}{\partial x}
&=\frac{1}{\mu_{i+1}}\gamma_{i+1}^2
\frac{\partial\Pi_{(i+1)x}}{\partial z},\\
\frac{1}{\mu_i}\gamma_i^2\Pi_{iz}
&=\frac{1}{\mu_{i+1}}\gamma_{i+1}^2\Pi_{(i+1)z},\\
\frac{\partial\Pi_{ix}}{\partial x}+\frac{\partial\Pi_{iz}}{\partial z}
&=\frac{\partial\Pi_{(i+1)x}}{\partial x}
+\frac{\partial\Pi_{(i+1)x}}{\partial z}.
\end{aligned}
\qquad\text{(18)}
```

The printed interface assignment is ``i=0`` at ``z=0``, ``i=1`` at ``z=-d_1``, ``i=2`` at ``z=-d_2``. These are preserved boundary statements, not a complete independent solution for every Hertz component.

Required coefficient of the parent (11), unnumbered lines following (19), p. 1527:

```math
B_1=\{(C_1+C_2)/\mu_0((\beta_0+\beta_1)C_1
+(\beta_0-\beta_1)C_2)\}.
```

This deliberately retains the original slash followed by a product/grouping. **Unresolved:** whether the parenthesized factor after ``\mu_0`` belongs inside the denominator. No stacked-fraction reinterpretation is substituted from the clear final (14).

```math
C_1=(\beta_1+\beta_2)(\beta_2+\beta_3)
+(\beta_1-\beta_2)(\beta_2-\beta_3)
\exp\{2\alpha_2(d_1-d_2)\},
```

```math
C_2=\bigl[(\beta_1-\beta_2)(\beta_2+\beta_3)
+(\beta_1+\beta_2)(\beta_2-\beta_3)
\exp\{2\alpha_2(d_1-d_2)\}\bigr]\exp(-2\alpha_1d_1),
```

```math
\alpha_i=\sqrt{\lambda^2+\gamma_i^2},\qquad
\beta_i=\alpha_i/\mu_i,\qquad
\gamma_i^2=-k_i^2=j\omega\mu_i(\sigma_i+j\omega\varepsilon_i).
```

Do not replace capital ``C`` with lower-case ``c`` or ``\alpha`` with ``a`` before the source's spectral transformation. The appendix dipole identity (20), p. 1527, is

```math
\int_0^\infty
\frac{J_0(r\lambda)\exp\{-\sqrt{\lambda^2+\gamma_0^2}(h+z)\}}
{\sqrt{\lambda^2+\gamma_0^2}}\lambda\,d\lambda
=\frac{\exp(-\gamma_0R_2)}{R_2},
\qquad R_2=\sqrt{r^2+(h+z)^2}.
\qquad\text{(20)}
```

Here ``h`` and ``z`` are the original dipole source/observation heights; ``h_1,h_2`` replace them in the two-conductor problem. No additional convergence/branch convention is printed with this identity.

### Original reduction in appendix §8.2

Printed p. 1528, (21):

```math
\int_{-\infty}^{\infty}
\left\{\frac{\exp[-jk\sqrt{x^2+D_1^2}]}{\sqrt{x^2+D_1^2}}
-\frac{\exp[-jk\sqrt{x^2+D_2^2}]}{\sqrt{x^2+D_2^2}}\right\}
\exp(-jkx)\,dx=W.
\qquad\text{(21)}
```

The intervening prose says to replace ``x+\sqrt{x^2+D^2}`` by ``t``; it does **not** include ``k`` in that verbal substitution. The subsequent equations nonetheless print:

```math
W=\lim_{s\to\infty}\left\{
\int_{q_1}^{\infty}\frac{\exp(-jt)}{t}\,dt
-\int_{q_2}^{\infty}\frac{\exp(-jt)}{t}\,dt\right\},
\qquad\text{(22)}
```

```math
q_1=k[\sqrt{s^2+D_1^2}-s],\qquad
q_2=k[\sqrt{s^2+D_2^2}-s],
```

```math
W=\lim_{s\to\infty}\left\{
\ln\left(\frac{q_2}{q_1}\right)
-\sum_{n=1}^{\infty}\frac{(-1)^n}{n\times n!}(q_1^n-q_2^n)
\right\}
=2\ln(D_2/D_1).
\qquad\text{(23)}
```

The coefficient is printed ``(-1)^n``, not ``(-j)^n``. Both the substitution and series issue are suspected published defects, not repaired identities. ``s`` is the limiting longitudinal coordinate here, **not** the 1975 spectral variable; ``t,q_1,q_2`` are the printed integration/end-point variables with unresolved substitution scaling.

The separate Bessel-transform identity (24) is

```math
\int_{-\infty}^{\infty}J_0(\sqrt{x^2+y^2}\lambda)\exp(-jkx)\,dx
=\begin{cases}
0,&\lambda<k,\\[2pt]
\displaystyle\frac{2\cos(y\sqrt{\lambda^2-k^2})}{\sqrt{\lambda^2-k^2}},
&\lambda>k.
\end{cases}
\qquad\text{(24)}
```

The source does not give a value at ``\lambda=k``. No endpoint value is manufactured. Equation (25) then prints

```math
Z_{12}=j\omega\frac{\mu_0}{2\pi}
\left(\ln\frac{D_2}{D_1}
+2\int_k^{\infty}B_1\exp\{-\alpha_0(h_1+h_2)\}
\frac{\cos[y\sqrt{\lambda^2-k^2}]}{\sqrt{\lambda^2-k^2}}
\lambda\,d\lambda\right).
\qquad\text{(25)}
```

The stated replacement ``\lambda^2-k^2=v^2`` produces the final unnumbered repetition of (12)–(14). This is the published reduction; it does not remove the annotated ambiguity in ``B_1`` or correct the preceding transform steps.

### Separately attributed 1975 restatement

Equations (1)–(2), pp. 500–501, repeat the external/earth decomposition with a different final spectral symbol:

```math
Z_{12}=j\omega(\mu_0/2\pi)\ln(D_2/D_1)+Z_e,
\qquad\text{(1975:1)}
```

```math
Z_e=j\omega(\mu_0/\pi)
\int_0^{\infty}A_3\exp\{-(h_1+h_2)s\}\cos(ys)\,ds,
\qquad\text{(1975:2)}
```

with ``D_1=\sqrt{y^2+(h_1-h_2)^2}``, ``D_2=\sqrt{y^2+(h_1+h_2)^2}``. Here and only here ``s`` is the inverse-length spectral variable. The 1975 page does not independently print the 1973 longitudinal exponential prescription.

Equation (3) and following definitions, p. 501:

```math
A_3=(c_1+c_2)/\{(s+\mu_0b_1)c_1+(s-\mu_0b_1)c_2\},
\qquad\text{(1975:3)}
```

```math
c_1=(b_1+b_2)(b_2+b_3)
+(b_1-b_2)(b_2-b_3)\exp\{2a_2(d_1-d_2)\}.
```

The printed ``c_2`` line visibly starts with ``[(b_1-b_2)(b_2+b_3)+(b_1+b_2)(b_2-b_3)`` and continues with ``\times\exp\{2a_2(d_1-d_2)`` followed by ``\exp(-2a_1d_1)``. **Unresolved source delimiter:** the closing square bracket/grouping before the last exponential cannot be verified in the inspected 500-dpi crop. This incomplete delimiter evidence is not silently completed with the clear 1973 ``c_2``. It blocks a completely resolved 1975 coefficient witness, not the independently transcribed 1973 kernel.

The constitutive lines themselves are clear:

```math
b_i=a_i/\mu_i,\qquad a_i=\sqrt{s^2+k_0^2-k_i^2},
\qquad i=1,2,3,
```

```math
k_0^2=-j\omega\mu_0(1/\rho_0+j\omega\varepsilon_0),
\qquad
k_i^2=-j\omega\mu_i(1/\rho_i+j\omega\mu_i).
```

The final term really prints **``j\omega\mu_i``**, not permittivity. The adjacent prose identifies ``\varepsilon,\mu,\rho`` as permittivity, permeability and resistivity. This is a suspected published constitutive/dimensional defect, retained rather than replaced by the 1973 bulk relation.

For ``d_2\to\infty``, the 1975 two-layer ``A_2`` equation (4) is split between the bottom left and top right columns of p. 501. Its first printed fragment is

```math
A_2=\frac{(b_1+b_2)+(b_1-b_2)\exp(-2a_1d_1)}
{(s+\mu_0b_1)(b_1+b_2)+(s-\mu_0b_1)(b_1-b_2)}\times
```

and the right-column continuation displays an asterisk alongside a fraction rule, ``\exp(-2a_1d_1)`` below that rule, no visible numerator, then label (4). **Unresolved layout interpretation:** the complete scope of the continued denominator cannot be asserted as a verified single fraction. The blank numerator is not changed into 1, nor is the expression silently replaced by original (15). This is a source-fragment transcription, not an executable formula.

The homogeneous limit ``d_1\to\infty`` and the correction decomposition are clear:

```math
A_1=(s+\mu_0b_1)^{-1},\qquad\text{(1975:5)}
```

```math
Z_e=(\omega\mu_0/2\pi)(P+jQ).\qquad\text{(1975:6)}
```

The authors state ``\varepsilon_i=\varepsilon_0`` for Fig. 2 and subsequent calculations (p. 501). Their stated near-negligibility below 1 MHz belongs to that study's justification, not a universal validity bound. Sections 3–5, pp. 502–506, study attenuation, propagation speed, step responses and switching surges using these earlier-attributed parameters. Their conditional homogeneous-earth response comparisons do not supply a new closed-form effective-medium impedance. No potential-coefficient/admittance kernel or arbitrary-layer recurrence occurs in the inspected seven pages.

## Notation map

No display renaming is performed. The following map distinguishes context-dependent reuse rather than treating matching glyphs as interchangeable inputs.

| Source symbol | Display symbol | Meaning / units / convention |
| --- | --- | --- |
| ``Z_{12},Z_e`` | Unchanged | External mutual impedance and earth-return correction respectively; Ω/m implied by (3)–(4). Neither is internal conductor impedance. |
| ``P,Q`` | Unchanged | Dimensionless correction quantities in (13)/(1975:6); not a potential-coefficient matrix. |
| ``I,V,\mathbf E,\mathbf H,\boldsymbol\Pi`` | Unchanged | Current A, potential V, electric field V/m, magnetic field A/m, Hertz vector (V·m implied by (1)); subscripts specify components/media. |
| ``x,y,z,h,h_1,h_2`` | Unchanged | Longitudinal/lateral/upward coordinates and source/target heights in m. ``h,z`` in dipole (20) become ``h_1,h_2`` in line interaction. |
| ``d_1,d_2`` | Unchanged | Depths from top interface in m; ``d_2`` is cumulative, not second-layer thickness. |
| ``D_1,D_2,r,R_1,R_2`` | Unchanged | Defined positive geometrical distances in m; own-image ``D_2=2h`` for self earth correction. No new radius substitution for ``D_1``. |
| ``\omega,j`` | Unchanged | Angular frequency rad/s and imaginary unit; printed harmonic factors retained. |
| ``\sigma_i,\rho_i,\varepsilon_i,\mu_i`` | Unchanged | Conductivity S/m, resistivity Ωm, permittivity F/m, permeability H/m. Independent scalar layers; 1975 prints a suspect permeability inside its conductivity-plus-displacement bracket. |
| ``\gamma,\gamma_0,\gamma_i,k,k_i`` | Unchanged | Bulk propagation/wavenumber quantities with source's separately prescribed longitudinal role. Inverse-length units for the consistent 1973 definitions. Main unindexed ``k`` is not replaced by a guessed independent input. |
| ``\lambda,v`` | Unchanged | 1973 original/final spectral variables, inverse m; source transform ``\lambda^2-k^2=v^2``. |
| ``\alpha_i,\beta_i,C_1,C_2,B_1`` | Unchanged | 1973 pre-transform factors. ``\alpha_i`` inverse m, ``\beta_i`` inverse H, ``C_1,C_2`` inverse H²; intended units of ambiguous ``B_1`` cannot select its grouping. |
| ``a_i,b_i,c_1,c_2,B_2,B_3`` | Unchanged | 1973 final factors, inverse m, inverse H, inverse H², and length for the clear coefficient fractions. ``B_3`` is the two-layer formula, despite its subscript. |
| ``J_0`` | Unchanged | First-kind Bessel function, order zero; dimensionless function of the printed argument. |
| ``W,s,t,q_1,q_2,n`` in 1973 appendix | Unchanged | ``W`` dimensionless transform integral; ``s`` limiting length coordinate; ``t,q_1,q_2`` printed transformed variables with scaling inconsistency; ``n`` positive integer and series extends to infinity. |
| ``w,\theta,\alpha_1,\alpha_2`` in numerical section | Unchanged | ``w=D_2\sqrt{\omega\mu_0/\rho_1}``; ``\theta`` image angle; numerical ``\alpha_1=\rho_1/\rho_2,\alpha_2=\rho_1/\rho_3`` are dimensionless ratios, not the appendix's spectral roots. |
| ``s,A_3,A_2,A_1`` in 1975 | Unchanged | 1975 spectral variable and three-/two-/one-layer coefficients; not the 1973 appendix limit coordinate or unrelated snow ``A`` kernels. |

## Evidence and approximation sources

The main 1973 fraction, all lower-case coefficients, material terms and stated layer limits were checked on original pp. 1522–1523. Appendix pp. 1527–1528 supplies the original spectral definitions and declared transform chain; the source's separate integrals remain separate. The prerequisite longitudinal approximation is explicitly introduced in §2.1 and imposed in §2.3, rather than inferred from an implementation's propagation route.

Section 3 defines ``w=D_2\sqrt{\omega\mu_0/\rho_1}``, ``\alpha_1=\rho_1/\rho_2``, ``\alpha_2=\rho_1/\rho_3``. Figs. 3–8 test roughly ``0.1\le w\le10``, heights 10 or 25 m and selected layer resistivity ratios. These test intervals do not establish universal accuracy. Section 3.1 says the plane-wave assumption may explain difficult low-``w`` behaviour and calls for further investigation. Section 4's modified FFT is a response-evaluation method, not a new physical impedance or admittance formulation, and is not implemented or registered here.

For the homogeneous/two-layer comparisons, p. 1523 first names Carson, then separately discusses Sunde, Wedepohl–Wasley, Iwamoto and Wise. Its Iwamoto agreement is expressly restricted to earth permittivities/permeabilities equal to air, whereas the Wise homogeneous comparison retains displacement current. These attributed comparisons are not independent verification of unseen originals or a claim that arbitrary ``\varepsilon_i,\mu_i`` is Carson's original model.

The earlier-source chain is preserved rather than closed by assumption:

- Iwamoto's reference [8], *Use of the travelling waves on the measurement of earth resistivity*, 78, pp. 1038–1049 (1958), is available from the [original publisher](https://www.jstage.jst.go.jp/article/ieejjournal1888/78/839/78_839_1038/_pdf/-char/ja). The first page verifies the Japanese title *進行波による大地固有抵抗の測定法*, K. Iwamoto/岩本国三, August 1958, vol. 78, issue 839. Its English title footnote spells the words **Traveling** and **Measurment**; these source words differ from the 1973 reference. The [separate secondary-witness record](../../1958/two-layer-earth-overhead-self-integral/Iwamoto1958b.md) contains its two-layer self integral, finite-depth expression, and normalization. Iwamoto credits the same form to Sunde p. 119 (1949). The inspected Sunde source is the corrected 1968 edition, and algebraic equivalence to the 1973 coefficient has not been established.
- Wedepohl–Wasley (1966) is not used as equation evidence. Sunde's 1968 second edition is not equated to a different edition without page evidence. Wise and Carson retain their own source identities.

The 1975 paper is separately classified as application-only **with a useful later-author transcription**, because it explicitly credits the 1973 sources for (1)–(3) and applies them to surge studies. Its main variable/coefficient renaming does not create a new contribution. Its ambiguous layout and suspect constitutive term are not scientifically validated modifications. Although the [snow potential-coefficient paper](../../../external-admittance/2000/snow-layer-overhead-potential-coefficient/Ametani2001.md) cites 1975 for dipole/interface methodology, the inspected 1975 article does not print the snow admittance/potential kernel. Method citation alone cannot backdate that later expression.

## Limitations and discrepancies

1. **Original model restriction, stated:** ``\gamma_0=jk`` and plane-wave longitudinal dependence are imposed. The three-layer final integral is not an arbitrary-layer recursion or an independently solved dispersion model. Its self earth-correction geometry is supported; a general finite-radius full-external self logarithm is not supplied.
2. **Original appendix grouping, unresolved:** ``B_1`` is printed with slash/product adjacency, unlike the clear final ``B_2`` fraction. The record does not choose a denominator scope by expected dimensions or by back-substitution.
3. **Original published field defects, suspected:** (18) differentiates ``\Pi_{ix}`` with respect to ``x`` on the left and ``\Pi_{(i+1)x}`` with respect to ``z`` on the right; its last row ends in ``\partial\Pi_{(i+1)x}/\partial z`` rather than the ``z`` component used on the left. Those tokens are image-confirmed. The pre-integration Hertz expression (10), p. 1522, also prints ``\exp(-\gamma_0x)`` both outside and inside the longitudinal integral and an unindexed ``\gamma^2`` denominator, whereas (11) is printed separately. This is not silently reconciled.
4. **Original published transform defects, suspected:** the substitution prose preceding (22) omits the ``k`` scaling used in the subsequent endpoints/exponential, and (23) prints ``(-1)^n`` for a series following ``\exp(-jt)``. The transform chain is retained as published, including its final result; no repaired coefficients are supplied.
5. **Original missing definitions:** root branch and the exact endpoint interpretation at ``\lambda=k`` are not printed. No complex-material generality beyond the stated scalar constitutive model, convergence extension or finite-radius regularization is inferred.
6. **Later 1975 defects/gaps:** the clear ``k_i^2`` line has ``j\omega\mu_i`` in its inner bracket. The ``c_2`` closing delimiter and the cross-column denominator scope in (4) are unresolved in the inspected scan. These are explicitly incomplete later witnesses; the readable 1973 expression is not presented as their corrected text.
7. **LCM/source identity disagreement, directly observed:** the existing `:Nakagawa1973` route is described as recursive ``N``-layer and returns zero impressed propagation; the cited original's eight pages provide neither that recurrence nor a zero prescription. This flags a later comparison task, not proof of the numerical implementation's correctness or incorrectness. LCM files are untouched.
8. **Citation/history limits:** Iwamoto's August publication supplies a Sunde-attributed secondary kernel witness. His January paper supplies a distinct logarithmic-integral evaluation and documents the operator discrepancies separately. The inspected Sunde source is the corrected 1968 edition.
