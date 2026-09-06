# Overhead-conductor potential coefficient with a finite snow layer

## Identity and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite parallel overhead conductors above snow; heights ``h_i,h_j`` and lateral spacing ``y``. The source gives no general finite-radius self substitution. |
| Calculated quantities | Parallel-conductor mutual coefficient, snow/earth correction, no-snow reduction, and contextual capacitance comparisons. |
| Earth structure | Air above a finite snow layer of thickness ``d`` and an earth half-space. |
| Model and approximation | Integral representation with ``\Gamma=\gamma_0``. The no-snow result is a parameter substitution. Equations (20)–(21) use an approximate two-capacitance interpretation and do not replace kernel (15). |
| Main source | A. Ametani, N. Nagaoka, and R. Koide (2001 English translation; 2000 Japanese original). |
| Citation key(s) | Primary English publication: `:Ametani2001`; Japanese original: `:Ametani2000` |
| Evidence status | Both publications checked against page images. The ``A_2`` grouping, constitutive mapping, branch, and general self prescription remain unresolved. |

**Description.** Potential coefficient between infinitely long parallel overhead conductors in air above a finite snow layer and an imperfectly conducting earth half-space, with a logarithmic geometric term and a spectral snow/earth correction.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Prescribed ``I=I_0\exp(-\gamma_0x)``: the imposed longitudinal constant is the air bulk constant, not an unknown dispersion solution and not zero. | Stated — JP p. 272/PDF 2 and EN p. 28/PDF 3 immediately before (10). |
| Air propagation constant ``γ_air`` | Source ``\gamma_0=j\omega\sqrt{\varepsilon_0\mu_0}`` after setting ``\sigma_0=0``; retained in the final ``a_i`` definitions. | Stated — JP p. 272/EN p. 28 immediately before (14); definitions after (7). |
| Earth propagation constant ``γ_earth`` | Source layer 2: ``\gamma_2^2=j\omega\mu_2(\sigma_2+j\omega\varepsilon_2)``. Finite snow is layer 1 with the same indexed definition. ``a_i=\sqrt{s^2-\gamma_0^2+\gamma_i^2}``, ``i=1,2``. Root branch not stated. | Stated — JP p. 272/EN p. 27 definitions after (5) and (7); Not stated — branch. |
| Earth permittivity and displacement current | ``\varepsilon_2`` and ``j\omega\varepsilon_2`` retained alongside ``\sigma_2``. Snow has a separately printed two-relaxation complex-permittivity relation (17), including ``\sigma_1/(j\omega)``; mapping that entire complex quantity into ``\varepsilon_1`` while keeping separate ``\sigma_1`` is unresolved, not an authorized double-counting prescription. No-snow formula sets ``\varepsilon_1=\varepsilon_0,\sigma_1=0``. | Equation-implied — indexed bulk definition, JP p. 272/EN p. 27, retains both terms. Stated — (17) and conditions before (16), JP p. 273/EN p. 28. Unresolved — relative/absolute permittivity and loss bookkeeping between these contexts. |
| Range of validity | Integral representation under the printed infinite-line and prescribed-propagation model. Authors state application to two-layer insulating surroundings beyond snow; no arbitrary stratification or quantified universal bound is supplied. Example plots span 10 Hz–1 MHz; snow depths 1–5 m and Table 1 states are tests, not proved limits. | Stated — JP p. 273 opening / EN p. 28 after (15); JP pp. 273–276/EN pp. 28–32 examples. Not stated — universal frequency/electrical-size bound. |
| Earth permeability ``μ_earth`` | ``\mu_2`` is retained independently in coefficients and in (16), not replaced with ``\mu_0``. General finite-layer ``\mu_1`` is retained; physical snow has relative permeability unity. No-snow formula prescribes ``\mu_1=\mu_0``. | Equation-implied — coefficient definitions and (16). Stated — JP p. 273 §3 explicitly says relative snow permeability; EN p. 28 §3 says permeability unity. |
| Arrangement | Both conductors above the air/snow interface; mutual coefficient uses distinct conductor/image distances. Single-wire and multiconductor self responses are evaluated in examples, but no general radius substitution for the integral self formula is stated. No underground or mixed source/target pair is supplied. | Stated — JP Figs. 1–3, (8), §§4–5 / EN Figs. 1–3, (8), §§4–5. Unresolved — full self regularization. |
| Earth structure | Air 0 above snow 1 of thickness ``d``, above earth half-space 2; parallel flat interfaces. The source's “two-layer insulator” counts air and snow, not two finite soil layers. | Stated — JP p. 272/EN p. 27 Figs. 2–3. Equation-implied — one finite-depth exponential ``\exp(-a_1d)`` and terminal layer-2 coefficient. |
| Conductor and insulation geometry | Infinite parallel conductors along ``x``; heights ``h_i,h_j`` measured above snow, lateral spacing ``y``. No conductor-attached insulation shell or internal surface term is derived. Single-wire example has ``h+d=25\,\mathrm m`` and ``a=2\,\mathrm{cm}``; its radius is not itself a source-prescribed substitution in (15). | Stated — JP p. 271 §2, p. 272 Fig. 3, p. 273 §4 / EN pp. 26–28. |
| Constitutive and field assumptions | Scalar medium parameters, homogeneous within each depicted layer; Hertz-vector line construction with ``\Pi_y=0``. Isotropic scalar treatment is equation-implied, not a general anisotropic result. Snow influence on series earth impedance is neglected and Carson is used separately; this does not set every ``\sigma_1`` or dielectric-loss term to zero in the admittance model. No internal skin/proximity formula is supplied here. | Stated — JP p. 271 §§1–2 / EN p. 26 §§1–2. Equation-implied — scalar definitions and flat-layer construction, JP p. 272/EN p. 27. |
| Conventions | ``x`` longitudinal; ``z`` upward from the air/snow surface, ``y`` horizontal. Explicit longitudinal factor ``\exp(-\gamma_0x)``. ``j\omega`` corresponds to the source's frequency-domain time derivative, but an explicit full time exponential is not printed. ``V=-\operatorname{div}\boldsymbol\Pi``, ``Q_0`` charge per unit length, ``V=Q_0P_{ij}``; ``P`` consequently has units m/F. No additional voltage-reference transformation is introduced. | Stated — JP Fig. 1 and (9)–(12) / EN Fig. 1 and (9)–(12). Equation-implied — (1)–(2) time derivative and (11) units. Not stated — explicit time exponential or root branch. |

**Expression.** The common main output is ``P_{ij}``, JP p. 272/PDF 2 and EN p. 28/PDF 3, (14)–(15). It includes the printed logarithmic term; ``M+jN`` alone is the snow/earth correction.

```math
P_{ij}=\frac{1}{2\pi\varepsilon_0}
\left\{\ln\left(\frac{D_{ij}}{d_{ij}}\right)+(M+jN)\right\},
\qquad
M+jN=2\int_0^\infty(A_1-sA_2)
\exp\{-(h_i+h_j)s\}\cos(ys)\,ds.
```

JP coefficient witness, p. 272/PDF 2 below (7). **The three-line slash/product layout of ``A_2`` is retained; whether the last factor belongs to the denominator is unresolved.** The line break is evidentiary, not a selected computational grouping.

```math
A_1=\frac{c_1+c_2}{(s+\mu_0b_1)c_1+(s-\mu_0b_1)c_2},
```

```math
\begin{aligned}
A_2={}&\{4b_1c_0c_5(1-\tau_2^2)
 +(c_1+c_2)(c_3-c_4)(1-\tau_1^2)\}\\
&/\{2\mu_1c_4c_5+(c_3-c_4)(\mu_1c_5+\mu_0s)\}\\
&\cdot\{(s/\mu_0+b_1)c_1+(s/\mu_0-b_1)c_2\}.
\end{aligned}
```

```math
\begin{gathered}
c_0=\exp(-a_1d),\qquad c_1=b_1+b_2,\qquad
c_2=(b_1-b_2)\exp(-2a_1d),\\
c_3=(c_6+a_1)/c_0,\qquad c_4=(c_6-a_1)c_0,\\
c_5=a_1\tau_1^2,\qquad c_6=\mu_2a_2\tau_2^2/\mu_1,
\qquad \tau_i^2=\gamma_{i-1}^2/\gamma_i^2,\\
b_i=a_i/\mu_i,\qquad a_i=\sqrt{s^2-\gamma_0^2+\gamma_i^2},
\qquad i=1,2,\\
\gamma_i^2=j\omega\mu_i(\sigma_i+j\omega\varepsilon_i),\qquad i=0,1,2,
\qquad \sigma_0=0,\quad\gamma_0=j\omega\sqrt{\varepsilon_0\mu_0},\\
d_{ij}=\sqrt{y^2+(h_i-h_j)^2},\qquad
D_{ij}=\sqrt{y^2+(h_i+h_j)^2}.
\end{gathered}
```

EN coefficient witness, p. 27/PDF 2 below (7), **separate from JP**: the second ``A_1`` factor uses ``b_2``, ``A_2`` explicitly divides by ``A_3A_4``, and the ratio definition has an unsquared left-hand ``\tau_i``. These are not substituted into the JP witness.

```math
\begin{gathered}
A_1=\frac{c_1+c_2}{(s+\mu_0b_1)c_1+(s-\mu_0b_2)c_2},\\
A_2=\frac{4b_1c_0c_5(1-\tau_2^2)
 +(c_1+c_2)(c_3-c_4)(1-\tau_1^2)}{A_3\cdot A_4},\\
A_3=\{2\mu_1c_4c_5+(c_3-c_4)(\mu_1c_5+\mu_0s)\},\\
A_4=\{(s/\mu_0+b_1)c_1+(s/\mu_0-b_1)c_2\},\\
c_0=\exp(-a_1d),\qquad c_1=b_1+b_2,\qquad c_2=(b_1-b_2)c_0^2,\\
c_3=(c_6+a_1)/c_0,\qquad c_4=(c_6-a_1)c_0,\\
c_5=a_1\tau_1^2,\qquad c_6=\mu_2a_2\tau_2^2/\mu_1,
\qquad \tau_i=\gamma_{i-1}^2/\gamma_i^2,\\
b_i=a_i/\mu_i,\qquad a_i=\sqrt{s^2-\gamma_0^2+\gamma_i^2},\qquad i=1,2.
\end{gathered}
```

The EN bulk-medium, air and distance definitions are the same expressions printed in the preceding JP block. In both versions ``s`` is the final real integration variable from zero to infinity; lengths are in metres, ``\mu_i`` in H/m, ``\sigma_i`` in S/m and the bulk definition requires an absolute ``\varepsilon_i`` in F/m. The separate snow-data convention below is unresolved. Square-root branches are not supplied. Both coefficient paragraphs confusingly also describe ``\tau_1,\tau_2`` as snow-permittivity parameters; the dimensionless interface ratios cannot silently be identified with the relaxation times in microseconds in Table 1.

Source-stated no-snow formula, JP p. 273/PDF 3 and EN p. 28/PDF 3, (16), with ``d=0,\mu_1=\mu_0,\varepsilon_1=\varepsilon_0,\sigma_1=0``:

```math
M+jN=2\int_0^\infty
\frac{(s+a_2\mu_2/\mu_0)\exp\{-(h_i+h_j)s\}\cos(ys)}
{(s+a_2\mu_0/\mu_2)(s/\tau_2^2+a_2\mu_2/\mu_0)}\,ds.
```

Each witness retains its own printed interface-ratio definition. The authors attribute agreement of (16) to Nakagawa [6]; this is not independent original verification of Nakagawa from the snow paper.

Source-provided admittance operation in the **perfect-earth, no-snow comparison context**, JP p. 273/EN p. 28, (18):

```math
[Y]=j\omega[P]^{-1},\qquad
P_{ij}=\frac{1}{2\pi\varepsilon_0}\ln\left(\frac{D_{ij}}{d_{ij}}\right).
```

This preserves the printed matrix inverse and the particular adjacent definition of ``P``. The paper computes snow admittances using its potential coefficient, but does not print a separate fully assembled finite-radius snow matrix here. This record does not turn ``P_{ij}^{-1}`` into a scalar reciprocal admittance or insert a new full-matrix formula attributed to (18).

Snow constitutive dependency, JP p. 273/PDF 3, (17):

```math
\varepsilon_s=\varepsilon_{s\infty}
 +\frac{C_1}{1+j\omega\tau_1}
 +\frac{C_2}{1+j\omega\tau_2}
 +\frac{\sigma_1}{j\omega}
 =\varepsilon'-j\varepsilon''.
```

EN p. 28/PDF 3, (17), instead ends in ``\varepsilon'-j\omega\varepsilon''``; the extra ``\omega`` is printed and is retained as a disagreement. Both versions use the same rational terms and Table 1 values. ``\tau_1,\tau_2`` in **this constitutive context** are relaxation times; ``C_1,C_2`` are the dielectric strengths, not the physical layer capacitances of (20). Their original determination procedure is not supplied in the snow paper.

Table 1, JP p. 273/EN p. 28, as printed. The conductivity column has the scale ``10^{-9}``; strengths are headed **pF**, not pF/m. No silent unit conversion or missing ``\varepsilon_0`` is supplied.

| Case | Density (g/cm³) | T (°C) | ``\varepsilon_{s\infty}`` | ``\tau_1`` (µs) | ``\tau_2`` (µs) | ``\sigma`` (S/m), ×10⁻⁹ | ``C_1`` (pF) | ``C_2`` (pF) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 0.0768 | −14 | 1.17 | 13.2 | 413 | 13 | 9.37 | 31.8 |
| 2 | 0.127 | −24 | 1.35 | 9.6 | 434 | 15 | 14.7 | 44.3 |
| 3 | 0.127 | −14 | 1.35 | 11.3 | 420 | 27 | 19.0 | 97.9 |
| 4 | 0.127 | −11 | 1.58 | 11.4 | 384 | 39 | 24.6 | 82.2 |
| 5 | 0.127 | −6 | 1.57 | 12.1 | 466 | 89 | 24.1 | 122. |
| 6 | 0.317 | −11 | 1.84 | 9.0 | 417 | 10 | 72.8 | 91.6 |

**Approximation.** The principal integral is a representation within the prescribed ``\Gamma=\gamma_0`` infinite-line model, not a fitted closed-form approximation; no expansion order for that propagation prescription is supplied. The no-snow reduction is a source-stated parameter substitution, not a series truncation. For interpretation, JP p. 274/EN p. 30 replaces the geometry by two series-connected layer capacitances: ``C=C_1C_2/(C_1+C_2)`` (20), then ``C\simeq C_1`` (21) under ``C_1\ll C_2``. JP explicitly describes an approximate parallel-plate replacement. No expansion series, retained numerical order, or remainder bound is given; it is not an algebraic replacement for kernel (15).

**Limitations.** Original ``A_2`` grouping, unprovided general self prescription and branches, interface/relaxation ``\tau`` reuse, and the relative/absolute snow-permittivity and conductivity mapping prevent an unambiguous numerical interpretation. EN differs from JP in ``A_1``, ``\tau_i`` definition, Hertz-field sign, the printed integration transformation and the imaginary-permittivity convention. These are source discrepancies, not conversion repairs. No arbitrary-layer recurrence, mixed/buried pair, new conductor internal impedance or new snow-modified series-earth kernel is supplied.

**Reference.** [Ametani2001](@cite), published English translation, pp. 26–30, equations (7)–(18) and Table 1; [Ametani2000](@cite), Japanese original, pp. 271–274, with the same equation labels and table.

**Transcription source.** Original JP and published EN page images inspected separately throughout. Main kernel, coefficients, constitutive dependency, comparisons and the particular conflicting tokens were checked at their stated locators, including magnified coefficient/permittivity crops. An image-verified ambiguity is not a resolved dependency. The English translation is not an erratum. No unverified Markdown enrichment is used as mathematical authority.

## Source transcription

### Hertz convention and imposed longitudinal dependence

JP p. 271/PDF 1, (1)–(2); EN p. 26/PDF 1, same labels. These are convention evidence for the retained line result, not independent harvested formulations. Common definition:

```math
\mathbf A=\mu\left(\sigma+\varepsilon\frac{\partial}{\partial t}\right)
\boldsymbol\Pi.\qquad\text{(1)}
```

JP:

```math
\mathbf E=-\gamma^2\boldsymbol\Pi+\operatorname{grad}\operatorname{div}\boldsymbol\Pi,
\qquad
\mathbf H=\frac{\gamma^2}{j\omega\mu}\operatorname{rot}\boldsymbol\Pi.
\qquad\text{(2)}
```

EN:

```math
\mathbf E=\gamma^2\boldsymbol\Pi+\operatorname{grad}\operatorname{div}\boldsymbol\Pi,
\qquad
\mathbf H=\frac{\gamma^2}{j\omega\mu}\operatorname{rot}\boldsymbol\Pi.
\qquad\text{(2)}
```

Both then define ``\gamma^2=j\omega\mu(\sigma+j\omega\varepsilon)``. The material-specific version below (5), JP p. 272/EN p. 27, is

```math
\gamma_i^2=j\omega\mu_i(\sigma_i+j\omega\varepsilon_i),\qquad i=0,1,2.
```

The source obtains the infinite-line result by solving the dipole interface problem (3)–(6) and integrating longitudinally from ``-\infty`` to ``+\infty``. The explanatory transformation before (7) is printed ``\lambda^2+\gamma^2=s^2`` in JP (unindexed ``\gamma``), but ``\lambda^2+r^2=s^2`` in EN, where EN earlier defined ``r=\sqrt{x^2+y^2}`` as a radial distance. These tokens are recorded as a conflicting derivation statement; no index or substitute transformation is invented. The final coefficient dependencies below are independently printed and do not require reconstructing the unspecified dipole amplitudes ``f,g,F,G``.

### Infinite-line vector and interface coefficients

JP p. 272/PDF 2 and EN p. 27/PDF 2, (7):

```math
\begin{aligned}
\Pi_{0x}&=\frac{j\omega\mu_0 I}{2\pi\gamma_0^2}
\left[\ln\left(\frac{D_{ij}}{d_{ij}}\right)
+2\int_0^\infty A_1\exp\{-(h_i+h_j)s\}\cos(ys)\,ds\right],\\
\Pi_{0z}&=\frac{j\omega\mu_0 I}{2\pi\gamma_0^2}
\int_0^\infty A_2\exp\{-(z+h_j)s\}\cos(ys)\,ds.
\end{aligned}
\qquad\text{(7)}
```

The two integrals remain separate, with a factor 2 in the ``\Pi_{0x}`` integral and no added factor 2 in ``\Pi_{0z}``. JP definitions immediately following (7), in their printed order:

```math
A_1=(c_1+c_2)/\{(s+\mu_0b_1)c_1+(s-\mu_0b_1)c_2\},
```

```math
\begin{aligned}
A_2={}&\{4b_1c_0c_5(1-\tau_2^2)+(c_1+c_2)(c_3-c_4)(1-\tau_1^2)\}\\
&/\{2\mu_1c_4c_5+(c_3-c_4)(\mu_1c_5+\mu_0s)\}\\
&\cdot\{(s/\mu_0+b_1)c_1+(s/\mu_0-b_1)c_2\},
\end{aligned}
```

```math
\begin{gathered}
c_0=\exp(-a_1d),\quad c_1=b_1+b_2,\quad c_2=(b_1-b_2)\exp(-2a_1d),\\
c_3=(c_6+a_1)/c_0,\quad c_4=(c_6-a_1)c_0,\\
c_5=a_1\tau_1^2,\quad c_6=\mu_2a_2\tau_2^2/\mu_1,\quad
\tau_i^2=\gamma_{i-1}^2/\gamma_i^2,\\
b_i=a_i/\mu_i,\quad a_i=\sqrt{s^2-\gamma_0^2+\gamma_i^2},\quad i=1,2.
\end{gathered}
```

The ``A_2`` slash/product scope remains unresolved; the original does not name ``A_3,A_4``. EN's own complete coefficient witness, not a correction applied to JP, is:

```math
\begin{gathered}
A_1=(c_1+c_2)/\{(s+\mu_0b_1)c_1+(s-\mu_0b_2)c_2\},\\
A_2=\{4b_1c_0c_5(1-\tau_2^2)+(c_1+c_2)(c_3-c_4)(1-\tau_1^2)\}/(A_3\cdot A_4),\\
A_3=\{2\mu_1c_4c_5+(c_3-c_4)(\mu_1c_5+\mu_0s)\},\\
A_4=\{(s/\mu_0+b_1)c_1+(s/\mu_0-b_1)c_2\},\\
c_0=\exp(-a_1d),\quad c_1=b_1+b_2,\quad c_2=(b_1-b_2)c_0^2,\\
c_3=(c_6+a_1)/c_0,\quad c_4=(c_6-a_1)c_0,\\
c_5=a_1\tau_1^2,\quad c_6=\mu_2a_2\tau_2^2/\mu_1,\quad
\tau_i=\gamma_{i-1}^2/\gamma_i^2,\\
b_i=a_i/\mu_i,\quad a_i=\sqrt{s^2-\gamma_0^2+\gamma_i^2},\quad i=1,2.
\end{gathered}
```

Both coefficient paragraphs also call ``\tau_1,\tau_2`` parameters of snow complex permittivity; that wording conflicts with treating the interface definitions as the microsecond relaxation times in (17). It is retained as a printed symbol collision.

### Potential coefficient and no-snow formula

JP p. 272/PDF 2; EN p. 27/PDF 2, distances (8), then p. 28/PDF 3, (9)–(15):

```math
d_{ij}=\sqrt{y^2+(h_i-h_j)^2},\quad
D_{ij}=\sqrt{y^2+(h_i+h_j)^2}.\qquad\text{(8)}
```

```math
V=-\operatorname{div}\boldsymbol\Pi
=-\left(\frac{\partial\Pi_x}{\partial x}
+\frac{\partial\Pi_y}{\partial y}
+\frac{\partial\Pi_z}{\partial z}\right).\qquad\text{(9)}
```

Immediately before (10), both witnesses prescribe ``I=I_0\exp(-\gamma_0x)``. They then print:

```math
V=\frac{j\omega\mu_0I}{2\pi\gamma_0}
\left[\ln\left(\frac{D_{ij}}{d_{ij}}\right)
+2\int_0^\infty(A_1-sA_2)\exp\{-(h_i+h_j)s\}\cos(ys)\,ds\right].
\qquad\text{(10)}
```

```math
V=Q_0P_{ij},\qquad\text{(11)}
\qquad I=j\omega Q_0/\gamma_0,\qquad\text{(12)}
\qquad P_{ij}=j\omega V/(\gamma_0 I).\qquad\text{(13)}
```

``Q_0`` is explicitly charge per unit conductor length. With source-stated ``\sigma_0=0`` and ``\gamma_0=j\omega\sqrt{\varepsilon_0\mu_0}``:

```math
P_{ij}=\frac{1}{2\pi\varepsilon_0}
\left\{\ln(D_{ij}/d_{ij})+(M+jN)\right\},\qquad\text{(14)}
```

```math
M+jN=2\int_0^\infty(A_1-sA_2)\exp\{-(h_i+h_j)s\}\cos(ys)\,ds.
\qquad\text{(15)}
```

The first term is identified as the geometric perfect-earth/no-snow potential coefficient; the second accounts for snow and imperfect earth. The coefficients must be read from **one selected published witness**, with its unresolved definitions, not assembled from favorable tokens in both.

JP p. 273/PDF 3 and EN p. 28/PDF 3 impose ``d=0,\mu_1=\mu_0,\varepsilon_1=\varepsilon_0,\sigma_1=0`` and give:

```math
M+jN=2\int_0^\infty
\frac{(s+a_2\mu_2/\mu_0)\exp\{-(h_i+h_j)s\}\cos(ys)\,ds}
{(s+a_2\mu_0/\mu_2)(s/\tau_2^2+a_2\mu_2/\mu_0)}.
\qquad\text{(16)}
```

All permeability ratios and the ``\tau_2^2`` power are retained. Neither version explicitly supplies the general integral's finite-radius self substitution. The later Papadopoulos self prescription must not be imported to complete this gap.

### Snow dielectric dependency and comparison assemblies

JP p. 273/PDF 3 and EN p. 28/PDF 3, §3, (17), respective witnesses:

```math
\begin{aligned}
\text{JP: }\quad\varepsilon_s
&=\varepsilon_{s\infty}+C_1/(1+j\omega\tau_1)
+C_2/(1+j\omega\tau_2)+\sigma_1/(j\omega)
=\varepsilon'-j\varepsilon'',\\
\text{EN: }\quad\varepsilon_s
&=\varepsilon_{s\infty}+C_1/(1+j\omega\tau_1)
+C_2/(1+j\omega\tau_2)+\sigma_1/(j\omega)
=\varepsilon'-j\omega\varepsilon''.
\end{aligned}
\qquad\text{(17)}
```

The complete Table 1 is transcribed at its source locator and applies to both witnesses. Its strengths and time constants are context-specific; the source's relative-permittivity terminology and printed pF strengths are retained instead of being changed to an absolute F/m model. The bulk ``\varepsilon_1,\sigma_1`` input mapping is not specified.

JP p. 273/EN p. 28, (18), perfect-earth/no-snow comparison:

```math
[Y]=j\omega[P]^{-1},\qquad
P_{ij}=\frac{1}{2\pi\varepsilon_0}\ln(D_{ij}/d_{ij}).\qquad\text{(18)}
```

EN alone, p. 29/PDF 4, unnumbered formula in the left column following the explanation of Fig. 5's ideal curve:

```math
C_i=2\pi\varepsilon_0/\ln(2h/r)=7.11\;[\mathrm{pF/m}].
```

Here ``r`` is used in the conductor-radius position of the ideal self-capacitance comparison, not the earlier dipole radial coordinate. The preceding example defines radius as ``a=2\,\mathrm{cm}``; the local ``a/r`` notation change does not define a general self rule for (15). This unnumbered formula appears only in the English publication, not on the corresponding Japanese p. 274.

JP p. 274/PDF 4 and EN p. 29/PDF 4, (19):

```math
\begin{aligned}
\text{JP: }\quad\mathrm{DY}&=(Y_s-Y_0)/Y_0\times100
\simeq(C_s-C_0)/C_0\times100\qquad(\%),\\
\text{EN: }\quad\mathrm{DY}&=(Y_s-Y_0)\times100/Y_0
\simeq(C_s-C_0)\times100/C_0\qquad[\%].
\end{aligned}
\qquad\text{(19)}
```

``Y_0,C_0`` refer to no snow; ``Y_s,C_s`` to snow. The capacitance approximation is explained by conductance being much smaller than susceptance in the illustrated case, not a replacement for all complex admittances. The source's depth trend in §4.1 (JP p. 274; EN pp. 29–30) is ``\mathrm{DY}=k\cdot d``, with ``k\simeq0.6\%/\mathrm m`` below 10 kHz and ``k\simeq0.2\%/\mathrm m`` above 100 kHz, for Table 1 case 5 and Fig. 6's 1–5 m depths. This observed near-proportionality has no printed fit residual or universal error bound and is not a new stand-alone kernel.

JP p. 274/PDF 4 and EN p. 30/PDF 5, §4.2, (20)–(21):

```math
C=C_1C_2/(C_1+C_2),\qquad\text{(20)}
\qquad C\simeq C_1.\qquad\text{(21)}
```

JP introduces approximate replacement by air/snow parallel-plate capacitances and ``C_1\ll C_2``; EN prints ``C_2\gg C_1``. In this context ``C_1`` is air-layer capacitance and ``C_2`` snow-layer capacitance, **not** the strengths named identically in (17). No geometric formula for these equivalent plate capacitances or a correction series is supplied here. Sections 5–6 contain modal/characteristic/step-response applications and conclusions, not additional line-parameter kernels.

## Notation map

Source notation is unchanged in the displayed mathematics. JP/EN labels identify witnesses; they do not rename physical quantities. The assumptions field ``\Gamma`` names the role filled by source ``\gamma_0`` in its prescribed longitudinal factor. Context-specific reuse is not a one-to-one identification across contexts.

| Source symbol | Display symbol | Meaning and units | Convention/evidence |
| --- | --- | --- | --- |
| ``P_{ij}``, ``[P]`` | unchanged | Potential coefficient, m/F; coefficient matrix | (11), (14); bracketed matrix in comparison (18) |
| ``M+jN`` | unchanged | Dimensionless snow/earth correction inside (14) | Combined complex quantity defined by (15); not independently declared resistance/reactance |
| ``[Y]``, ``Y_s,Y_0`` | unchanged | Per-length admittance, S/m | Matrix (18); snow/no-snow comparison (19) |
| ``I,I_0,Q_0,V`` | unchanged | Current and its amplitude (A), charge per length (C/m), scalar potential (V) | Prescribed negative longitudinal exponent; (9)–(13) |
| ``\gamma_i``, ``\gamma`` | unchanged | Bulk propagation constants, m⁻¹; unindexed generic constant in (2) | Layers 0 air, 1 snow, 2 earth; unindexed transformation before (7) remains unresolved |
| ``\mu_i,\sigma_i,\varepsilon_i`` | unchanged | Layer permeability H/m, conductivity S/m, bulk-definition permittivity F/m | Do not conflate ``\varepsilon_i`` with unresolved relative snow convention in (17) |
| ``j,\omega,t`` | unchanged | Imaginary unit, angular frequency rad/s, time s | Frequency derivative ``j\omega``; explicit full time exponential not printed |
| ``x,y,z,h_i,h_j,d`` | unchanged | Longitudinal/lateral/vertical coordinates, heights above snow, snow thickness, m | Figs. 1–3; upward ``z``; no earth-surface height renaming |
| ``d_{ij},D_{ij}`` | unchanged | Conductor separation and separation to image in air/snow plane, m | (8); no supplied general radius regularization |
| ``s,\lambda`` | unchanged | Final and preceding spectral integration variables | Final ``s`` in m⁻¹ by exponent; printed transformation conflicts retained |
| ``a_i,b_i`` | unchanged | Transverse spectral root (m⁻¹) and root/permeability coefficient | After (7); not conductor radius ``a`` |
| ``A_1,A_2,A_3,A_4,c_0,\ldots,c_6`` | unchanged | Algebraic interface coefficients; individual units not separately stated | Definitions after (7); ``A_3,A_4`` appear only in EN; original grouping unresolved |
| ``\tau_i^2`` (JP) / ``\tau_i`` (EN), interface context | unchanged, witness-labelled | Dimensionless ratio as printed | Both used with squared ``\tau_1,\tau_2`` in coefficients; left-hand-power difference retained |
| ``\varepsilon_s,\varepsilon_{s\infty},\varepsilon',\varepsilon''`` | unchanged | Snow complex/high-frequency/real/imaginary parameters | (17); absolute/relative units and EN extra ``\omega`` unresolved |
| ``\tau_1,\tau_2`` in (17) | unchanged, constitutive context | Relaxation times, Table 1 µs | Not silently equated with interface ratios |
| ``C_1,C_2`` in (17) | unchanged, constitutive context | Dielectric strengths, Table 1 pF as printed | Dimensional mapping and original parameter-determination procedure unresolved; not plate capacitances |
| ``C_1,C_2,C`` in (20)–(21) | unchanged, equivalent-capacitance context | Air/snow layer and total equivalent capacitances | Separate physical meaning; per-length/plate-area normalization not explicitly specified |
| ``C_i,r,h`` in EN ideal comparison | unchanged | Ideal self-capacitance pF/m, radius and no-snow height m | Unnumbered EN p. 29; not a new definition of preceding ``r=\sqrt{x^2+y^2}`` |
| ``C_s,C_0,\mathrm{DY},k`` | unchanged | Snow/no-snow capacitance per length, relative change percent, depth slope %/m | (19), §4.1; case-dependent observation, not universal accuracy measure |
| ``m,T,a,\rho_e`` | unchanged | Example snow density g/cm³, temperature °C, wire radius cm, earth resistivity Ωm | Table 1 and §4; ``m`` is density here, not a layer index |
| ``\mathbf A,\boldsymbol\Pi,\mathbf E,\mathbf H`` | unchanged | Vector potential, Hertz vector, electric field, magnetic field | (1)–(2); JP/EN electric-field sign differs; no new vector-potential normalization |

## Evidence and approximation sources

The Japanese original's title, byline, volume and pages agree with the published translation's explicit source header and the J-STAGE metadata. Its received/revised dates on p. 277 are March 31/August 2, 1999; they are not publication-year substitutions. The two versions are preserved as distinct mathematical witnesses, not as two independently original snow contributions.

The source chain is the Hertz dipole interface model (3)–(6), longitudinal integration with the stated imposed current, potential definition (9) and charge/current relations (11)–(13), leading to (14)–(15). The printed final kernels have all their named interface coefficients here, but original ``A_2`` grouping and the collisions explicitly prevent claiming a fully resolved evaluation recipe. No numerical quadrature, kernel merge, root choice or repair was performed. No derivation from modern LCM code was used as publication evidence.

The no-snow formula removes the finite insulating layer using all four source-stated substitutions, not just ``d=0``. The source explicitly identifies (16) with Nakagawa [6], *Admittance correction effects of a single overhead line* (1981). Wise [3] is likewise an attributed predecessor; the separate [Wise record](../../1948/homogeneous-earth-overhead-potential-coefficient/Wise1948.md) is original evidence for its own formula, not a repair of (16). Carson [2] is used for the **series impedance** after neglecting snow effects on that quantity; no new external-impedance formula follows from this use.

The dipole/interface construction cites Ametani, Nakagawa, and Iwamoto, *多層大地におけるサージ伝搬特性*, 95-B(10), pp. 500–506 (1975), reference [4]. The [official J-STAGE issue listing](https://www.jstage.jst.go.jp/browse/ieejpes1972/95/10/_contents/-char/ja) identifies DOI `10.1541/ieejpes1972.95.500`. Equations (1)–(6) restate three-, two-, and one-layer earth impedance. Sections 3–5 apply that impedance to propagation and switching surges; no snow potential or admittance kernel is printed there. The 1973 original and 1975 restatement are retained in the [three-layer impedance record](../../../external-impedance/1973/three-layer-earth-overhead-integral/Nakagawa1973.md). This methodology citation does not place the snow potential expression in the 1975 paper.

The approximate series-capacitance interpretation has two distinct operations: geometrical replacement by air/snow plate capacitances, then dominant-capacitance reduction under ``C_1\ll C_2``. Neither an expansion order nor a remainder is supplied. The approximate equality of relative admittance/capacitance changes in (19) uses small conductance relative to susceptance in the illustrated case. The printed near-linear depth rule is confined to the stated example. The reported roughly 3% low-frequency/1% high-frequency capacitance changes compare snow/no-snow models in Fig. 5; they are not measured transcription accuracy or integral-approximation errors.

The later [Papadopoulos–Papagiannis–Labridis record](../../2009/two-layer-earth-overhead-potential-correction/Papadopoulos2009.md) cites the 2001 snow version as a similar predecessor. The independently inspected original now establishes the snow witness itself, but this record does not prove complete algebraic equivalence, import the later self prescription, or transfer a later matrix assembly to the earlier publication.

## Limitations and discrepancies

1. **Unresolved published grouping:** JP ``A_2`` has a slash before one braced factor and a following dot/product line with no encompassing denominator bracket. Its literal layout is retained. EN explicitly has ``/(A_3A_4)``; this is a separate published witness, not an erratum resolving JP.
2. **Directly observed witness differences:** JP ``A_1`` has ``(s-\mu_0b_1)c_2``; EN has ``(s-\mu_0b_2)c_2``. JP defines ``\tau_i^2=\gamma_{i-1}^2/\gamma_i^2``; EN prints ``\tau_i=\gamma_{i-1}^2/\gamma_i^2`` while still using squared ``\tau`` in the coefficients. JP (2) has ``-\gamma^2\boldsymbol\Pi`` in the electric field; EN has a plus. JP (17) ends ``\varepsilon'-j\varepsilon''``; EN ends ``\varepsilon'-j\omega\varepsilon''``. All were verified in page images; none is classified as an OCR/conversion defect.
3. **Unresolved transformation token:** the statement before (7) uses unindexed ``\gamma`` in JP and ``r`` in EN, the latter already a length. This is a suspected published derivation/notation defect, not a license to substitute ``\gamma_0`` or reconstruct a transform. No isolated-term singularity claim is made.
4. **Unresolved material dependencies:** interface ratios and microsecond relaxation times reuse ``\tau_1,\tau_2``; dielectric strengths and equivalent plate capacitances reuse ``C_1,C_2``. Snow permittivity is discussed as relative, yet (17) includes a conductivity term and Table 1 labels strengths pF. The source does not explicitly reconcile these units with bulk ``\varepsilon_1`` or specify how the loss contribution is partitioned from separate ``\sigma_1``. No missing ``\varepsilon_0`` or subtraction of conductivity is introduced.
5. **Missing self and branch evidence:** a circular-wire numerical example and an EN ideal self formula do not specify the general finite-radius self limit of the snow integral. No square-root branch is printed. These dependencies remain open rather than imported from later authors or implementation convention.
6. **Output/model restriction:** (14) is a potential coefficient. The explicitly printed matrix operation (18) is adjacent to the perfect-earth/no-snow coefficient. Neither author-labelled admittance prose nor the family name authorizes an elementwise inverse, invented snow self matrix, buried/mixed extension, arbitrary-layer recurrence or new internal impedance. The equivalent plate model is explanatory, not the exact spectral result.
7. **Citation and historical limits:** the 2000 original and 2001 English translation remain separate bibliography entries for one formula. The 1975 paper supplies an earlier impedance restatement, not this snow potential kernel. The method reference does not establish complete algebraic equivalence.
