# Iwamoto logarithmic-integral evaluation of overhead earth impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Thin overhead line-current model with self radius ``a\ll h``; conductor skin impedance and insulation are excluded. |
| Calculated quantities | Normalized earth-return resistance and reactance evaluated by a logarithmic-coordinate integral and Duhamel operator. |
| Earth structure | Uniform scalar resistivity ``\rho_e`` in a half-space below a plane surface. No finite layer or arbitrary stratification is included. |
| Model and approximation | The appendix transforms and numerically evaluates Carson's kernel without changing its material model. It gives no uniform error estimate, and its printed inconsistencies prevent algebraic verification of every equality. |
| Main source | K. Iwamoto (January 1958), evaluating the Carson kernel. |
| Citation key(s) | Primary: `:Iwamoto1958a`; physical-kernel source: `:Carson1926` |
| Evidence status | Original pages checked. The exponential factor, radius mapping, operator arguments, and resistance normalization remain unresolved. |

**Description.** Logarithmic-coordinate integral and Duhamel evaluation representation for normalized self and mutual earth-return impedance of parallel overhead wires above homogeneous conducting earth.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | No independent imposed longitudinal parameter is printed in the earth kernel. The paper computes a line propagation constant from distributed parameters for its wave-response application; that solved quantity is not an imposed-current prescription. Do not assign author-stated zero from its absence in (2). | Not stated — §2 (2)–(3). Stated — §3 (7)–(10), p. 24, distinguishes the solved line quantity. |
| Air propagation constant ``γ_air`` | Not independently defined in the extracted kernel. Its geometric term uses the perfect-earth external-field reference, not an explicit retained air bulk wavenumber. | Not stated — §2. Stated — decomposition paragraph preceding (4), p. 24. |
| Earth propagation constant ``γ_earth`` | No bulk gamma definition is supplied for this representation. Conduction is encoded in ``r_e`` and ``F(\lambda)=\sqrt{\lambda^2+j}-\lambda``. The square-root branch is not specified. | Stated — definitions following (2)–(3) and appendix, pp. 23–24,31. Unresolved — explicit branch. |
| Earth permittivity and displacement current | Neglected in the adopted Carson model. The introduction discusses when omitted permittivity might matter, but the extracted expression does not reinstate it. | Stated — §2 first paragraphs, p. 23, homogeneous-earth/displacement-current-omission description; introductory discussion on the same page. |
| Range of validity | Homogeneous-earth, thin overhead wire model; self radius ``a\ll h`` is explicit. The charted ``\theta`` values and wave-response test cases are tested ranges, not a universal bound. No convergence/error guarantee for every logarithmic-integral parameter is printed. | Stated — prose before (3), p. 23; appendix sampling list, p. 31. |
| Earth permeability ``μ_earth`` | No independent earth-permeability input is present. The source's MKS coefficient ``10^{-7}`` is fixed and must not be replaced by an arbitrary permeability parameter. | Equation-implied — (2)–(3) and ``r_e`` definitions. Not stated — separately named earth permeability. |
| Arrangement | Parallel overhead source and target at heights ``h,h'`` and horizontal separation ``D``. Mutual output ``\dot Z'``; self by ``D=0,h'=h``, with source-prescribed finite radius in (3). | Stated — Fig. 3 and §2, p. 23. |
| Earth structure | Uniform scalar resistivity ``\rho_e`` in a half-space below a plane surface. No finite layer or arbitrary stratification is included. | Stated — §2 and Fig. 3, p. 23. |
| Conductor and insulation geometry | Infinitesimally thin external line-current construction; self radius ``a\ll h`` used in the logarithmic geometric term. The conductor's separately discussed Bessel skin impedance is not part of the transformed earth output. No insulation formula is supplied. | Stated — §2 geometry and self prescription, p. 23; separately introduced conductor term (5), p. 24. |
| Constitutive and field assumptions | Adopted Carson homogeneous-earth calculation with neglected displacement current; scalar resistivity and linear harmonic response. No new boundary-value solution or nonlinear/anisotropic soil model is supplied. | Stated — §2 attribution and model discussion, p. 23. Equation-implied — (2)–(4). |
| Conventions | ``z`` longitudinal, ``y`` vertical, plane earth at ``y=0`` in Fig. 3. The source evaluates at ``p=j\omega``; no separate full time exponential is printed beside the kernel. MKS external impedance is per length; ``R_e,R'_e`` in the main text are geometry-dependent Ω/m scales. Chart ``R_e`` is separately labelled Ω/km and must not be silently identified with the main-text numerical value. | Stated — Fig. 3, (2)–(6), and appendix chart. Equation-implied — MKS dimensions and resistance definitions. |

**Expression.** The appendix's **unnumbered** logarithmic-integral representation, p. 31/PDF page 10. The actual output is normalized earth impedance ``P_e+jQ_e``, not the complete external impedance:

```math
P_e+jQ_e=
\frac{r_e}{4}\int_{-\infty}^{\infty}
F_l(u-\zeta)e^{-\cos\theta\,e^\zeta}
\cos(\sin\theta\cdot e^\zeta)e^\zeta\,d\zeta
=\frac{r_e}{4}\{F_l(u)\star E_{lc\theta}(u)\}.
```

Its published helpers and substitution prose:

```math
F(\lambda)=\sqrt{\lambda^2+j}-\lambda,\qquad
F_l(u)=F(e^{-u}),
```

```math
r=e^u,\qquad \lambda=e^{-\zeta},
```

```math
E_{lc\theta}(u)=
\cos\theta-e^{-e^u\cos\theta}\cos(e^u\sin\theta+\theta).
```

**Unresolved mapping:** the substitution prose prints unindexed ``r``, while the surrounding equations use ``r_e``. The record preserves that difference rather than inserting an unprinted equality. The final helper's angular shift is explicitly **plus** ``\theta``. Its compatibility with the printed operator representation is not silently repaired.

The geometry and normalized scale supplied with the Carson-attributed mutual parent (2), p. 23:

```math
\theta=\tan^{-1}\frac{D}{h+h'},\qquad
r_e=4\sqrt{\frac{\omega}{R'_e}\times10^{-7}},\qquad
R'_e=\frac{4\rho_e}{\pi\{D^2+(h+h')^2\}}.
```

For the expressly stated self case ``D=0,h'=h``, (3), pp. 23–24:

```math
r_e=4\sqrt{\frac{\omega}{R_e}\times10^{-7}},\qquad
R_e=\frac{\rho_e}{\pi h^2}.
```

The printed impedance assembly relation, p. 24, is retained **with its unprimed right-hand scale**, including when the left lists mutual and self quantities:

```math
\dot Z'_e,\ \dot Z_e
=R_e\{P_e(r_e)+jQ_e(r_e)\}.\qquad\text{(4)}
```

The mutual normalization's relation to ``R'_e`` above is unresolved in this compressed printed line; no prime is inserted. The distinct parent (2) and self (3) are given in Source transcription without altering them.

The crossed 相乗 operator is displayed as ``\star``, a notation-only alias. Section 4 defines it on p. 25:

```math
\varphi(u)=\varphi_0(u)H(u-\infty),\qquad
\psi(u)=\psi_0(u)H(u-\infty),
```

```math
\varphi(u)\star\psi(u)=\varphi_0(-\infty)\psi_0(\infty)
+\int_{-\infty}^{\infty}
\varphi'_0(u-\zeta)\psi_0(\zeta)\,d\zeta.\qquad\text{(18)}
```

Here ``H`` is the unit step (description after (12), p. 24); ``\varphi,\psi`` are generic operands, the subscript-0 functions are their printed helpers, and the prime is differentiation. The prose says the discontinuity is at **negative** infinity, but both printed step arguments are ``u-\infty``. Preserve the source conflict and endpoint term. The source says that the endpoint term vanishes if one or both operands do not have the stated discontinuity, and that the subscript 0 can then be removed; this does not authorize discarding it unconditionally.

**Approximation.** The appendix presents a change of integral representation and a graphical/numerical evaluation of the Carson kernel, not a new material model or an analytical truncation of that kernel. This is more than renaming an integration variable: the explicit logarithmic integration, new helper and Duhamel assembly are retained together as the published evaluator. The associated waveform approximations in §§3–6 have different outputs and are not counted as independent earth-impedance formulations. No uniform error estimate is supplied for the appendix integral. Its literal internal inconsistencies prevent a claim that all printed equalities are mathematically verified.

**Limitations.** The appendix's opening integral omits ``\lambda`` from its exponential, whereas the main parent and chart include it. The helper prints ``+\theta``, and its derivative/operator relationship remains unresolved. Unindexed ``r``, resistance-scale/prime usage, the infinity-sign conflict and root branch remain visible. No corrected evaluator is provided.

**Reference.** [Iwamoto1958a](@cite).  Iwamoto, January 1958, appendix pp. 31–32, with §2 (2)–(4), pp. 23–24, and §4 (18), p. 25. [Carson1926](@cite) identifies the earlier physical kernel expressly credited by Iwamoto; the original c.g.s. transcription is in the separate Carson record.

**Transcription source.** Original January 1958 publication images. Magnified views verify the main exponential, appendix's missing variable, helper's plus sign, logarithmic domain and measure, normalization, and operator endpoints. Main, appendix, and chart variants remain separate. The earlier physical kernel remains attributed to Carson.

## Source transcription

The main parents precede the transformed representation. Source labels and quantities are retained; the appendix integrals are unnumbered and receive prose locators only.

### §2, mutual parent, p. 23/PDF page 2

```math
\dot Z'=\omega\left\{
j\,2\log\sqrt{\frac{D^2+(h+h')^2}{D^2+(h-h')^2}}
+4\int_0^\infty(\sqrt{\lambda^2+j}-\lambda)
e^{-r_e\cos\theta\cdot\lambda}
\cos(r_e\sin\theta\cdot\lambda)\,d\lambda
\right\}\times10^{-7}\quad\text{(MKS)}.\qquad\text{(2)}
```

```math
\theta=\tan^{-1}\frac{D}{h+h'},\qquad
r_e=4\sqrt{\frac{\omega}{R'_e}\times10^{-7}},\qquad
R'_e=\frac{4\rho_e}{\pi\{D^2+(h+h')^2\}}.
```

The logarithm remains outside the separate correction integral. The source's fixed MKS factor is copied, not derived by converting Carson's original c.g.s. units. ``\dot Z'`` is mutual **external** impedance; the conductor-internal term introduced in §3 is separate.

### §2, self parent, pp. 23–24/PDF pages 2–3

With ``D=0,h'=h`` and ``a\ll h`` as explicitly stated:

```math
\dot Z=\omega\left\{
j\,2\log\frac{2h}{a}
+4\int_0^\infty(\sqrt{\lambda^2+j}-\lambda)
e^{-r_e\lambda}\,d\lambda
\right\}\times10^{-7}\quad\text{(MKS)}.\qquad\text{(3)}
```

```math
r_e=4\sqrt{\frac{\omega}{R_e}\times10^{-7}},\qquad
R_e=\frac{\rho_e}{\pi h^2}.
```

```math
\dot Z'_e,\ \dot Z_e=R_e\{P_e(r_e)+jQ_e(r_e)\}.\qquad\text{(4)}
```

The source describes ``R'_e,R_e`` as unit-length resistances of conceptual circular earth conductors with radii ``\sqrt{D^2+(h+h')^2}/2`` and ``h`` respectively. They are reference resistance scales, not the real parts of the complex earth correction. The printed common unprimed scale in (4) is not silently replaced for the mutual case.

### §4, operator, p. 25/PDF page 4

```math
\varphi(u)=\varphi_0(u)H(u-\infty),\qquad
\psi(u)=\psi_0(u)H(u-\infty),
```

```math
\varphi(u)\star\psi(u)=\varphi_0(-\infty)\psi_0(\infty)
+\int_{-\infty}^{\infty}
\varphi'_0(u-\zeta)\psi_0(\zeta)\,d\zeta.\qquad\text{(18)}
```

The source calls this the Duhamel operation, extended to the stated negative-infinity discontinuity. The step-sign discrepancy and endpoint convention remain as printed.

### Appendix, p. 31/PDF page 10, opening integral

**Literal appendix witness — no ``\lambda`` appears in this exponential:**

```math
P_e+jQ_e=\frac{r_e^2}{4}\int_0^\infty
F(\lambda)e^{-r_e\cos\theta}
\cos(r_e\sin\theta\cdot\lambda)\,d\lambda,\qquad
F(\lambda)=\sqrt{\lambda^2+j}-\lambda.
```

This is a demonstrated difference from main (2) and the p. 32 chart, not a conversion error. The record does not insert the missing variable or claim convergence/failure from inspecting one term alone.

### Appendix, p. 31, substitution and final representation

The prose prints unindexed ``r=e^u`` and ``\lambda=e^{-\zeta}`` and describes ``F(e^{-u})`` by ``F_l(u)``:

```math
r=e^u,\qquad \lambda=e^{-\zeta},\qquad F_l(u)=F(e^{-u}),
```

```math
E_{lc\theta}(u)=\cos\theta-e^{-e^u\cos\theta}
\cos(e^u\sin\theta+\theta),
```

```math
P_e+jQ_e=\frac{r_e}{4}\int_{-\infty}^{\infty}
F_l(u-\zeta)e^{-\cos\theta\,e^\zeta}
\cos(\sin\theta\cdot e^\zeta)e^\zeta\,d\zeta
=\frac{r_e}{4}\{F_l(u)\star E_{lc\theta}(u)\}.
```

The final measure is ``e^\zeta d\zeta`` as printed, not ``d\zeta`` alone. The shift is ``u-\zeta``. The helper's plus angle is literal and is not changed using an inferred antiderivative.

### Appendix, p. 32/PDF page 11, independent chart witness

Fig. 付第1 repeats the normalized parent **with** ``\lambda`` in the exponential:

```math
P_e+jQ_e=\frac{r_e^2}{4}\int_0^\infty
(\sqrt{\lambda^2+j}-\lambda)e^{-r_e\cos\theta\cdot\lambda}
\cos(r_e\sin\theta\cdot\lambda)\,d\lambda.
```

The same inset prints:

```math
r_e=4\sqrt{\frac{\omega}{R_e}\times10^{-7}}
=0.10026\sqrt{f/R_e}.
```

Its geometry box labels the earth resistance in Ω/km with a factor of 1000. The same unprimed ``R_e`` is reused in the main text for Ω/m; under one common numerical resistance scale the two displayed frequency prefactors do not agree. This is an unresolved published unit/notation switch, not a silently inserted prime or coefficient repair. The August source repeats a related issue; its two-layer physical kernel is a separate record.

The appendix's calculation list on p. 31 prints ``\theta=0,\pi/8,\pi/4,3\pi/8,7\pi/16,15\pi/32,0``, ending with a repeated zero. The p. 32 charts also label ``\theta=\pi/2``. Both are retained as different printed sampling evidence; the final listed zero is not silently changed to a right angle.

Fig. 付第2's inset separately prints the **primed** earth ratio:

```math
P_e+jQ_e=\dot Z'_e/R'_e.
```

That primed ratio remains distinct from the compressed unprimed scale in main (4); it is not used to rewrite (4). The shared charts also show the known cylindrical-conductor function, which is not part of this earth-only evaluator.

### Related asymptotic parent used for waveform approximation

The unnumbered earth series immediately before (11), p. 24/PDF page 3, is:

```math
P_e(r_e)+jQ_e(r_e)\simeq
\frac14\left(\sqrt j\,r_e-1+
\frac{1}{\sqrt j\,r_e}-\frac{3}{j r_e^2}+\cdots\right).
```

The printed inverse-square term and ellipsis are preserved. The source then uses two terms in its propagation/waveform approximation; no uniform earth-impedance error bound is supplied. This is recorded as a parent-series witness, not a separately claimed new approximation by Iwamoto. It was compared with the inventoried Carson source at the level needed to establish parentage; an expected different order/coefficient is not inserted.

## Notation map

The crossed 相乗 glyph is displayed one-to-one as ``\star``; no other renaming is performed. In particular, ``r`` is not silently renamed ``r_e`` and the prime on ``R'_e`` is not dropped or inserted to fix a witness.

| Source → display | Meaning / units | Convention and locator |
| --- | --- | --- |
| ``\dot Z',\dot Z`` → same | Mutual/self external impedance, Ω/m in main MKS expressions | Dot is complex notation; prime distinguishes mutual, not differentiation |
| ``\dot Z'_e,\dot Z_e`` → same | Mutual/self earth-return correction | Does not include the logarithmic geometric term |
| ``P_e,Q_e`` → same | Dimensionless normalized resistance/reactance | Not potential or admittance coefficients |
| ``\rho_e`` → same | Homogeneous earth resistivity, Ωm | Uniform scalar input |
| ``D,h,h',a`` → same | Lateral separation, heights and self wire radius, m | Self ``a\ll h``; prime on height is target identification |
| ``x,y,z`` → same | Cartesian coordinates | Fig. 3: vertical y, longitudinal z |
| ``R'_e,R_e`` → same | Mutual/self conceptual-earth-conductor resistance scales, Ω/m | (2)–(4); common unprimed scale in (4) unresolved |
| Chart ``R_e`` and ratio ``R'_e`` → same | Earth resistance labelled Ω/km, with separately primed ratio notation | Distinct numerical scale from main text; not harmonized |
| ``\theta`` → same | Geometry angle, radians | ``\tan^{-1}(D/(h+h'))``; sampled values are not a universal domain theorem |
| ``\omega,f,j`` → same | Angular frequency, cyclic frequency, imaginary unit | Positive j retained; fixed MKS coefficient |
| ``r_e`` → same | Dimensionless frequency/geometry parameter | Formula-specific resistance denominator |
| Unindexed ``r`` → same | Variable in appendix logarithmic substitution | Exact identification with ``r_e`` is not printed |
| ``\lambda`` → same | Dimensionless spectral integration variable | Positive real half-line in parent |
| ``u,\zeta`` → same | Logarithmic coordinate and integration variable | Full real-line integral, source substitutions preserved |
| ``F,F_l,E_{lc\theta}`` → same | Dimensionless spectral and transformed functions | Exact nested exponentials and plus-angle helper retained |
| Crossed 相乗 glyph → ``\star`` | Defined Duhamel operation | Endpoint and derivative terms preserved, (18) |
| ``\varphi,\psi,\varphi_0,\psi_0,H`` → same | Generic operands, helper functions, unit step | Printed infinity-sign conflict remains visible |
| Prime on ``\varphi'_0`` → same | Derivative of generic helper | Different from mutual and target-height primes |

## Evidence and approximation sources

1. **Distinct evaluator, inherited physics.** Section 2 expressly adopts Carson and contrasts it with earlier approximate wave treatments. The appendix describes evaluating the Carson integral by the author's ``E_{lc\theta}`` graphical transformation. The record is for that published integral/evaluation representation, not an additional count of Carson's physical kernel for a notation change. The [Carson original record](../../1926/homogeneous-earth-overhead-integral/Carson1926.md) remains independent and in original c.g.s. units.
2. **No implicit repaired equality.** The main kernel and chart supply a decaying spectral exponential; the appendix opening line lacks its variable. The transformed integral, printed helper and defined Duhamel operation are all independently preserved, even though their consistency is unresolved. No unprinted antiderivative, sign change, root branch or normalized resistance replacement is supplied.
3. **Approximation scope.** Section 3 expands the adopted earth and standard solid-cylinder impedances and keeps two terms to approximate the solved propagation constant and wave response. Sections 5–6 compare finite-conductivity and zero-conductor-resistance waveforms. Their error curves and graph-integration estimates concern those applications, not a new general earth model or a global bound for this appendix evaluator.
4. **Self and mutual are not inferred extensions.** Main (2) explicitly supplies mutual height/separation geometry, and (3) explicitly supplies self radius with ``a\ll h``. These are earlier-attributed parents retained with the evaluator. No arbitrary-radius or proximity correction is added to either.
5. **Other expressions.** The paper describes its solid-cylinder Bessel expression as known, so it is not a new internal-impedance contribution. Its scalar lossless space capacitance is an input to a wave model, not a new earth-admittance expression. Kikuchi's 1957 high-frequency potential formulation is documented separately.
6. **Earlier transformation history.** January §4 cites Iwamoto's 1954 JIEE 74 pp. 293 and 284 for general functional transformations. The present source itself defines the operator needed here. Those generic-method references do not establish an unseen earlier publication of this specific earth-impedance evaluator; no priority claim is transferred to them.
7. **Later use remains separate.** The [August two-layer witness](../two-layer-earth-overhead-self-integral/Iwamoto1958b.md) applies related machinery to a finite-depth layer and explicitly credits Sunde's kernel. It is not a correction of the January helper, and its simpler ``E_l`` must not replace the January angle-dependent function.

## Limitations and discrepancies

1. **Demonstrated printed witness difference:** appendix p. 31 omits ``\lambda`` in the exponential that is present in main (2) and the p. 32 chart. The source scan is legible; this is not an OCR-only defect.
2. **Unresolved evaluator relationship:** the appendix prints ``r=e^u`` without the subscript of surrounding ``r_e``, and the angle-dependent helper contains ``+\theta``. Their full relationship to the final integral and operator definition is not established by merely claiming a change of variables. No repair is selected.
3. **Demonstrated operator wording conflict:** negative-infinity discontinuity in prose versus literal ``H(u-\infty)`` in (18)'s helpers. The endpoint factor ``\psi_0(\infty)`` is also preserved exactly, not changed to a negative-infinity value.
4. **Normalization:** the main mutual scale ``R'_e``, common unprimed assembly scale (4), and chart Ω/km scale are kept separate. The chart's decimal prefactor does not follow from its first equality under one common numerical resistance convention. Evaluate the main integral using the Carson parent in SI units.
5. **Sampling-list discrepancy:** the appendix's angle list ends in a repeated zero, while the chart labels include a right angle. The list is retained literally and is not presented as a universal validity range.
6. **Parent series not repaired:** the earth expansion's printed ``-3/(jr_e^2)`` term remains as printed; original-series comparison is pending. The ellipsis is not completed from memory.
7. **Evidence limit:** source-page transcription does not establish physical validity, convergence for all parameters, consistency of the printed equalities, or earliest historical priority. Missing definitions are not replaced with software conventions.
