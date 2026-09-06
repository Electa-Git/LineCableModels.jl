# Papadopoulos–Papagiannis–Labridis two-layer overhead earth potential correction

## Identity and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Two infinite, electrically thin overhead conductors; the self term uses the source-prescribed radius substitution. Conductor internal and insulation terms are excluded. |
| Calculated quantities | Mutual earth-return potential correction, self correction, homogeneous reduction, and admittance-matrix assembly. |
| Earth structure | Air above a finite earth layer of thickness ``d`` and a lower earth half-space. |
| Model and approximation | Integral representation under the thin-conductor quasi-TEM model. The 2010 derivation replaces the unknown longitudinal constant by the air value and applies a Bessel transform; no additional matrix inversion is applied to (2b). |
| Main source | T. A. Papadopoulos, G. K. Papagiannis, and D. A. Labridis (2009); the 2010 paper supplies the auxiliary derivation. |
| Citation key(s) | Primary: `:Papadopoulos2009`; auxiliary derivation: `:Papadopoulos2010a` |
| Evidence status | Both publications checked against page images. Prime/index conflicts, the root branch, and scalar-inverse interpretation remain unresolved. |

**Description.** Per-unit-length earth-related potential-coefficient correction for self and mutual overhead thin-conductor interactions above two horizontal earth layers, including a radial-displacement-current kernel in addition to the impedance kernel. The source prints an associated admittance relation; a later witness explicitly supplies full potential-matrix assembly.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | 2009 states quasi-TEM but gives no explicit imposed longitudinal constant or exponential; ``a_m=\sqrt{u^2+\gamma_m^2-\gamma_0^2}`` is retained without turning it into an author-stated ``Γ`` prescription. 2010 parent (4b) uses unknown ``\gamma_x`` and ``e^{-\gamma_xx}``; final (6) uses the stated ``\gamma_x=\gamma_0=jk_0=j\omega\sqrt{\varepsilon_0\mu_0}``, not lossy earth propagation. | Stated — 2009 I–II.A, appendix; 2010 p. 1162 between (4) and (5). Unresolved — explicit original prescription. |
| Air propagation constant ``γ_air`` | 2009 air has ``\mu_0,\varepsilon_0``; indexed bulk definition gives ``\gamma_0^2=j\omega\mu_0(\sigma_0+j\omega\varepsilon_0)`` without separately specifying ``\sigma_0`` there. 2010 explicitly sets the longitudinal prescription to free-space ``\gamma_0=j\omega\sqrt{\varepsilon_0\mu_0}``. Air terms remain in ``G`` and the interface coefficients. | Stated/equation-implied — 2009 II.A, (2c), appendix; 2010 (6d) and longitudinal prescription. |
| Earth propagation constant ``γ_earth`` | ``\gamma_m^2=j\omega\mu_m(\sigma_m+j\omega\varepsilon_m)``, with earth indices 1 and 2; conductivity and permittivity retained. Final 2009 ``a_m`` and 2010 ``\alpha_k`` definitions are preserved separately; no explicit root branch. | Stated — 2009 appendix p. 1067; 2010 appendix p. 1169. |
| Earth permittivity and displacement current | ``\varepsilon_1,\varepsilon_2`` arbitrary scalar layer values; the authors attribute ``G`` to radial displacement currents and retain axial displacement effects in both media. No independent loss tangent, complex permittivity, or loss-conductivity separation is supplied. The claim about lossless propagation when omitting ``G`` is retained only as an author statement. | Stated/equation-implied — 2009 II.A,III, (2c); 2010 section 3.1 and (6d). |
| Range of validity | Electrically thin conductors and quasi-TEM TL-mode approximation, not all full-wave modes. 2010 states ``f\ll c/h`` as the TL limit and discusses conditional use beyond it. Its 50 Hz–10 MHz study is not a universal validity or accuracy bound. No kernel-series truncation order is specified. | Stated — 2009 I–II.A; 2010 pp. 1161–1163,1168. |
| Earth permeability ``μ_earth`` | Arbitrary scalar ``\mu_1,\mu_2`` in the actual kernel; unit relative permeability is only the numerical-example choice. | Stated — 2009 II.B,IV; equation-implied — (2c), (A.1)–(A.4); 2010 sections 3.3,4.2. |
| Arrangement | Both conductors above the upper air/earth interface. Mutual ``i,j`` and self by replacing ``y_{ij}`` with outer conductor radius and ``h_j`` with ``h_i``. No underground or mixed arrangement is supplied. | Stated — 2009 Fig. 1, II.A; 2010 Fig. 2/prose after (6). |
| Earth structure | Air half-space; finite earth layer of thickness ``d``; infinite lower earth half-space. Homogeneous earth obtained by identical-layer prescription. No explicit arbitrary-layer recursion. | Stated — 2009 II.A,II.C; 2010 section 3.1. |
| Conductor and insulation geometry | Uniform electrically thin perfect parallel conductors; original 2009 uses infinite lengths. Later derivation integrates infinite source ``i`` against per-length target ``j``. The self radius is geometric; no conductor insulation or dielectric-layer admittance is part of this earth correction. | Stated — 2009 II.A; 2010 section 3.1. |
| Constitutive and field assumptions | Scalar uniform layer properties and linear Hertz-vector superposition; quasi-TEM final model. Potential correction is not a shunt conductance alone and not an insulation capacitance. No additional bonding constraint or internal skin/proximity model is inferred. | Equation-implied — (1)–(2) and appendix; stated — 2009 I–II.A, 2010 sections 2–3.1. |
| Conventions | ``j`` imaginary, ``\omega=2\pi f``; no explicitly printed time exponential. Heights upward from air/earth surface. Later appendix uses upward ``z`` with lower interface at 0, upper interface at ``d``. ``P_g`` is a potential-coefficient correction; source-labelled ``Y'_g`` is per-length admittance with an inverse. Only the 2010 section 2 explicitly calls ``\mathbf P`` an ``N\times N`` matrix. | Stated — 2009 p. 1065/Fig. 1; 2010 (3), Fig. 18, appendix layer ranges. |

**Expression.** The 2009 source prints the following **admittance correction relation and potential correction**, (2a)–(2c), p. 1065:

```math
Y'_{g_{ij}}=j\omega P_{g_{ij}}^{-1}.\qquad\text{(2a)}
```

```math
P_{g_{ij}}=\frac{1}{\pi\varepsilon_0}
\int_0^\infty [F(u)+G(u)]e^{-u(h_i+h_j)}\cos(y_{ij}u)\,du.\qquad\text{(2b)}
```

```math
G(u)=u\frac{
\mu_0\mu_1(\gamma_0^2-\gamma_1^2)
(s_{12}+d_{12}e^{-2\alpha_1d})(S_{12}+D_{12}e^{-2\alpha_1d})
-4\mu_0\mu_1^2\mu_2\alpha_1^2\gamma_0^2
(\gamma_2^2-\gamma_1^2)e^{-2\alpha_1d}
}{\Delta_2\cdot\Delta}.\qquad\text{(2c)}
```

The leading ``u`` multiplies the **whole** fraction. Both parenthesized interface factors multiply one another. The complete shared impedance kernel (1b) is repeated as a dependency, without folding its height exponential into it:

```math
F(u)=\mu_1\frac{s_{12}+d_{12}e^{-2\alpha_1d}}
{s_{01}s_{12}+d_{01}d_{12}e^{-2\alpha_1d}}.\qquad\text{(1b)}
```

All six original appendix definitions, p. 1067, retain their printed Latin/prime notation:

```math
s_{mn}=(a_m\mu_n+a_n\mu_m),\qquad\text{(A.1)}
```

```math
d_{mn}=(a_m\mu_n-a_n\mu_m),\qquad\text{(A.2)}
```

```math
S_{mn}=(\mu_m\gamma_n^2a_m+\mu_n\gamma_m^2a_n),\qquad\text{(A.3)}
```

```math
D_{mn}=(\mu_m\gamma_n^2a_m-\mu_n\gamma_m^2a_n),\qquad\text{(A.4)}
```

```math
\Delta=s_{01}s_{12}+d_{01}d_{12}e^{-2a'_1d},\qquad\text{(A.5)}
```

```math
\Delta_2=S_{01}S_{12}+D_{01}D_{12}e^{-2a_1d}.\qquad\text{(A.6)}
```

```math
a_m=\sqrt{u^2+\gamma_m^2-\gamma_0^2},\qquad
\gamma_m^2=j\omega\mu_m(\sigma_m+j\omega\varepsilon_m),
\qquad m,n=0,1,2,\qquad \omega=2\pi f.
```

The source does **not** explicitly identify Greek ``\alpha_1`` in the main kernel or primed ``a'_1`` in (A.5) with its defined Latin ``a_1``. These are specific unresolved dependencies, not omitted definitions to be reconstructed. ``F,G`` have length units inferred from the interface factors, so ``P_g`` has ``\mathrm{m/F}`` and the printed ``Y'_g`` has ``\mathrm{S/m}``. Layer indices are 0 air, 1 finite earth and 2 lower earth.

Self correction: substitute the conductor's outer radius for ``y_{ij}`` and ``h_i`` for ``h_j`` as stated below (2). Homogeneous earth: section II.C prescribes replacing ``a_2`` by ``a_1`` and ``\gamma_2`` by ``\gamma_1``. This is a source-prescribed identical-earth-layer reduction, not a separately derived closed form. The scalar/subscripted inverse in (2a) is retained; the record does not independently choose entrywise inversion or a matrix convention for this equation.

**Approximation.** Integral representation within the prescribed thin-conductor quasi-TEM model; no stated analytical-series order or truncation. The later 2010 parent equation (4b) contains initially unknown longitudinal propagation. Its replacement by the free-space value, retaining the TL mode, is the explicitly stated approximation; appendix (A.12) performs the longitudinal Bessel transform. Neither the source's claim that ``G`` may be omitted under another treatment nor a matrix inversion invented by the reviewer is applied to (2b).

**Limitations.** Potential correction, associated printed admittance correction, and full external admittance are different quantities. Original undefined ``\alpha_1,a'_1``, later mixed-prime denominator notation, and root branches remain unresolved. The 2010 explicit matrix relation is reported only as a later source-provided assembly, not as a license to reinterpret every scalar inverse or add admittances in parallel.

**Reference.** [Papadopoulos2009](@cite), (1b), (2a)–(2c), p. 1065; (A.1)–(A.6) and definitions, p. 1067. Later witness [Papadopoulos2010a](@cite), (3a)–(3c), (4b), (5c), (6a)–(6d), pp. 1161–1162; appendix A, pp. 1168–1170.

**Transcription source.** Original 2009 PDF pages 1–4 inspected, equations image-verified on pages 2 and 4. The 2010 equations below are independently image-verified against that publication. No converted text or LCM executable supplies a missing coefficient or repairs a sign/prime.

## Source transcription

The formula section preserves the original 2009 (2a)–(2c), its (1b) dependency and original appendix sequence. The following is the **distinct 2010 witness**, retained in its own main-text/appendix order.

2010 (3a)–(3c), p. 1161, explicitly uses the potential-coefficient **matrix**:

```math
\mathbf Y'(\omega)=\mathbf Y'_e(\omega)
=j\omega\mathbf P_e^{-1}=j\omega(\mathbf P_{pg}+\mathbf P_g)^{-1},\qquad\text{(2010:3a)}
```

```math
\mathbf Y'_{pg}=j\omega\mathbf P_{pg}^{-1},\qquad\text{(2010:3b)}
```

```math
\mathbf Y'_g=j\omega\mathbf P_g^{-1}.\qquad\text{(2010:3c)}
```

``\mathbf P`` is expressly an ``N\times N`` potential-coefficient matrix; ``pg`` is perfectly conducting ground and ``g`` the imperfect-earth correction. The paper's one-conductor equivalent circuit places the corresponding branches in series. It does **not** print ``\mathbf Y'_e=\mathbf Y'_{pg}+\mathbf Y'_g``. No such addition is supplied here.

2010 (4b), p. 1162, before the prescribed longitudinal approximation:

```math
Y_{e_{ij}}^{\prime-1}(\gamma_x)
=\frac{j\omega\mu_0}{4\pi\gamma_0^2}
\int_0^\infty
\left[\frac{u}{a_0}
\left(e^{-a_0|h_j-h_i|}+e^{-a_0(h_i+h_j)}\cdot T'_1\right)
+2a_0u\cdot T'_2\cdot e^{-a_0(h_i+h_j)}\right]
\times\left[\int_{-\infty}^{\infty}
J_0\!\left(u\sqrt{x^2+y_{ij}^2}\right)e^{-\gamma_xx}\,dx\right]du.
\qquad\text{(2010:4b)}
```

The leading ``j\omega\mu_0/(4\pi\gamma_0^2)`` and the inverse on ``Y'`` are exactly retained; neither is turned into a different output. ``J_0`` is the first-kind order-zero Bessel function. The nested integrals remain separate. Parent ``a_0`` versus appendix ``a'_0`` is a printed mismatch; no arbitrary-complex-``\gamma_x`` convergence prescription is supplied.

2010 (6a)–(6d), p. 1162, after the stated ``\gamma_x=\gamma_0=jk_0=j\omega\sqrt{\varepsilon_0\mu_0}``:

```math
Y'_{e_{ij}}=j\omega P_{e_{ij}}^{-1},\qquad\text{(2010:6a)}
```

```math
P_{e_{ij}}=P_{pg_{ij}}+P_{g_{ij}}
=\frac{1}{2\pi\varepsilon_0}\ln\frac{D_{ij}}{d_{ij}}
+\frac{1}{\pi\varepsilon_0}(M+jN),\qquad\text{(2010:6b)}
```

```math
M+jN=\int_0^\infty[F_{strat}(\lambda)+G_{strat}(\lambda)]
e^{-\lambda(h_i+h_j)}\cos(y_{ij}\lambda)\,d\lambda,\qquad\text{(2010:6c)}
```

```math
G_{strat}(\lambda)=\lambda\frac{
\mu_0\mu_1(\gamma_0^2-\gamma_1^2)
(s_{12}+d_{12}e^{-2\alpha_1d})(S_{12}+D_{12}e^{-2\alpha_1d})
-4\mu_0\mu_1^2\mu_2\alpha_1^2\gamma_0^2
(\gamma_2^2-\gamma_1^2)e^{-2\alpha_1d}
}{\Delta_2\cdot\Delta}.\qquad\text{(2010:6d)}
```

The shared dependency, 2010 (5c), is:

```math
F_{strat}(\lambda)=\mu_1
\frac{s_{12}+d_{12}e^{-2\alpha_1d}}
{s_{01}s_{12}+d_{01}d_{12}e^{-2\alpha_1d}}.\qquad\text{(2010:5c)}
```

The indexed ``P_{e_{ij}}^{-1}`` in (6a) is not independently replaced by an entrywise reciprocal or by new matrix notation. The explicitly stated matrix operation is (3a), above. ``M+jN`` is the dimensionless integral, whereas ``P_e`` has potential-coefficient units.

2010 original spectral definitions, p. 1169 below (A.3):

```math
a'_k=\sqrt{u^2+\gamma_k^2},\qquad
\gamma_k^2=j\omega\mu_k(\sigma_k+j\omega\varepsilon_k),\qquad k=0,1,2.
```

The complete parent reflection dependencies (A.6a)–(A.7g), with the one unprimed exponential in (A.7c) deliberately preserved:

```math
T'_1=\frac{\Delta'_1}{\Delta'}
=\frac{d'_{01}s'_{12}+s'_{01}d'_{12}e^{-2a'_1d}}
{s'_{01}s'_{12}+d'_{01}d'_{12}e^{-2a'_1d}},\qquad\text{(2010:A.6a)}
```

```math
T'_2=\frac{
\mu_0\mu_1(\gamma_0^2-\gamma_1^2)
[s'_{12}+d'_{12}e^{-2a'_1d}][S'_{12}+D'_{12}e^{-2a_1d}]
-4\mu_0\mu_1^2\mu_2a_1^{\prime2}\gamma_0^2e^{-2a_1d}
(\gamma_2^2-\gamma_1^2)
}{\Delta'_2\cdot\Delta'}.\qquad\text{(2010:A.6b)}
```

```math
\Delta'=s'_{01}s'_{12}+d'_{01}d'_{12}e^{-2a'_1d},\qquad\text{(2010:A.7a)}
```

```math
\Delta'_1=d'_{01}s'_{12}+s'_{01}d'_{12}e^{-2a'_1d},\qquad\text{(2010:A.7b)}
```

```math
\Delta'_2=S'_{01}S'_{12}+D'_{01}D'_{12}e^{-2a_1d},\qquad\text{(2010:A.7c)}
```

```math
s'_{mn}=(a'_m\mu_n+a'_n\mu_m),\qquad\text{(2010:A.7d)}
```

```math
d'_{mn}=(a'_m\mu_n-a'_n\mu_m),\qquad\text{(2010:A.7e)}
```

```math
S'_{mn}=(\mu_m\gamma_n^2a'_m+\mu_n\gamma_m^2a'_n),\qquad\text{(2010:A.7f)}
```

```math
D'_{mn}=(\mu_m\gamma_n^2a'_m-\mu_n\gamma_m^2a'_n),\qquad m,n=0,1,2.
\qquad\text{(2010:A.7g)}
```

The later (A.6b) also has mixed primed/unprimed exponent factors; their literal forms are retained above. Appendix (A.12) supplies the longitudinal-transform identity:

```math
\int_{-\infty}^{\infty}
J_0\!\left(u\sqrt{x^2+y_{ij}^2}\right)e^{-jk_0x}\,dx
=\begin{cases}
0,&u<k_0,\\
2\dfrac{\cos\!\left(y_{ij}\sqrt{u^2-k_0^2}\right)}
{\sqrt{u^2-k_0^2}},&u>k_0.
\end{cases}\qquad\text{(2010:A.12)}
```

No endpoint value at ``u=k_0`` is provided. Immediately following (A.12), the source prescribes ``u^2-k_0^2=\lambda^2`` and transforms ``a'_k`` to:

```math
\alpha_k=\sqrt{\lambda^2+\gamma_k^2+k_0^2},\qquad k=0,1,2,
\qquad T'_1\longmapsto T_1,\quad T'_2\longmapsto T_2.
```

The corresponding final coefficient dependency rendering, applying that prescription to the primed definitions, is:

```math
s_{mn}=(\alpha_m\mu_n+\alpha_n\mu_m),\qquad
d_{mn}=(\alpha_m\mu_n-\alpha_n\mu_m),
```

```math
S_{mn}=(\mu_m\gamma_n^2\alpha_m+\mu_n\gamma_m^2\alpha_n),\qquad
D_{mn}=(\mu_m\gamma_n^2\alpha_m-\mu_n\gamma_m^2\alpha_n),
```

```math
\Delta=s_{01}s_{12}+d_{01}d_{12}e^{-2\alpha_1d}.
```

For ``\Delta_2``, the source points to (A.7c) and the same transform, but its exponential is **already unprimed** ``a_1`` there, not the defined ``a'_1``. That specific exponent mapping remains unresolved; the full original (A.7c) is present above and no corrected final denominator is fabricated. The 2009 (A.6) is a distinct earlier witness, not permission to overwrite the 2010 prime convention.

The logarithmic contribution in (6b) uses the original definitions immediately below 2010 (A.13), p. 1170:

```math
D_{ij}=\sqrt{y_{ij}^2+(h_i+h_j)^2},\qquad
d_{ij}=\sqrt{y_{ij}^2+(h_i-h_j)^2}.
```

The associated identity is separately preserved because its sign and factors are suspect as printed:

```math
\int_0^\infty
\left(\frac{e^{-a'_0|h_j-h_i|}}{a'_0}
+\frac{e^{-a'_0(h_i+h_j)}}{a'_0}\right)
\cos(y_{ij}\lambda)\,d\lambda
=\ln\frac{D_{ij}}{d_{ij}}.\qquad\text{(2010:A.13)}
```

No minus is inserted, and no branch or convergence conclusion is imposed on this ambiguous post-transform notation. 2010 self substitution is applied to the complete (6): ``y_{ij}`` replaced by outer conductor radius, ``h_j`` by ``h_i``. Identical electromagnetic properties of the two earth layers give its stated homogeneous reduction ``\gamma_2=\gamma_1,a_2=a_1``; no new derived limiting kernel is substituted.

## Notation map

Notation is unchanged within each publication. Cross-publication correspondences are comparisons, not a normalization that resolves conflicting tokens. Units below are SI and equation-implied unless expressly stated by the source.

| Source symbol | Display symbol | Meaning | Units/convention |
| --- | --- | --- | --- |
| ``P_{g_{ij}},P_{pg_{ij}},P_{e_{ij}}`` | unchanged | Imperfect-earth correction, perfect-ground and full external potential coefficients | ``\mathrm{m/F}``; not admittances |
| ``Y'_{g_{ij}},Y'_{e_{ij}}`` | unchanged | Source-printed per-length admittance relations | ``\mathrm{S/m}``; indexed inverse notation preserved |
| ``\mathbf P_g,\mathbf P_{pg},\mathbf P_e,\mathbf Y'_g,\mathbf Y'_{pg},\mathbf Y'_e,\mathbf Y'`` | unchanged | Explicit 2010 ``N\times N`` potential/admittance matrices | Inversion only as supplied by (3); ``N`` here is conductor count |
| ``F,G,F_{strat},G_{strat}`` | unchanged | Impedance and radial-displacement kernels in each witness | Length; leading ``u``/``\lambda`` in ``G`` retained |
| ``M,N`` in 2010 (6c) | unchanged | Real/imaginary names of integral correction | Dimensionless; ``N`` is not conductor count in this equation |
| ``s,d,S,D`` with indices and primes | unchanged | Two distinct interface coefficient families | Lowercase coefficients proportional to permeability/length; uppercase additionally contain bulk-constant squares |
| ``\Delta,\Delta_2,\Delta',\Delta'_1,\Delta'_2,T'_1,T'_2,T_1,T_2`` | unchanged | Denominator combinations and parent/transformed reflection functions | Defined products/quotients; prime is source-stage notation, not differentiation |
| ``a_m,\alpha_1,a'_1`` (2009) | unchanged | Defined Latin spectral factors versus other printed factors | ``\mathrm{m^{-1}}``; ``\alpha_1`` and ``a'_1`` not separately identified |
| ``a'_k,\alpha_k,a_0,a_1`` (2010) | unchanged | Original, transformed, and mixed main/appendix factors | ``\mathrm{m^{-1}}``; mixed-prime occurrences remain unresolved |
| ``u`` (2009), ``\lambda`` (2010 final), ``u`` (2010 parent) | unchanged | Final or original spectral integration variables | ``\mathrm{m^{-1}}``; source-specific roles retained |
| ``\gamma_m,\gamma_k,\gamma_x,k_0`` | unchanged | Bulk-medium constants; later longitudinal parameter and air wavenumber | ``\mathrm{m^{-1}}``; imposed ``\gamma_x`` distinct from solved TL propagation |
| ``\mu_m,\varepsilon_m,\sigma_m`` | unchanged | Layer scalar permeability, permittivity, conductivity | ``\mathrm{H/m},\mathrm{F/m},\mathrm{S/m}``; media 0,1,2 |
| ``h_i,h_j,y_{ij},d,D_{ij},d_{ij},x,z`` | unchanged | Heights, horizontal separation, layer thickness, image/direct distances and longitudinal/upward coordinates | Metres; ``d`` and ``d_{ij}`` are different quantities |
| ``J_0,j,\omega,f,c,h`` | unchanged | First-kind order-zero Bessel function; imaginary unit; angular frequency; frequency; light speed; height in TL criterion | ``\omega=2\pi f``; no explicit time exponential |
| ``f_{cr},f_{cr,min},\varepsilon_{r1}`` | unchanged | Source-attributed critical frequency, one-tenth threshold, and first-earth-layer relative permittivity | Hz, Hz, dimensionless; earth-behaviour classification, not a universal kernel bound |

## Evidence and approximation sources

The 2009 main ``G`` and 2010 final ``G_{strat}`` have the same two numerator contributions after matching the final spectral variables ``u`` and ``\lambda``. The literal dependency conflicts prevent claiming a completely resolved common evaluable expression. They are retained as separately labelled witnesses in one source-attributed formulation record rather than counted as a new contribution from variable renaming. The 2010 main/full-matrix assembly, parent integral, transform and appendix are additional source evidence, not a repaired 2009 edition.

The source-stated approximation is quasi-TEM TL-mode propagation. The original article cites the Hertz-vector methodology to Sunde; the later paper explicitly imposes the free-space longitudinal constant on (4b), applies (A.12), and includes the resulting ``G`` term. No Taylor/Bessel asymptotic order, finite series, ignored-tail magnitude or general error bound is supplied. The 2009 paper expressly says omitting ``G`` produces purely imaginary propagation and lossless propagation over imperfect earth; 2010 repeats it on p. 1162. This remains an **author assertion**, not a reviewer-proved consequence of discarding a kernel or a basis for changing the formula.

2009 IV–V tests the horizontal 150 kV three-phase geometry (10 m height, 6.5 m adjacent spacing), earth resistivity ratios 10,5,0.1,0.2 with lower resistivities 10,100,1000 ``\mathrm{\Omega\,m}``, relative permittivities 1–20, upper-layer depth 5–20 m, and unit relative permeability. The plotted cases extend into MHz frequencies. 2010 expressly studies 50 Hz–10 MHz. Its model comparisons and transient simulations test differences among selected assumptions, not exactness against a universally valid full-wave solution. No numerical implementation or regression test was made for this transcription.

The source's earth-behaviour classification, 2009 (3)/2010 (7), is ``f_{cr}=\sigma_1/(2\pi\varepsilon_0\varepsilon_{r1})`` in Hz, with ``f_{cr,min}=0.1f_{cr}``: conducting below this minimum, comparable displacement/resistive currents between ``0.1f_{cr}`` and ``2f_{cr}``, and displacement-dominated above ``2f_{cr}``. This is attributed to Semlyen, not a new independent admittance formulation or a universal admissibility bound. Boundary equality cases are not assigned by the printed strict inequalities.

The authors identify their two-layer admittance expression as similar to Ametani–Nagaoka–Koide's 2001 overhead-conductor-above-snow model. That English translation and its 2000 Japanese original are transcribed in the [snow-layer record](../../2000/snow-layer-overhead-potential-coefficient/Ametani2001.md). Its imposed air propagation, potential-coefficient output, and finite surface layer are documented separately, with the original/translation coefficient and constitutive conflicts retained. The authors' similarity statement does not establish full algebraic equivalence or permit importing this later self prescription into the snow original. The [impedance record](../../../external-impedance/2009/two-layer-earth-overhead-integral-correction/Papadopoulos2009.md) documents the source's differing Nakagawa citation identities and the distinction between earth correction and complete external impedance.

## Limitations and discrepancies

1. **Original definition gaps:** Greek ``\alpha_1`` in (1b)/(2c) versus Latin ``a_m`` below the appendix, and primed ``a'_1`` only in (A.5), are actually printed. The main ``F`` denominator and the separately printed ``\Delta`` therefore do not have an explicitly identical exponent notation. None is corrected.
2. **Later mixed-prime factors:** 2010 (A.6b) mixes ``a'_1`` and ``a_1`` in exponential factors and (A.7c) prints unprimed ``a_1``. The announced ``a'_k\to\alpha_k`` transform does not explicitly resolve those already-unprimed occurrences. The record keeps the full parent coefficients and marks the particular final ``\Delta_2`` dependency unresolved.
3. **Potential versus admittance:** original (2a) and later (6a) have indexed inverse notation; later (3a) explicitly provides matrix inversion of a sum of potential matrices. No entrywise interpretation or parallel sum of correction admittances is inferred. ``P_g`` alone is not full external ``P_e``.
4. **Suspected published identity defect:** 2010 (A.13) prints a plus between direct and image exponentials with post-transform primed factors. The image, not OCR, confirms it. No corrected sign, substituted spectral factor or isolated-term singularity claim is supplied.
5. **Conventions and parent domain:** root branches and an original explicit longitudinal exponent are absent. Later free-space propagation explains its own final kernel only; the parent complex-``\gamma_x`` convergence domain is not supplied. No zero-propagation or restored displacement-current assumption is invented.
6. **Bibliographic and priority limits:** the 2009 printed third-author name conflicts with the faithfully copied existing entry. The two publications and their equation locators remain separate. The snow predecessor is independently verified in its own 2000/2001 record, not retrospectively certified from this paper; full equivalence and earliest priority are not established. No arbitrary-layer recursion is extracted from extension prose, and no separately published erratum has been inspected.
