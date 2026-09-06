# de Lima–Portela buried impedance with frequency-dependent soil

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Cable external radius ``r`` for self evaluation; burial depths ``h_i,h_j``; horizontal separation ``d_{ij}``; direct and reflected distances ``d,D,D_c``. |
| Calculated quantities | Buried-cable self and mutual ground-return impedance per length |
| Earth structure | Homogeneous isotropic linear half-space under air. |
| Model and approximation | Integral representation within the source's reduced field model, combined with an empirical joint soil law. No analytical approximation of these integrals is introduced; Gauss–Kronrod is their numerical evaluation method. The stated conduction-only limit sets ``\epsilon_s=0`` and ``\sigma_s=\sigma_0``, so ``\eta_s=\sqrt{j\omega\mu_0\sigma_0}``. |
| Main source | Antonio Carlos Siqueira de Lima and Carlos Portela (2007), extending the Pollaczek-type earth representation to complex frequency-dependent soil |
| Citation key(s) | `:DeLima2007` |
| Evidence status | Original PDF equations checked against page images; main-text/appendix Bessel-order disagreement and appendix dependencies unresolved |

**Description.** Self and mutual earth-return impedances of horizontal cables buried in homogeneous earth with frequency-dependent conductivity and permittivity, expressed as Bessel terms and separate interface integrals.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Appendix specifies ``\exp(-kz)`` with unknown ``k``; no ``k`` remains in (3)–(4). The source does not state a corresponding imposed ``Γ=0`` prescription. | Stated — appendix opening, p. 497; equation-implied — (3)–(4). |
| Air propagation constant ``γ_air`` | No air bulk propagation constant in final interface kernel; appendix uses ``\nabla^2E_a=0``. Its reduction from the stated longitudinal wave dependence is unresolved. | Equation-implied — (3)–(4), p. 493; (15), p. 497. |
| Earth propagation constant ``γ_earth`` | Main text: ``\eta_s=\sqrt{j\omega\mu_0(\sigma_s+j\omega\epsilon_s)}``. Appendix prints an inconsistent ``\eta_g`` definition without the radical. | Stated — definition after (2), p. 493; after (16), p. 497. |
| Earth permittivity and displacement current | Retained together with frequency-dependent conductivity in ``\kappa'\simeq\sigma_s+j\omega\epsilon_s``. | Stated — (7)–(8), p. 493, and section II-B, p. 494. |
| Range of validity | TEM/quasi-TEM treatment; approximately 1 MHz discussed as dependent on geometry, soil and line length, not a universal cutoff. Soil fitting through 2 MHz is a different evidential claim. | Stated — section II opening, p. 493; section II-B, p. 494. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0`` in main (3)–(4). Appendix boundary conditions retain medium permeabilities, but they do not establish arbitrary permeability in the main expressions. | Equation-implied — (3)–(4); (18), p. 497. |
| Arrangement | Underground self and mutual, parallel horizontal cables; distinct depths retained in the mutual distance. | Stated — section II-A, p. 493. |
| Earth structure | Homogeneous isotropic linear half-space under air. | Stated — section II opening and appendix. |
| Conductor and insulation geometry | Cable external radius ``r`` for self evaluation; burial depths ``h_i,h_j``; horizontal separation ``d_{ij}``; direct and reflected distances ``d,D,D_c``. | Stated — definitions following (4), p. 493. |
| Constitutive and field assumptions | Source TEM/quasi-TEM approximation, linear isotropic homogeneous media, and scalar frequency-dependent soil response. Conductor/insulation terms are assembled separately using earlier cable formulations. | Stated — section II and section III, pp. 493–494. |
| Conventions | Time dependence ``\exp(j\omega t)``, longitudinal dependence ``\exp(-kz)``. Main burial depths positive; appendix uses negative ``y`` below the interface and an opposed depth coordinate ``t``. Ground impedance defined by ``Z_g=-E_g/I``. | Stated — appendix opening, definitions after (16), and (32), pp. 497–498. |

**Expression.** Main-text equations (3)–(4), preserving the printed ``K_1`` image term and both independent integration domains.

```math
z_{p_{ii}}=\frac{j\omega\mu_0}{2\pi}
\left[K_0(\eta_s r)-K_1(\eta_sD_c)
+\int_{-\infty}^{\infty}
\frac{\exp\left(-2h_i\sqrt{\xi^2+\eta_s^2}\right)}
{|\xi|+\sqrt{\xi^2+\eta_s^2}}
\exp(jr\xi)\,d\xi\right],
\tag{3}

z_{p_{ij}}=\frac{j\omega\mu_0}{2\pi}
\left[K_0(\eta_s d)-K_1(\eta_sD)+{}
+\int_{-\infty}^{\infty}
\frac{\exp\left(-(h_i+h_j)\sqrt{\xi^2+\eta_s^2}\right)}
{|\xi|+\sqrt{\xi^2+\eta_s^2}}
\exp(jd_{ij}\xi)\,d\xi\right],
\tag{4}

d=\sqrt{d_{ij}^2+(h_i-h_j)^2},\qquad
D_c=\sqrt{r^2+4h_i^2},\qquad
D=\sqrt{d_{ij}^2+(h_i+h_j)^2},

\eta_s=\sqrt{j\omega\mu_0(\sigma_s+j\omega\epsilon_s)}.
```

Equation (4) actually prints a plus at the end of the first line and another before the integral; both are retained above. ``r,h_i,h_j,d_{ij},d,D,D_c`` are in meters, ``\xi,\eta_s`` inverse meters, and the outputs are Ω/m. The main text calls ``K_0,K_1`` Bessel functions without explicitly specifying their kind; this dependency remains unresolved rather than silently inferred from the notation.

The source's constitutive dependency is

```math
\sigma_s+j\omega\epsilon_s\simeq\kappa'
=\sigma_0+\delta_{\sigma_s}+j\delta_{\omega\epsilon_s},
\tag{7}

\delta_{\sigma_s}+j\delta_{\omega\epsilon_s}
=\Delta_i\left(\frac{f}{10^6}\right)^\alpha
\left(\cot(\alpha\pi/2)+j\right).
\tag{8}
```

``\sigma_0`` is low-frequency soil conductivity, ``\Delta_i`` a fitted amplitude in S/m, and ``\alpha`` the dimensionless fitted exponent; ``\omega=2\pi f``. ``\sigma_s`` is in S/m, ``\epsilon_s`` F/m, and ``\kappa'`` S/m. The source's conductor/insulation assembly refers to Wedepohl–Wilcox and Ametani (1980); no new insulation equation is supplied here.

**Approximation.** Integral representation within the source's reduced field model, combined with an empirical joint soil law. No analytical approximation of these integrals is introduced; Gauss–Kronrod is their numerical evaluation method. The stated conduction-only limit sets ``\epsilon_s=0`` and ``\sigma_s=\sigma_0``, so ``\eta_s=\sqrt{j\omega\mu_0\sigma_0}``.

**Limitations.** The ``K_1`` main-text image terms disagree with the appendix's ``K_0`` image term. The appendix also has unresolved definition and exponential-layout defects; its alternative witness is reproduced below without using it to correct (3)–(4). The paper does not publish a matching frequency-dependent ground-admittance correction.

**Reference.** [DeLima2007](@cite), equations (3)–(4), (7)–(8), p. 493; appendix (30)–(32), p. 498.

**Transcription source.** Original PDF page images, pp. 493–494 and 497–498. Both Bessel subscripts, each exponential, the real-line integral domains and absolute-value denominator, and all geometric definitions were checked visually. The main-text/appendix disagreement is in the printed publication, not merely in conversion text.

## Source transcription

The main source equations and geometry are retained above. The appendix prints an alternative field witness and a source-prescribed conversion to impedance:

```math
E_g=-\frac{j\omega\mu_g I}{2\pi}
\left(K_0(\eta_gD)-K_0(\eta_gD')+2\chi\right),
\tag{30}

\chi=\int_0^\infty
\frac{\exp(y-h)\sqrt{\alpha^2+\eta^2}}
{\dfrac{\mu_g}{\mu_a}\alpha+\sqrt{\alpha^2+\eta^2}}
\cos(\alpha x)\,d\alpha,
\tag{31, printed layout}

D=\sqrt{x^2+(y+h)^2},\qquad
D'=\sqrt{x^2+(y-h)^2},\qquad
Z_g=-\frac{E_g}{I}.
\tag{32 and adjacent definitions}
```

Here ``h`` is positive burial depth, ``y`` is negative in the ground, and the text states ``x_p=r,y_p=-h`` for the self evaluation. The appendix's ``D`` is not the main text's reflected distance ``D``: the local definitions are retained separately. The radical in printed (31) is outside ``\exp(y-h)``; it has not been moved into its argument. The bare ``\eta`` appearing there is not defined consistently with ``\eta_g`` or main-text ``\eta_s``. This appendix witness is transcription-verified but mathematically unresolved.

## Notation map

No renaming or coordinate transformation is applied. Main and appendix symbols with conflicting uses remain local to their displayed equations.

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``z_{p_{ii}},z_{p_{ij}}`` | unchanged | Underground self and mutual ground impedances | Ω/m |
| ``r,h_i,h_j,d_{ij}`` | unchanged | External radius, depths and horizontal separation | m; main depths positive |
| ``d,D,D_c`` in (3)–(4) | unchanged | Main direct, reflected and self-image distances | m |
| ``\eta_s,\xi`` | unchanged | Main soil propagation constant and spectral variable | inverse m |
| ``K_0,K_1`` | unchanged | Source's Bessel functions | kind not explicitly defined in main text |
| ``\sigma_s,\epsilon_s,\mu_0`` | unchanged | Soil conductivity, permittivity, permeability | S/m, F/m, H/m |
| ``\kappa',\sigma_0,\delta_{\sigma_s},\delta_{\omega\epsilon_s},\Delta_i,\alpha`` in (7)–(8) | unchanged | Joint soil response and fitted parameters | S/m except dimensionless exponent ``\alpha`` |
| ``E_g,I,Z_g`` | unchanged | Ground longitudinal electric field, current, converted ground impedance | V/m, A, Ω/m |
| ``D,D',x,y,h`` in appendix | unchanged | Appendix direct/image geometry and coordinates | m; ``y<0`` below ground |
| ``\alpha`` in (31) | unchanged | Appendix spectral variable, not soil exponent in (8) | inverse m |
| ``\eta_g,\eta,\chi`` in appendix | unchanged | Appendix propagation symbols and interface integral | inconsistent definitions/dimensions unresolved |
| ``\mu_a,\mu_g`` | unchanged | Appendix air/ground permeabilities | H/m |
| ``j,\omega,f,k,z`` | unchanged | Imaginary unit, angular frequency, frequency, unknown longitudinal propagation function, axial coordinate | dimensionless, rad/s, Hz, inverse m, m |

## Evidence and approximation sources

The overhead companion record contains the [sample soil fit and source assembly context](../homogeneous-earth-overhead-frequency-dependent-soil/DeLima2007.md). The paper calls (1)–(4) extensions of earlier Carson/Pollaczek expressions with complex soil parameters. Its material fit comes from earlier Portela studies; neither a quadrature rule nor a later soil-model implementation is substituted for the printed mathematical expression.

## Limitations and discrepancies

- Directly observable source conflict: main (3)–(4) use ``K_1`` for the image term; appendix (30) uses ``K_0``. Both are preserved, with no selected correction.
- Equation (4) repeats a plus at the line break; (3) does not.
- Appendix (31) prints the radical outside the exponential, and switches from ``\eta_g`` to undefined ``\eta``. The corresponding main-text exponent is different. No reconstruction of the appendix kernel is made.
- Appendix (16) and the definition immediately below it have symbol/definition inconsistencies, including a propagation definition without the main-text radical. Those lines are not used to replace ``\eta_s``.
- The unknown longitudinal ``k`` is introduced but not explicitly removed or prescribed. The final formula alone does not justify reporting an author-stated ``Γ=0``.
- The main text's generic Bessel-function identification and missing square-root branch leave dependencies unresolved even though the printed expressions themselves are image-verified.
