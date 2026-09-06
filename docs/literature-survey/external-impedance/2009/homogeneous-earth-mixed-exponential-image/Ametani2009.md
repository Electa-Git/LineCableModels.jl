# Ametani–Yoneda–Baba–Nagaoka mixed-conductor approximation

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel horizontal isolated conductors with heights/depths ``h_1,h_2`` and transverse separation ``y``. Their radii are drawn in Fig. 1 but do not enter this mutual earth term. |
| Calculated quantities | Mutual earth-return impedance per unit length between one overhead and one buried conductor |
| Earth structure | Homogeneous conducting half-space below homogeneous air with a planar interface. |
| Model and approximation | The parent is the Pollaczek mixed integral as restated in (5)–(6). The authors first substitute ``\sqrt{s^2+m^2}\simeq s+m`` in the exponent, (17), then differentiate with respect to ``y``, change variables ``s=mt``, deform the contour under a stated no-pole assumption, and substitute ``t/(\sqrt{t^2+1}+t)\simeq(1-e^{-2t})/2``, (24). Integration and the subsequent antiderivative in ``y`` lead to (27). No controlled retained/discarded series order or uniform expansion parameter is assigned to these two substitutions. A further power-frequency approximation is recorded separately. |
| Main source | Akihiro Ametani, Tetsuzo Yoneda, Yoshihiro Baba, and Naoto Nagaoka (2009) |
| Citation key(s) | `:Ametani2009` |
| Evidence status | Original PDF equations checked against page images; source coordinate convention unresolved. |

**Description.** Exponential-image approximation to the mutual earth-return impedance between parallel horizontal conductors on opposite sides of a homogeneous air–earth interface.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | The parent is described as a TEM-mode approximation; no independent imposed longitudinal constant appears in (5)–(6) or (27). Its longitudinal phase prescription is not stated. | Stated — sections II–III, pp. 860–862; equation-implied — (5)–(6), (27). |
| Air propagation constant ``γ_air`` | ``m_1=0`` in the reduced air model: zero air conductivity and omitted displacement current. | Stated/equation-implied — (4), section II-B, and the ``|s|`` air kernel in (6), p. 861. |
| Earth propagation constant ``γ_earth`` | ``m=\sqrt{j\omega\mu_0/\rho_e}=|m|\exp(j\pi/4)``, with ``h_e=1/m``. | Stated — (22) and definitions following (27), p. 864. |
| Earth permittivity and displacement current | Omitted from this approximation. The more general medium definition in (3) includes ``j\omega\varepsilon_i``, but is reduced to conduction-only (4); it must not be substituted back into (27). | Stated — (3)–(4) and section III, p. 861; (22), p. 864. |
| Range of validity | Restricted to the parent TEM approximation. The authors propose a burial-depth condition ``h\leq\lambda/8`` and a critical-frequency formula (15), whose printed consistency is unresolved. The reported errors below 10% overall and below 3% at low frequency refer to the Fig. 10 examples, not a universal bound. | Stated — (11), p. 862; (15), p. 863; section IV-B and Fig. 10, p. 865. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0``. | Stated — definitions after (3), p. 861; (22), (27). |
| Arrangement | Mixed mutual term: one overhead conductor and one buried conductor; no self formula in this new approximation. | Stated — title, Fig. 1, (5), and section IV. |
| Earth structure | Homogeneous conducting half-space below homogeneous air with a planar interface. | Stated — Fig. 1 and section II, p. 861. |
| Conductor and insulation geometry | Parallel horizontal isolated conductors with heights/depths ``h_1,h_2`` and transverse separation ``y``. Their radii are drawn in Fig. 1 but do not enter this mutual earth term. | Stated/equation-implied — Fig. 1 and (27). |
| Constitutive and field assumptions | Scalar uniform resistivity, nonmagnetic media, conduction-dominated earth, and the source's TEM assumption. Internal conductor and insulation contributions are separate from this output. | Stated/equation-implied — sections II–IV and (5). |
| Conventions | ``j`` is the imaginary unit; ``\omega=2\pi f``. Fig. 1 uses an upward ``z`` axis and section III writes ``h_2=-h<0``, but section IV uses ``\exp(-h_2m)`` and sums ``h_1+h_2`` without explicitly changing the depth convention. The record retains this ambiguity. ``H`` is redefined after (27). Output is per unit length, plotted in ``\Omega/\mathrm m``. | Stated — Fig. 1, section III-B, definitions after (18) and (27), Fig. 9. |

**Expression.** Source equation (27), including the definitions immediately below it.

```math
Z_m=j\omega\left(\frac{\mu_0}{2\pi}\right)
\exp\left(\frac{-h_2}{h_e}\right)\ln\left(\frac{S}{D}\right),
\qquad\text{(27)}
```

```math
\begin{aligned}
h_e&=\frac{1}{m} \\
S&=\sqrt{H^2+y^2} \\
D&=\sqrt{(h_1+h_2)^2+y^2} \\
H&=h_1+h_2+2h_e, \\
m&=\sqrt{\frac{j\omega\mu_0}{\rho_e}}
=|m|\exp\left(\frac{j\pi}{4}\right).
\end{aligned}\qquad\text{(22)}
```

``\rho_e`` is earth resistivity in ``\Omega\,\mathrm m``; ``h_e`` is a complex length. ``h_1,h_2,y,S,D,H`` have dimensions of length. The paper does not separately prescribe branches for ``S`` and the logarithm. In the derivation before (27), ``H`` instead means ``h_1+h_2``; the two definitions are kept at their original locations below.

**Approximation.** The parent is the Pollaczek mixed integral as restated in (5)–(6). The authors first substitute ``\sqrt{s^2+m^2}\simeq s+m`` in the exponent, (17), then differentiate with respect to ``y``, change variables ``s=mt``, deform the contour under a stated no-pole assumption, and substitute ``t/(\sqrt{t^2+1}+t)\simeq(1-e^{-2t})/2``, (24). Integration and the subsequent antiderivative in ``y`` lead to (27). No controlled retained/discarded series order or uniform expansion parameter is assigned to these two substitutions. A further power-frequency approximation is recorded separately.

**Limitations.** The printed depth convention and reuse of ``H`` prevent an unambiguous single-coordinate interpretation of the entire derivation. The source's critical-frequency and accuracy statements apply to its stated cases. The model neglects earth displacement current and does not provide multilayer or self-interaction extensions of the new mixed approximation.

**Reference.** [Ametani2009](@cite), (5)–(6), p. 861; (17)–(27), p. 864; section IV-B, p. 865. The existing bibliography entry conflicts with the four-author byline on p. 860.

**Transcription source.** Original publication PDF images, pp. 860–865. The source equation, exponential sign, logarithmic distances, explicit branch of ``m``, parent integral, and approximation substitutions were visually checked. No equation-preserving Markdown was used as authority.

## Source transcription

Source parent, attributed in this article to Pollaczek (not newly attributed to the 2009 authors):

```math
Z_m=Z(1,2)=j\omega\left(\frac{\mu_0}{2\pi}\right)
\int_{-\infty}^{\infty}F_c(s)\exp(jys)\,ds,
\qquad\text{(5)}
```

```math
F_c(s)=\frac{\exp\{-h_1|s|+h_2\sqrt{s^2+m^2}\}}
{\sqrt{s^2+m^2}+|s|}.
\qquad\text{(6)}
```

Here ``s`` is a real spectral integration variable, with units of inverse length; in this parent, buried ``h_2`` is negative. The authors write the following substitutions in section IV:

```math
\sqrt{s^2+m^2}\simeq s+m,
\qquad\text{(17)}
```

```math
\begin{aligned}
H&=h_1+h_2>0,\qquad\text{definition following (18)}, \\
s&=mt \\
ds&=m\,dt,
\end{aligned}\qquad\text{(22, change of variable)}
```

```math
\frac{t}{\sqrt{t^2+1}+t}\simeq\frac{1-e^{-2t}}{2}.
\qquad\text{(24)}
```

The source then prints

```math
Z_m'=-j\omega\left(\frac{\mu_0}{2\pi}\right)e^{-h_2m}
\left[\frac{y}{y^2+H^2}-\frac{y}{y^2+(H+2/m)^2}\right].
\qquad\text{(26)}
```

``Z_m'`` is the derivative with respect to ``y``. The final (27) and its redefined ``H`` are reproduced above without notation changes. Equation (28) is explicitly a secondary transcription of Lucca's 1994 approximation.

## Notation map

Notation is unchanged. No sign or positive-depth normalization has been applied.

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_m,Z(1,2)`` | unchanged | Mixed mutual earth-return impedance | ``\Omega/\mathrm m`` |
| ``Z_m'`` | unchanged | Derivative of ``Z_m`` with respect to separation ``y`` | ``\Omega/\mathrm m^2`` |
| ``h_1,h_2`` | unchanged | Overhead height and buried coordinate/depth | m; conflicting sign treatment noted above |
| ``h=-h_2`` | unchanged | Positive burial depth used in section III | m |
| ``y`` | unchanged | Horizontal separation | m |
| ``\rho_e,\mu_0`` | unchanged | Earth resistivity and vacuum permeability | ``\Omega\,\mathrm m``, H/m |
| ``m,h_e`` | unchanged | Earth propagation constant and its reciprocal | inverse m, m |
| ``H`` | unchanged | ``h_1+h_2`` in (18)–(26); ``h_1+h_2+2h_e`` after (27) | m; local redefinition |
| ``S,D`` | unchanged | Complex image distance and geometric distance in (27) | m |
| ``s,t`` | unchanged | Original spectral variable and variable after ``s=mt`` | inverse m, dimensionless |
| ``F_c`` | unchanged | Parent mixed-interface kernel | m |
| ``j,\omega,f`` | unchanged | Imaginary unit, angular frequency, frequency | dimensionless, rad/s, Hz |
| ``\lambda`` | unchanged | Approximate soil wavelength in the source validity discussion | m |

## Evidence and approximation sources

The new contribution is (27), not the earlier Pollaczek or Lucca kernels reproduced for comparison. The Fig. 10 tests vary overhead heights 5, 10, 50 m, burial depths 0.1, 1, 3 m, resistivities 10, 100, 1000 ``\Omega\,\mathrm m``, and separations 0, 10, 50 m in the captioned geometries. The graphs extend to high frequencies where their errors exceed the prose's blanket 10% statement. Error (29) is printed as ``100(Z_{\mathrm{app}}-Z_P)/Z_P``; the panels separately label resistance and inductance errors.

The article states ``h\leq\lambda/8`` in (11), with ``\lambda\simeq2\pi\sqrt{2\rho_e/(\omega\mu_0)}`` in (10), and prints ``f_0=1.407\times10^5(\rho_e/h)`` in (15). These are author-stated claims, with an apparent inconsistent dependence on ``h``; no corrected critical frequency is supplied here.

## Limitations and discrepancies

- Source convention conflict: (6) uses the negative buried coordinate; (17)–(27) use a decaying factor for a positive depth without declaring the change. Existing LCM uses positive height/depth magnitudes; that implementation is not authority for repairing this source.
- Directly observable source redefinition: ``H`` gains ``2h_e`` only in the definitions after (27).
- Suspected published defect: (26) lacks the factor needed to obtain the printed (27) by an antiderivative as claimed. Both expressions are retained exactly.
- Source layout defects: (18), (20), (21), and (23) place the integral sign in the numerator of a displayed fraction and the differential in its denominator. Those intermediate display layouts are not silently rewritten as purportedly source-exact integral equations; the well-formed parent (5)–(6) and actual kernel operations are retained above. The prose immediately before (21) also prints identical exponential terms in its sine identity.
- Suspected source inconsistency: the claimed frequency criterion (15) has a different burial-depth dependence from the immediately preceding ``ah/\sqrt2=\pi/4`` with ``a^2=\omega\mu_0/\rho_e``. No repair is made.
