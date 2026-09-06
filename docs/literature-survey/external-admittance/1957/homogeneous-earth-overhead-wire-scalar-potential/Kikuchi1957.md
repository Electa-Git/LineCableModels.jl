# Kikuchi overhead-wire scalar potential above homogeneous earth

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite thin circular overhead conductor of radius ``a`` and height ``h``; no insulation. The internal current distribution is approximated by that of an isolated conductor. |
| Calculated quantities | Scalar potentials ``V_1(x,y)`` in air and ``V_2(x,y)`` in earth due to one overhead wire; source-defined self terminal voltage ``V=V_1(0,h-a)``; charge/current relation retained separately |
| Earth structure | One flat homogeneous half-space ``y<0`` below air ``y>0``; surface ``y=0``. |
| Model and approximation | Hankel-plus-integral representation under the thin-wire, homogeneous-medium, and principal-mode assumptions. The small-Hankel and ``\Gamma\simeq jk_1`` reductions are auxiliary evaluators. |
| Main source | H. Kikuchi (1957); the 1955 Japanese paper is an earlier-language witness, and priority remains unresolved. |
| Citation key(s) | `:Kikuchi1957` |
| Evidence status | All 13 original pages checked; the auxiliary reductions are not independent physical formulations. |

**Description.** Scalar electric potential produced by one infinitely long thin cylindrical conductor above a plane homogeneous earth, with observations in air or earth and a source-defined wire terminal-voltage reference. The representation retains longitudinal dependence and earth conductivity and permittivity. It is recorded in the external-admittance family as potential evidence, not relabelled as a self/mutual admittance matrix.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | ``\Gamma`` remains in the potential prefactor and ``\lambda_n^2=k_n^2+\Gamma^2``; dependence is ``e^{-\Gamma z+j\omega t}``. The author selects the propagating principal mode. An imposed zero or numerical modal closure is not supplied for this general equation set. The later prescription ``\Gamma\simeq jk_1`` belongs to (5·9), not the unreduced formulas recorded here. | Stated — §3, p. 722; (5·1)–(5·4), pp. 723–724; distinct (5·9), p. 725. |
| Air propagation constant ``γ_air`` | Source symbol is ``k_1``, not ``\gamma``. The common definition is ``k_n=\omega\sqrt{\mu_0(\varepsilon_n-j\sigma_n/\omega)}``. The scalar air source equation uses ``\varepsilon_1`` and ``Q=\Gamma I/(j\omega)``. No nonzero air-conductivity extension or replacement of ``\varepsilon_1`` with a complex permittivity is introduced. Air is region 1; the symbol list identifies vacuum/air constants ``\varepsilon_0,\mu_0``. | Stated — symbol list, p. 721; (3·5), definition below (3·7) and Fig. 1, p. 722. Equation-implied — (3·5) and (5·1) use the lossless-air scalar-source factor; an independently nonzero ``\sigma_1`` is not resolved by that factor. |
| Earth propagation constant ``γ_earth`` | Source symbol ``k_2=\omega\sqrt{\mu_0(\varepsilon_2-j\sigma_2/\omega)}``; ``\lambda_2=\sqrt{k_2^2+\Gamma^2}`` and ``\kappa_2^2=u^2-\lambda_2^2``. These are three different quantities. ``\mathop{\mathrm{Re}}\kappa_2>0`` is printed. | Stated — §§2–4, pp. 721–722, especially (4·4). |
| Earth permittivity and displacement current | Both ``\varepsilon_2`` and ``\sigma_2`` remain in ``k_2`` and hence the material-weighted denominator. No extra dielectric loss or frequency model is added. | Stated — definition below (3·7), p. 722; Equation-implied — (5·3)–(5·4), pp. 723–724, retain ``k_2^2``. |
| Range of validity | Conductor radius is sufficiently small compared with height that the earth's influence on internal field/current distribution is neglected and the exterior source is a line current. Only the propagating principal mode is treated. No universal numerical frequency interval or explicit error bound is provided for this general set. Later small-Hankel-argument and near-air-speed approximations are not applied here. | Stated — end of §3 opening, pp. 721–722; §11, p. 732; separate approximation conditions below (5·8), p. 724, and (5·9), p. 725. |
| Earth permeability ``μ_earth`` | The source explicitly writes ``\mu_1=\mu_2`` and uses ``\mu_0`` in ``k_n``. This is not an arbitrary-permeability earth expression. | Stated — definition split between columns at the bottom/top of p. 722. |
| Arrangement | One overhead source wire at ``(x,y)=(0,h)``. Air and earth observation points are both supplied. The self terminal prescription evaluates the lower wire surface ``(0,h-a)``. There is no buried source or multiple-wire admittance assembly in this set. | Stated — Fig. 1, p. 722; (5·1)–(5·4), pp. 723–724; §9, p. 729. |
| Earth structure | One flat homogeneous half-space ``y<0`` below air ``y>0``; surface ``y=0``. | Stated — Fig. 1 and boundary conditions (3·6)–(3·7), p. 722. |
| Conductor and insulation geometry | Infinite cylindrical conductor, radius ``a``, centre height ``h``, with radius small relative to height. No finite insulation layer is supplied. Internal current distribution is assumed the same as for the isolated conductor and represented externally by total line current ``I``. | Stated — pp. 721–722, §3 and Fig. 1. |
| Constitutive and field assumptions | Scalar homogeneous material parameters; harmonic Maxwell fields and coupled scalar/vector potentials. Conductor-internal earth/proximity perturbation is neglected. The general scalar potential is not obtained by declaring a pure TEM field: the paper discusses hybrid-mode corrections to circuit power. | Stated — §§3–5, pp. 721–724; §8, pp. 728–729; §§11 and Appendix B, pp. 732–733. |
| Conventions | Rationalized MKS units; ``j^2=-1``; ``e^{j\omega t}`` and ``e^{-\Gamma z}``; vertical ``y`` positive into air, longitudinal ``z`` along the wire. ``Q`` is charge per wire length, ``I`` total wire current. The terminal potential uses an infinitely remote return reference, not a silently imposed zero potential at every earth-surface point. | Stated — symbols and §3, p. 721; Fig. 1 and (3·5), p. 722; §9, p. 729. |

**Expression.** Scalar formulas of source (5·1)–(5·4), pp. 723–724. Primes on ``Q'_n,P'_n`` are part of the auxiliary names, not derivatives.

```math
V_1(x,y)=\frac{\Gamma I}{j\omega\varepsilon_1\pi}
\left[
\frac{j\pi}{4}
\left\{H_0^{(1)}(\lambda_1\rho)-H_0^{(1)}(\lambda_1\rho_0)\right\}
+(Q'_1-jP'_1)
\right],
```

```math
Q'_1-jP'_1
=k_1^2\int_0^\infty
\frac{e^{-(y+h)\sqrt{u^2-\lambda_1^2}}\cos xu}
{k_2^2\sqrt{u^2-\lambda_1^2}+k_1^2\sqrt{u^2-\lambda_2^2}}\,du,
```

```math
V_2(x,y)=\frac{\Gamma I}{j\omega\varepsilon_1\pi}(Q'_2-jP'_2),
```

```math
Q'_2-jP'_2
=k_1^2\int_0^\infty
\frac{e^{\sqrt{u^2-\lambda_2^2}\,y-\sqrt{u^2-\lambda_1^2}\,h}\cos xu}
{k_2^2\sqrt{u^2-\lambda_1^2}+k_1^2\sqrt{u^2-\lambda_2^2}}\,du.
```

Source definitions, §§2–4, pp. 721–722:

```math
\begin{aligned}
k_n&=\omega\sqrt{\mu_0\left(\varepsilon_n-j\frac{\sigma_n}{\omega}\right)}
=k_n^{(r)}+jk_n^{(i)} \\
\mu_1&=\mu_2 \\
\lambda_n&=\sqrt{k_n^2+\Gamma^2},
\end{aligned}
```

```math
\begin{aligned}
\kappa_1^2&=u^2-\lambda_1^2,\quad \mathop{\mathrm{Re}}\kappa_1>0 \\
\kappa_2^2&=u^2-\lambda_2^2,\quad \mathop{\mathrm{Re}}\kappa_2>0 \\
Q&=\frac{\Gamma}{j\omega}I.
\end{aligned}
```

Here ``H_0^{(1)}`` is the first-kind order-zero Hankel function, ``\rho`` is distance to the wire centre, and ``\rho_0`` distance to its image. Fig. 1 and the explicitly printed direct-distance square root in (4·3) give the geometrical identification ``\rho=\sqrt{x^2+(y-h)^2}`` and ``\rho_0=\sqrt{x^2+(y+h)^2}``; this is a notation map of the shown distances, not a new self-radius regularization. The sign of ``y`` is retained in the earth exponential. Both scalar potentials have ``\varepsilon_1``, not ``\varepsilon_2``, in their printed prefactors.

The symbol list gives ``k^{(r)}>0,\ k^{(i)}<0`` for the complex material constant. It gives no separate Hankel-argument branch prescription for ``\lambda``. Retain that limitation alongside the explicit ``\mathop{\mathrm{Re}}\kappa_n>0`` Fourier conditions; do not replace the printed first-kind Hankel function with another kind.

Section 9, p. 729, explicitly defines the terminal quantity by:

```math
V=V_1(0,h-a).
```

The source describes an infinitely remote return reference and notes that a practical earth-surface reference differs because the earth potential is not spatially constant. No division by ``Q``, multiplication by ``j\omega``, or matrix inversion has been performed here to manufacture an unprinted line-admittance expression.

**Approximation.** The Hankel-plus-integral representation is not an analytical approximation within the stated thin-wire, homogeneous-medium, and principal-mode model. Its physical line-current approximation is explicit. The author's separate first approximation (5·9) uses small Hankel arguments and ``\Gamma\simeq jk_1``, leading to the §6 evaluators and §7 logarithmic reductions. Those equations evaluate or specialize the same potential representation and are not substituted here.

**Limitations.** The source does not provide a complete independent dispersion closure for ``\Gamma`` in this equation set, arbitrary earth permeability, a buried source, multiple-wire assembly, or a source-defined admittance matrix. The general ``\lambda`` branch and the relation of the stated complex-``k`` sign to lossless limiting cases remain unresolved. Original page verification does not establish earliest priority.

**Reference.** [Kikuchi1957](@cite), (3·5), (3·7), (4·4), (5·1)–(5·4), pp. 722–724, and terminal-reference discussion in §9, p. 729.

**Transcription source.** Original June 1957 publication images, not the PDF's damaged OCR text or Iwamoto's 1958 citation. All thirteen pages were inspected; the retained scalar equations, column-spanning bulk definition, signed earth exponential, prefactors and terminal definition were additionally inspected in enlarged crops. The last part of PDF page 13 is an unrelated conference announcement and supplies no formulation.

## Source transcription

The original symbols remain unchanged. Source labels use a centred dot; labels below identify the corresponding printed equation group without assigning new equation numbers.

### Source equation and boundary coupling — p. 722, (3·5) and (3·7)

```math
\begin{aligned}
(\nabla^2+k_1^2)V_1(\mathbf r)&=-\frac{1}{\varepsilon_1}Q\delta(\mathbf r-\mathbf r_0) \\
(\nabla^2+k_2^2)V_2(\mathbf r)&=0 \\
Q&=\frac{\Gamma}{j\omega}I.
\end{aligned}
```

```math
(V_1)_{y=+0}=(V_2)_{y=-0},
```

```math
k_1^2\left(\frac{\partial V_1}{\partial y}\right)_{y=+0}
-k_2^2\left(\frac{\partial V_2}{\partial y}\right)_{y=-0}
=j\omega(k_2^2-k_1^2)(A_{1y})_{y=+0}.
```

The derivative weights are transcribed in their printed order, not interchanged to fit an independently expected boundary condition.

### General potentials — pp. 723–724, scalar lines of (5·1)–(5·4)

The four scalar expressions are transcribed in full in **Expression**, in source order: air potential (5·1), earth potential (5·2), air auxiliary (5·3), earth auxiliary (5·4). Their display above places each helper next to its dependent potential for readability; their labels and original order are explicitly preserved here.

The vector-potential dependency of the boundary condition is printed in (5·1) and (5·3):

```math
A_{1y}(x,y)=\frac{\mu_1 I}{\pi}(Q_{1y}-jP_{1y}),
```

```math
Q_{1y}-jP_{1y}
=\Gamma\int_0^\infty
\frac{\sqrt{u^2-\lambda_2^2}-\sqrt{u^2-\lambda_1^2}}
{k_2^2\sqrt{u^2-\lambda_1^2}+k_1^2\sqrt{u^2-\lambda_2^2}}
e^{-(y+h)\sqrt{u^2-\lambda_1^2}}\cos xu\,du.
```

This auxiliary is included to make the cited boundary condition readable, not counted as an independent external-impedance formulation. In particular, ``A_{1y}`` is not ``Z_e``.

### Terminal reference — p. 729, §9

The author's prescription is ``V=V_1(0,h-a)``. The subsequent approximate logarithmic voltage (9·2) omits an earth-potential auxiliary under a stated dominance condition. That omission is not made in the general terminal prescription above. Neither a new full finite-radius integral nor a surface-averaged self coefficient is derived here.

## Notation map

| Source symbol | Display symbol | Physical meaning, units or convention |
| --- | --- | --- |
| ``V_1,V_2`` | Unchanged | Scalar potentials in air and earth, volts in rationalized MKS; omitted common ``e^{-\Gamma z+j\omega t}`` |
| ``V`` | Unchanged | Wire terminal potential ``V_1(0,h-a)`` relative to remote return; not potential of all earth-surface points |
| ``I`` | Unchanged | Total wire current, amperes, directed longitudinally as in Fig. 1 |
| ``Q`` | Unchanged | Wire charge per unit length, C/m; ``Q=\Gamma I/(j\omega)`` |
| ``Q'_n,P'_n`` | Unchanged | Dimensionless scalar-potential auxiliary parts; primes are labels, not derivatives; not the total line charge ``Q`` |
| ``Q_{1y},P_{1y}`` | Unchanged | Vector-potential auxiliary parts, dimensionless by (5·3); subscript ``y`` is a component label |
| ``A_{1y}`` | Unchanged | Air vector-potential component; Wb/m in the source's MKS formulation |
| ``k_n`` | Unchanged | Source intrinsic bulk constant, m⁻¹; region 1 air, region 2 earth |
| ``k_n^{(r)},k_n^{(i)}`` | Unchanged | Printed real/imaginary parts; source gives positive real and negative imaginary sign conventions |
| ``\Gamma`` | Unchanged | Longitudinal modal constant, m⁻¹; not set to zero in this record |
| ``\lambda_n,\kappa_n,u`` | Unchanged | Transverse argument constant, Fourier decay constant and integration variable, respectively, m⁻¹; ``\kappa_n`` branch as printed |
| ``a,h,x,y,z,\rho,\rho_0`` | Unchanged | Metres; radius, height, Cartesian coordinates, direct/image distances; earth ``y<0`` |
| ``\mathbf r,\mathbf r_0,\delta`` | Unchanged | Observation/source-position vectors and Dirac delta in (3·5); no change of source dimensionality introduced |
| ``\varepsilon_n,\sigma_n`` | Unchanged | Permittivity (F/m) and conductivity (S/m); no added dispersion or dielectric-loss model |
| ``\mu_0,\mu_1,\mu_2`` | Unchanged | Permeabilities (H/m); common ``\mu_0`` in bulk definition, explicit ``\mu_1=\mu_2`` |
| ``\omega,j,\pi,H_0^{(1)}`` | Unchanged | Angular frequency (s⁻¹), imaginary unit, circle constant, first-kind order-zero Hankel function |
| ``\gamma`` in the source symbol list | Not used in the retained equations | ``1.78107\ldots=e^c``, with Euler constant ``c``; this is not a bulk propagation constant |

No mathematical symbols were renamed. Units of the auxiliary parts follow their displayed integral measures and prefactors (equation-implied), rather than a separate source units table.

## Evidence and approximation sources

1. **Identity and source scope.** First-page byline and English footnote, printed p. 721; running headers and final article material, pp. 722–733. The existing bibliography's key and metadata are copied faithfully. The actual print byline is 菊地, not a silently substituted OCR spelling.
2. **Physical reduction.** Source pp. 721–722 and §11, p. 732, state the small-radius/exterior-line-current and principal-mode treatment. Arbitrary ``\Gamma`` evaluation or a dispersion root-selection algorithm is not supplied merely because ``\Gamma`` occurs symbolically.
3. **Material and branch evidence.** Column-spanning definition on p. 722 uses ``\mu_0`` under the square root and separately writes ``\mu_1=\mu_2``. The Fourier roots have ``\mathop{\mathrm{Re}}\kappa_n>0`` in (4·4). These facts are stronger than guesses from later LCM formulas and are not replaced by arbitrary layer permeability.
4. **General versus approximate formulas.** The unreduced (5·1)–(5·4) are distinct from (5·9), which introduces small-argument Hankel expressions, ``\lambda_1\simeq0``, ``\Gamma\simeq jk_1`` and ``\lambda_2^2\simeq k_2^2-k_1^2``. Section 6 then expands auxiliary integrals; §7 gives special observation geometry; §9 gives the further terminal-voltage simplification. Their approximations must not be assigned to the unreduced equation set.
5. **Actual output and reference.** The scalar-potential kernel is preserved with its published ``\Gamma I/(j\omega\varepsilon_1\pi)`` prefactor. The source defines ``Q`` separately and defines terminal voltage in §9. No new ``V/Q`` coefficient, ``Y`` conversion, vector-potential-to-impedance relation or multiwire matrix is inferred.
6. **Earlier publications.** References (6)–(9), p. 732, identify Kikuchi's JIEE 75 p. 1176 (1955), E.T.J. of Japan 2 p. 73 (1956), Bulletin of the Electrotechnical Laboratory, Japan 21 p. 49 (1957), and the first subcommittee radio-interference report (1956). The introduction says earlier transmission work and a summary preceded this complete treatment. Those citations do not establish that the retained 1957 formulas first appeared in an earlier source.
7. **Source scope.** Sections 4–5 contain Fourier construction and field/potential equations; §§6–7 contain auxiliary expansions and observation limits; §8 studies circuit/Poynting power; §9 studies terminal voltage versus measured ground fields; and §10 gives examples. The paper does not supply a full ``Z_e`` or multiwire ``Y`` kernel.

## Limitations and discrepancies

- **Unresolved branch/closure.** The retained potential is symbolic in ``\Gamma``; the explicit ``\kappa_n`` decay conditions do not specify a separate ``\lambda_n`` Hankel branch or determine ``\Gamma`` numerically. Do not import a later propagation prescription or select a different Hankel kind.
- **Air constitutive restriction.** The common bulk definition contains ``\sigma_n``, while the scalar source/prefactor uses ``\varepsilon_1`` and ``j\omega`` without an air-conductivity term. The generality of a nonzero ``\sigma_1`` extension is unresolved; none is supplied here.
- **Potential is not admittance.** A below-surface observation in (5·2)/(5·4) does not supply a buried source, mutual cable coefficient, or mixed-arrangement line-admittance matrix. Mixed coverage is supplied by separate Pawlik (2018) and Martins-Britto et al. (2024) records; it is not inferred from Kikuchi's observation field.
- **Reference potential.** Section 9 explicitly notes nonconstant earth potential. Do not replace the remote return reference by a local equipotential-earth assumption.
- **Printed boundary ordering.** In (3·7), the air derivative is weighted by ``k_1^2`` and the earth derivative by ``k_2^2``. The retained final potential denominator has the distinct printed crossed material weights. No algebraic reconciliation or corrected boundary condition is claimed.
- **OCR navigation defects.** The PDF text layer omits many equations and flattens vector/operator typography. It is not used to choose the mathematical symbols in this record. No original conversion was edited.
- **Related formulas.** Source (5·9), the auxiliary expansions (6·1)–(6·7), logarithmic observation reduction (7·5), and terminal relations (9·1)–(9·2) evaluate or specialize the same scalar-potential construction. They are not independent physical formulations.
- **Priority.** The page-image evidence does not establish that the 1957 publication was the first source of the formula.
