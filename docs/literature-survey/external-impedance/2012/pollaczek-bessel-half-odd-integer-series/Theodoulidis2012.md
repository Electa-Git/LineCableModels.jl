# Theodoulidis first exact series for Pollaczek's integral

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel buried conductors at depths ``h_1,h_2`` and horizontal separation ``x``; the source uses a filamentary mutual geometry and the stated radius substitution for self. No insulation region appears. |
| Calculated quantities | Exact modified-Bessel series for the Pollaczek interface integral, used in buried-conductor self and mutual earth-return impedance |
| Earth structure | Homogeneous earth half-space below air. |
| Model and approximation | Not an analytical approximation of ``J_{\mathrm{Pollaczek}}``: the infinite series is derived exactly from (2). Finite truncation is a numerical approximation; the author's rule of thumb is roughly ``20x/H`` terms for ``x>H``. The physical parent remains the TEM Pollaczek model, with optional Sunde constitutive replacement. |
| Main source | Theodoros Theodoulidis (2012) |
| Citation key(s) | `:Theodoulidis2012` |
| Evidence status | Original publication page images checked |

**Description.** Exact convergent modified-Bessel series for the interface integral in the homogeneous-earth Pollaczek impedance of two parallel buried conductors. The parent impedance supplies mutual and source-prescribed self cases under the transmission-line form of the Pollaczek model.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | No independent ``\Gamma`` occurs. The author confines the analysis to the TEM/transmission-line limit of Pollaczek's buried-conductor formula; no full-wave modal closure is retained. | Stated — printed p. 807, paragraph preceding §II; Equation-implied — (1)–(2), pp. 807–808. |
| Air propagation constant ``γ_air`` | Not retained as an independent quantity in the parent formula. Air enters through the planar interface and image distance ``R``. | Equation-implied — (1)–(2), printed p. 807. |
| Earth propagation constant ``γ_earth`` | Source symbol ``k``. Classical formula: ``k=\sqrt{j\omega\mu_0\sigma}=(1+j)/\delta``. The source separately permits Sunde's higher-frequency replacement ``k=\sqrt{j\omega\mu_0\sigma+\omega^2\mu_0\epsilon_0\epsilon_r}``. | Stated — below (2), pp. 807–808. |
| Earth permittivity and displacement current | Omitted in the classical Pollaczek ``k``; optionally retained through the separately identified Sunde replacement containing ``\epsilon_0\epsilon_r``. | Stated — printed p. 808, paragraph above (3). |
| Range of validity | The infinite series converges for every ``x/H`` because ``|Z|>|z|``. A practical rule is about ``20x/H`` terms for ``x>H``; only a few terms are normally needed for ``x<H``. This mathematical convergence does not remove the parent impedance's high-frequency negative-resistance limitation. | Stated — printed p. 809 and validity discussion on pp. 812–813. |
| Earth permeability ``μ_earth`` | Fixed to ``\mu_0``. | Stated — §II and (1)–(2), printed p. 807. |
| Arrangement | Two underground parallel conductors; mutual impedance is printed. The self prescription is ``x`` equal to conductor radius with ``h_1=h_2``. A mixed overhead/underground extension is mentioned but not printed and is not claimed here. | Stated — Fig. 1 and text below (2), pp. 807–808. |
| Earth structure | Homogeneous earth half-space below air. | Stated — Fig. 1 and §II, printed p. 807. |
| Conductor and insulation geometry | Parallel buried conductors at depths ``h_1,h_2`` and horizontal separation ``x``; the source uses a filamentary mutual geometry and the stated radius substitution for self. No insulation region appears. | Stated — Fig. 1 and text below (2), pp. 807–808. |
| Constitutive and field assumptions | Linear isotropic homogeneous earth; TEM/transmission-line Pollaczek parent; conductor skin/proximity and finite insulation are outside this term. The exactness is mathematical exactness for the selected integral, not a full-wave correction. | Stated — printed pp. 807 and 812–813. |
| Conventions | ``j`` is the imaginary unit, ``\omega=2\pi f``, depths are positive downward in Fig. 1, and the displayed impedance is per unit length in SI. The chosen classical root has ``k=(1+j)/\delta``. | Stated — (1)–(2) and Fig. 1, printed p. 807; per-length use stated in the cable-constants discussion. |

**Expression.** First exact series for the auxiliary ``I_{\mathrm{Pollaczek}}`` in the parent buried-conductor impedance, equations (1)–(4), printed pp. 807–808.

```math
Z(j\omega)=\frac{j\omega\mu_0}{2\pi}
\left[K_0(kr)-K_0(kR)+2J_{\mathrm{Pollaczek}}\right],
\tag{1}
```

```math
J_{\mathrm{Pollaczek}}=\int_0^\infty
\frac{\exp(-Hu)}{\lambda+u}\cos(\lambda x)\,d\lambda,
\qquad
u=\sqrt{\lambda^2+k^2},
\tag{2}
```

```math
J_{\mathrm{Pollaczek}}
=\left(\frac{H}{R}\right)^2K_0(kR)
+\frac{1}{kR}\left[2\left(\frac{H}{R}\right)^2-1\right]K_1(kR)
-\frac{1}{k^2}I_{\mathrm{Pollaczek}},
\tag{3}
```

```math
I_{\mathrm{Pollaczek}}
=\frac{2}{xR}\sum_{n=0}^{\infty}(-1)^n(2n+1)^2
\left[
z I_{n+3/2}(z)K_{n+1/2}(Z)
+Z I_{n+1/2}(z)K_{n+3/2}(Z)
\right].
\tag{4}
```

The complete definitions are

```math
r=\sqrt{x^2+(h_1-h_2)^2},\qquad
R=\sqrt{x^2+H^2},\qquad H=h_1+h_2,
```

```math
k=\sqrt{j\omega\mu_0\sigma}=\frac{1+j}{\delta},
\qquad z=\frac{k(R-H)}{2},\qquad Z=\frac{k(R+H)}{2}.
```

``K_0,K_1`` in (1) and (3) and ``I_\nu,K_\nu`` in (4) are modified Bessel functions. ``Z`` used as the capital auxiliary argument in (4) is distinct from the impedance ``Z(j\omega)``. The source supplies the self formula by setting ``h_1=h_2`` and using the conductor radius for ``x``.

**Approximation.** Not an analytical approximation of ``J_{\mathrm{Pollaczek}}``: the infinite series is derived exactly from (2). Finite truncation is a numerical approximation; the author's rule of thumb is roughly ``20x/H`` terms for ``x>H``. The physical parent remains the TEM Pollaczek model, with optional Sunde constitutive replacement.

**Limitations.** The source warns that exact integration does not repair high-frequency failure of the parent Pollaczek impedance: the mutual resistance can become negative. For close conductors the direct ``K_0(kr)`` term gives a cited first-zero criterion near ``r/\delta\simeq2.8``; for interface-dominated coupling the cited estimate is ``f\simeq0.63\times10^6/(\sigma H^2)``. No mixed-case series is printed. Large orders or arguments can overflow without scaled Bessel evaluation; that numerical caution is not a change to (4).

**Reference.** [Theodoulidis2012](@cite), equations (1)–(4), printed pp. 807–809 (PDF pages 2–4), with validity discussion on pp. 812–813.

**Transcription source.** Original IEEE publication. Equations (1)–(4), every summation bound, half-odd-integer Bessel order, sign, factor, geometry definition and the ``|Z|>|z|`` convergence statement were checked against the PDF page images. The Markdown conversion was used only for navigation because it corrupts equation (4).

## Source transcription

The source order is parent impedance (1), integral (2), decomposition (3), and three alternative series (4)–(6). This record retains alternative (4). It does not merge it with alternatives (5) or (6), which are separate corpus records.

The derivation introduces

```math
I_H=\frac{2}{x}\sum_{n=0}^{\infty}(-1)^n(2n+1)^2
I_{n+1/2}(z)K_{n+1/2}(Z),
```

and obtains (4) by differentiating with respect to ``H``. This dependency records the source operation; it is not counted as another physical formula.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z(j\omega)`` | unchanged | Earth-return mutual or prescribed self impedance per unit length | ``\Omega/\mathrm m`` |
| ``J_{\mathrm{Pollaczek}}`` | unchanged | Pollaczek interface integral | ``\mathrm m^{-1}`` |
| ``I_{\mathrm{Pollaczek}}`` | unchanged | Auxiliary used in decomposition (3) | ``\mathrm m^{-2}`` |
| ``x,h_1,h_2,H,r,R`` | unchanged | Horizontal separation, depths, depth sum, direct and image distances | metres; depths positive downward |
| ``\lambda,u,k`` | unchanged | Spectral variable, vertical spectral root and earth bulk constant | ``\mathrm m^{-1}`` |
| ``z,Z`` | unchanged | Dimensionless Bessel arguments | capital ``Z`` is not impedance here |
| ``\sigma,\epsilon_0\epsilon_r,\mu_0`` | unchanged | Earth conductivity, optional permittivity and fixed permeability | SI |
| ``I_\nu,K_\nu`` | unchanged | Modified Bessel functions | order ``\nu`` |
| ``\delta`` | unchanged | Skin depth defined through ``k=(1+j)/\delta`` | metres |

No notation was renamed.

## Evidence and approximation sources

- Geometry, TEM restriction and parent equations: printed p. 807.
- Sunde full-``k`` option and self substitution: printed p. 808.
- Series (4) and definitions ``z,Z``: printed p. 808.
- Derivative construction, convergence proof and term-count guidance: printed p. 809.
- Negative-resistance diagnosis and critical-frequency estimates: printed pp. 812–813.
- The paper calls all three infinite series exact solutions. This label is retained only for mathematical evaluation of the selected parent integral.

## Limitations and discrepancies

- The source's mathematical ``Z`` argument and electrical ``Z(j\omega)`` share a letter. The record keeps both and states the distinction instead of renaming one.
- Exact series evaluation does not validate the physical Pollaczek formula outside its stated transmission-line regime.
