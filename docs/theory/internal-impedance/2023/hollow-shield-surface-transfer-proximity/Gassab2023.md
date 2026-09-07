# Gassab et al. hollow-shield surface and transfer impedances with proximity

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Circular solid or annular conductors; hollow shield has inner radius ``r_1`` and outer radius ``R_1``. |
| Calculated quantities | Inner/outer p.u.l. surface impedances and through-wall transfer impedance of a hollow circular shield with asymmetric proximity fields |
| Earth structure | No earth model; return is another conductor/shield. |
| Model and approximation | Exact coefficients come from the complete Bessel boundary solution. Equations (36a)–(36c) replace the exterior magnetic distributions by the closed proximity factors ``\Pi`` and ``\Upsilon_1``; the source states this is accurate when conductor thickness/radius is large compared with skin depth. |
| Main source | Oussama Gassab, Hao Xie, Yanning Chen, Ding-E Wen, Sichao Du, Kang Luo, Fangmin He, Jin Meng, Dongyan Zhao, Duo Xiao, and Wen-Yan Yin (2023) |
| Citation key(s) | `:Gassab2023` |
| Evidence status | Original publication page images checked |

**Description.** Semi-analytical harmonic solution for asymmetric skin and proximity fields in solid and hollow circular conductors, reduced to exact series and high-frequency closed approximations for a shield's inner/outer surface impedances and transfer impedance.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Longitudinally invariant conductor fields; no axial propagation term appears in (4). | Stated — §II. |
| Air propagation constant ``γ_air`` | Not an earth-return formulation; exterior fields use two-dimensional magnetostatic/equipotential constructions. | Stated — §§II–IV. |
| Earth propagation constant ``γ_earth`` | Not applicable. | Scope — §I. |
| Earth permittivity and displacement current | Conductor displacement current is neglected in ``\nabla\times H=J=\sigma_cE``. | Stated — (1b)–(1c). |
| Range of validity | Exact series applies within the stated 2-D conductor model; simplified formulas target the high-frequency region where skin depth is much smaller than conductor radius. | Stated — abstract, §§IV, VII-C. |
| Earth permeability ``μ_earth`` | Not applicable; conductor permeability ``\mu_c`` is retained. | Stated — (1a), (4). |
| Arrangement | Adjacent solid/hollow conductors of different radii; noncoaxial and twinax shield configurations. | Stated — Figs. 2, 3, 12. |
| Earth structure | No earth model; return is another conductor/shield. | Stated — scope. |
| Conductor and insulation geometry | Circular solid or annular conductors; hollow shield has inner radius ``r_1`` and outer radius ``R_1``. | Stated — Fig. 1. |
| Constitutive and field assumptions | Linear isotropic good conductors, harmonic 2-D fields; exact Bessel harmonic expansion, plus perfect-conductor exterior-potential approximation at high frequency. | Stated — §§II–IV. |
| Conventions | ``e^{j\omega t}``; currents at the inner/outer shield surfaces have opposite reference directions. | Stated — (1a), text after (34). |

**Expression.** The exact p.u.l. inner surface, outer surface and transfer impedances are

```math
Z_{s1}^{\mathrm{int}}=
\frac{\gamma_cG_0(R_1,r_1)}{2\pi r_1\sigma_c}
\left[1+\sum_{m=1}^{\infty}\frac{G_m(R_1,r_1)}{G_0(R_1,r_1)}
\left((a_m^{\mathrm{int}})^2+(b_m^{\mathrm{int}})^2\right)\right],
\qquad\text{(35a)}
```

```math
Z_{s1}^{\mathrm{ext}}=
\frac{\gamma_cG_0(r_1,R_1)}{2\pi R_1\sigma_c}
\left[1+\sum_{m=1}^{\infty}\frac{G_m(r_1,R_1)}{G_0(r_1,R_1)}(a_m^{\mathrm{ext}})^2\right],
\qquad\text{(35b)}
```

```math
Z_{t1}=\frac{1}{2\pi R_1r_1\sigma_c\Delta_0}
\left[1+\sum_{m=1}^{\infty}\frac{G_m(r_1,r_1)}{G_0(r_1,r_1)}a_m^{\mathrm{ext}}a_m^{\mathrm{int}}\right].
\qquad\text{(35c)}
```

The compact high-frequency transfer approximation is

```math
Z_{t1}=\frac{1}{4\pi^2R_1r_1\sigma_c\Delta_0}
\int_0^{2\pi}\Upsilon_1(\phi,\beta_1,\beta_2)\Pi(\phi,\alpha_1)\,d\phi.
\qquad\text{(36c)}
```

**Approximation.** Exact coefficients come from the complete Bessel boundary solution. Equations (36a)–(36c) replace the exterior magnetic distributions by the closed proximity factors ``\Pi`` and ``\Upsilon_1``; the source states this is accurate when conductor thickness/radius is large compared with skin depth.

**Limitations.** Circular, parallel and longitudinally invariant conductors only; no braid apertures, dielectric-transfer admittance, ground return, finite length, or general cable nodal assembly. The simplified formulas degrade in the low-frequency transition region.

**Reference.** [Gassab2023](@cite).  Gassab et al., *IEEE Transactions on Electromagnetic Compatibility* 65(6), 2023, DOI `10.1109/TEMC.2023.3298809`, equations (1)–(42), especially (35)–(40), printed pp. 1639–1642.

**Transcription source.** Original IEEE page images. Radii ordering, Bessel-series factors, current-direction sign explanation, summation terms and the product in (36c) were visually verified.

## Source transcription

The paper defines ``\gamma_c=\sqrt{j\omega\mu_c\sigma_c}`` in (4). Functions ``G_m`` and ``\Delta_m`` and harmonic coefficients ``a_m,b_m`` are printed in the field development and appendices; this record does not replace those dependencies with a guessed coaxial form.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``r_1,R_1`` | unchanged | inner and outer shield radii | m |
| ``Z_{s1}^{int},Z_{s1}^{ext}`` | unchanged | p.u.l. surface impedances | ``\Omega/\mathrm m`` |
| ``Z_{t1}`` | unchanged | through-wall p.u.l. transfer impedance | ``\Omega/\mathrm m`` |
| ``\Pi,\Upsilon_1`` | unchanged | approximate exterior proximity-field factors | dimensionless |

No notation was renamed.

## Evidence and approximation sources

Equations (35) are explicitly introduced as exact formulas and (36) as approximate formulas. The transfer result is used in the later twinax coupling equations but is itself a physical shield-impedance formula.

## Limitations and discrepancies

- The source's rule-of-thumb transition-frequency variable ``x_0`` in (42) is not a universal validity bound.
- One validation paragraph repeats ``r_1,R_1,t_1`` for conductor 2 where its context indicates the second conductor; no symbol repair is made here.
