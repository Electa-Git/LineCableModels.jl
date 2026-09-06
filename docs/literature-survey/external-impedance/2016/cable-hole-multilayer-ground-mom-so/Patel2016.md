# Patel–Triverio multilayer-ground MoM–SO cable impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Round conductors and circular cable hole; sheaths/armor may be hollow conductors or explicit strands. |
| Calculated quantities | Complete p.u.l. cable series matrix with arbitrary horizontal ground layers, cable hole, skin and proximity effects |
| Earth structure | Arbitrary number of flat horizontal layers; top and bottom are semi-infinite. |
| Model and approximation | Layer stacks are represented by exact equivalent input impedances in a spectral transmission-line circuit, then inverse-transformed and MoM-discretized. Numerical Fourier and boundary truncations remain. |
| Main source | U. R. Patel and P. Triverio (2016) |
| Citation key(s) | `:Patel2016` |
| Evidence status | Original publication page images checked |

**Description.** Extension of the cable-hole MoM–SO operator to an arbitrary stack of horizontal conductive/dielectric/magnetic layers through an equivalent spectral transmission-line Green function.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Longitudinally invariant cable fields; telegrapher reduction. | Stated — §II. |
| Air propagation constant ``γ_air`` | Air may be any layer and uses ``k_l=\sqrt{\omega\mu_l(\omega\varepsilon_l-j\sigma_l)}`` with ``\sigma_l=0``. | Stated — §II, (9)–(13). |
| Earth propagation constant ``γ_earth`` | Each layer has its own ``k_l`` and vertical spectral root ``\gamma_l=\sqrt{\beta_x^2-k_l^2}``. | Stated — §§II, IV. |
| Earth permittivity and displacement current | Retained independently in every layer. | Stated — Fig. 1 and (9)–(13). |
| Range of validity | Fourier/MoM truncations; examples span 1 Hz–1 MHz and include air–sea–seabed, but no universal bound is claimed. | Stated — §VI. |
| Earth permeability ``μ_earth`` | Arbitrary layer ``\mu_l`` retained. | Stated — Fig. 1 and (9)–(13). |
| Arrangement | Arbitrary parallel solid/hollow round conductors in a circular hole; multiple systems are stated as supported. | Stated — §II. |
| Earth structure | Arbitrary number of flat horizontal layers; top and bottom are semi-infinite. | Stated — §II and Fig. 1. |
| Conductor and insulation geometry | Round conductors and circular cable hole; sheaths/armor may be hollow conductors or explicit strands. | Stated — §II. |
| Constitutive and field assumptions | Linear isotropic layers, 2-D harmonic fields, equivalence theorem, Fourier/MoM boundary discretization. | Stated — §§III–V. |
| Conventions | ``k_l=\sqrt{\omega\mu_l(\omega\varepsilon_l-j\sigma_l)}``; Fourier kernel ``e^{-j\beta_xx}``; p.u.l. ``R+j\omega L``. | Stated — (9)–(17), (20)–(22). |

**Expression.** For source and observation in layer ``s``, the spectral Green function is

```math
\widetilde G_g(\beta_x,y)=\frac{Z_sI_s}{2}e^{-|y-y'|\gamma_s}
+\frac{Z_sI_s}{2}\frac{\Gamma_L e^{(2y_s-y-y')\gamma_s}+\Gamma_R e^{(-2y_{s-1}+y+y')\gamma_s}
+\Gamma_R\Gamma_L e^{(-2y_{s-1}+2y_s+y'-y)\gamma_s}}
{1-\Gamma_R\Gamma_L e^{-2(y_{s-1}-y_s)\gamma_s}},
\qquad\text{(14)}
```

with ``Z_l=(\beta_x^2-k_l^2)^{-1/2}``, ``\gamma_l=\sqrt{\beta_x^2-k_l^2}``, and the inverse transform

```math
G_g(x,y)=\frac1{2\pi}\int_{-\infty}^{\infty}\widetilde G_g(\beta_x,y)e^{-j\beta_xx}\,d\beta_x.
\qquad\text{(17)}
```

The final matrices are

```math
\mathbf R(\omega)=\Re\{[\mathbf U^T(\mathbf1-j\omega\mathbf Y_s\mathbf\Psi)^{-1}\mathbf Y_s\mathbf U]^{-1}\},
\quad
\mathbf L(\omega)=\omega^{-1}\Im\{[\mathbf U^T(\mathbf1-j\omega\mathbf Y_s\mathbf\Psi)^{-1}\mathbf Y_s\mathbf U]^{-1}\},
\qquad\text{(21,22)}
```

where ``\mathbf\Psi`` is printed in (23).

**Approximation.** Layer stacks are represented by exact equivalent input impedances in a spectral transmission-line circuit, then inverse-transformed and MoM-discretized. Numerical Fourier and boundary truncations remain.

**Limitations.** Horizontally stratified, laterally infinite media and longitudinally invariant circular cable-hole geometry. This is a numerical operator formulation, not a scalar self/mutual closed form.

**Reference.** [Patel2016](@cite), equations (1)–(23), especially (14)–(23), printed pp. 1235–1236.

**Transcription source.** Original IEEE page images. Equivalent-line definitions, four numerator exponentials, reflection-product denominator, inverse transform and final matrix ordering were visually checked.

## Source transcription

The layer-stack input impedances ``Z_{eq,s+1}`` and ``Z_{eq,s-1}`` generate ``\Gamma_L,\Gamma_R`` in (15)–(16). Their textbook input-impedance recursion is cited rather than printed; this dependency remains explicit.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\widetilde G_g,G_g`` | unchanged | spectral/spatial multilayer magnetic Green function | source normalized |
| ``Z_l,\gamma_l`` | unchanged | equivalent spectral line impedance/vertical constant | source convention |
| ``\Gamma_L,\Gamma_R`` | unchanged | lower/upper stack reflection coefficients | dimensionless |
| ``\mathbf\Psi`` | unchanged | full cable-hole/exterior potential operator | source normalized |

No notation was renamed.

## Evidence and approximation sources

The source identifies the arbitrary multilayer Green function as its extension of the 2015 homogeneous-ground method. Conductor surface operators are inherited and not re-credited.

## Limitations and discrepancies

- Equation (14) is coordinate/order sensitive; this record keeps the printed layer-boundary indices and does not convert to another recursion convention.

