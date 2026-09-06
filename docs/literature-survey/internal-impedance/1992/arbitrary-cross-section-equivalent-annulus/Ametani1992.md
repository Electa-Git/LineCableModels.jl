# Ametani–Fuse arbitrary-cross-section internal impedance approximation

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Isolated homogeneous conductor of arbitrary cross-section with area ``S`` and perimeter ``\ell``; this includes solid and annular sectors. |
| Calculated quantities | Approximate frequency-dependent p.u.l. internal impedance and the radii of an impedance-equivalent circular annulus |
| Earth structure | Not applicable. |
| Model and approximation | ``Z_i=R_{dc}\sqrt{1+j\omega\mu_cS/(R_{dc}\ell^2)}`` interpolates between ``R_{dc}=\rho_c/S`` and the high-frequency surface-layer limit ``\sqrt{j\omega\mu_c\rho_c}/\ell``. The equivalent annulus has ``r_o=\ell/(2\pi)`` and ``r_i=\sqrt{r_o^2-S/\pi}``. The approximation retains skin effect through area and perimeter but does not resolve corner current crowding or proximity effect. |
| Main source | Akihiro Ametani and Ikuko Fuse (1992 English translation; 1991 Japanese original) |
| Citation key(s) | Primary English publication: `:Ametani1992`; Japanese original: `:Ametani1991`; later equation witness: `:Ametani2021` |
| Evidence status | Japanese original equations (1)–(6) and (14)–(15), sector-cable examples, and the later book's equation (A1.59) checked against page images; English publication metadata and abstract checked |

**Description.** The approximation uses the conductor area and perimeter to join the direct-current and surface-layer limits. It applies to arbitrary homogeneous cross-sections, including sector-shaped conductors.

**Expression.** The approximation is

```math
\begin{aligned}
Z_i&=R_{dc}
\sqrt{1+\frac{j\omega\mu_c S}{R_{dc}\ell^2}} \\
R_{dc}&=\frac{\rho_c}{S}.
\end{aligned}\qquad\text{(1)}
```

Its high-frequency limit is

```math
Z_i\sim\frac{\sqrt{j\omega\mu_c\rho_c}}{\ell}.
\qquad\text{(2)}
```

The cross-section can instead be represented by an impedance-equivalent circular annulus with

```math
\begin{aligned}
r_o&=\frac{\ell}{2\pi} \\
r_i&=\sqrt{r_o^2-\frac{S}{\pi}}.
\end{aligned}\qquad\text{(3)}
```

For a solid sector with radius ``r`` and angle ``\theta`` in radians,

```math
\begin{aligned}
S&=\frac{\theta r^2}{2} \\
\ell&=r(2+\theta).
\end{aligned}\qquad\text{(4)}
```

For an annular sector with inner radius ``r_i^{(s)}`` and outer radius ``r_o^{(s)}``, use its physical area and full wetted perimeter in (1)–(3).

**Implementation.** Calculate ``S`` and ``\ell`` from the physical cross-section. Equation (1) is the direct scalar evaluator. Equations (3) define an annulus that can be passed to an existing circular-tube evaluator.

**Limitations.** The formula does not resolve corners, localized current crowding, proximity effect, composite conductors, or magnetic anisotropy. It gives an isolated-conductor approximation.

**Reference.** [Ametani1992](@cite), with the Japanese original [Ametani1991](@cite) and later equation witness [Ametani2021](@cite).

## Notation map

| Symbol | Meaning |
| --- | --- |
| ``S`` | Physical conductor cross-sectional area |
| ``\ell`` | Physical conductor perimeter |
| ``\rho_c,\mu_c`` | Conductor resistivity and permeability |
| ``R_{dc}`` | Direct-current resistance per unit length |
| ``r_i,r_o`` | Inner and outer radii of the impedance-equivalent annulus |

## Evidence and approximation sources

- The Japanese original supplies the defining arbitrary-section relations and sector examples in equations (1)–(6) and (14)–(15).
- The English publication is the primary index reference.
- The later cable-modeling book repeats the approximation as equation (A1.59).

