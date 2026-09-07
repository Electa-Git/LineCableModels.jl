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

## Numerical interpretation and backend scope

`InternalImpedance.Formula(:Ametani1992)` accepts the physical area and
perimeter through its `Val(:section)` call. The state also supplies the
equivalent-annulus radii in (3). The numerical expression is the
author-defined combination of the dc and high-frequency limits,
``\sqrt{R_{dc}^2+Z_{hf}^2}``, with
``Z_{hf}=\sqrt{j\omega\mu_c\rho_c}/\ell``. This avoids the
dimensionally inconsistent expanded expression in the later book.

The adapter preserves the original homogeneous material and the actual
area and perimeter of each supported isolated section. Sector contours,
including fillets, are measured before circular reduction. Temperature
correction applies to the physical resistivity. For a complete circular
annulus carrying isolated-conductor current, the relevant perimeter is
the outer circumference; an open sector uses its full contour.

This route supplies one scalar internal term, not inner-surface or
transfer coefficients for nested conductive terminals. Strand formations,
composite materials, and concentric terminal assemblies are rejected.
The selected external and insulation formulas still use the backend's
equivalent geometry; this internal approximation does not resolve their
noncircular fields.

Tests check the source interpolation, equivalent-annulus area and outer
perimeter, zero-frequency and negative-frequency limits, and complete
matrix assembly with solid, annular, sharp-sector, and filleted-sector
conductors. The sector test explicitly distinguishes its physical
perimeter from the circumference of an equal-area circle.

## Evidence and approximation sources

- The Japanese original supplies the defining arbitrary-section relations and sector examples in equations (1)–(6) and (14)–(15).
- The English publication is the primary index reference.
- The later cable-modeling book repeats the approximation as equation (A1.59).
