# Brandão Faria inhomogeneous Euler–Cauchy tubular-conductor impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Circular tube with inner radius ``r_1`` and outer radius ``r_2``; the solid limit is ``r_1\to0``; no insulation layer. |
| Calculated quantities | Frequency-dependent p.u.l. internal impedance of a radially inhomogeneous tubular conductor, plus its solid-cylinder limit |
| Earth structure | Not applicable. |
| Model and approximation | Not an analytical approximation after imposing the very-good-conductor diffusion model and the specific constitutive power laws. The Euler–Cauchy reduction is exact for (11)–(12); arbitrary radial profiles require numerical solution. |
| Main source | José António Brandão Faria (2011) |
| Citation key(s) | `:BrandaoFaria2011` |
| Evidence status | Original publication page images checked |

**Description.** Internal per-unit-length impedance of a circular tubular conductor whose scalar permeability and conductivity follow the coupled radial power laws that reduce the diffusion equation to an Euler–Cauchy equation.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Not applicable; the conductor fields are axially invariant and driven by total axial current ``\bar I``. | Stated — §2 and Fig. 1. |
| Air propagation constant ``γ_air`` | Not applicable. | The model is confined to the conductor region. |
| Earth propagation constant ``γ_earth`` | Not applicable. | No earth medium occurs. |
| Earth permittivity and displacement current | Not applicable. | No earth medium occurs. |
| Range of validity | Very-good-conductor approximation ``\sigma\gg\omega\epsilon``; exact closed form only for the paired power-law profiles in (11)–(12), ``r_1\le r\le r_2``. | Stated — §2, p. 90, and (8)–(13), pp. 92–93. |
| Earth permeability ``μ_earth`` | Not applicable. | No earth medium occurs. |
| Arrangement | Not applicable to earth placement; a single isolated current-carrying tubular conductor is treated. | Stated — §2 and Fig. 1. |
| Earth structure | Not applicable. | No earth medium occurs. |
| Conductor and insulation geometry | Circular tube with inner radius ``r_1`` and outer radius ``r_2``; the solid limit is ``r_1\to0``; no insulation layer. | Stated — §3 and following (19). |
| Constitutive and field assumptions | Linear isotropic conductor with radial scalar ``\mu(r)`` and ``\sigma(r)``; displacement current neglected, azimuthal magnetic field and axial current/electric field, Coulomb gauge. | Stated — (1)–(5), pp. 90–92. |
| Conventions | Time dependence ``e^{j\omega t}``; overbars denote complex amplitudes; impedance is ``\bar E(r_2)/\bar I`` in ``\Omega/\mathrm m``. | Stated — §2 and (19). |

**Expression.** The source defines the admissible Euler–Cauchy material class by

```math
\mu(r)=\mu_2\left(\frac r{r_2}\right)^p,
\qquad
\sigma(r)=\sigma_2\left(\frac{r_2}{r}\right)^{2+p},
\qquad r_1\le r\le r_2,
\qquad\text{(11,12)}
```

and

```math
m_{1,2}=\frac p2\pm
\sqrt{\left(\frac p2\right)^2-(\bar k_2r_2)^2},
\qquad
\bar k_2=\sqrt{-j\omega\mu_2\sigma_2},
\qquad\text{(6,13)}
```

with ``m_1`` selected in the first quadrant and ``m_2`` in the third. The resulting p.u.l. impedance is

```math
\bar Z=R+jX=\frac{\bar E(r_2)}{\bar I}
=\frac{m_2(r_1/r_2)^{m_2}-m_1(r_1/r_2)^{m_1}}
{2\pi\sigma_2r_2^2\left((r_1/r_2)^{m_1}-(r_1/r_2)^{m_2}\right)}.
\qquad\text{(19)}
```

For ``r_1\to0``, the source gives

```math
\bar Z=\frac{j\omega\mu_2}{2\pi m_1}.
```

**Approximation.** Not an analytical approximation after imposing the very-good-conductor diffusion model and the specific constitutive power laws. The Euler–Cauchy reduction is exact for (11)–(12); arbitrary radial profiles generally require numerical solution.

**Limitations.** This is not a general inhomogeneous-conductor solution: ``\mu(r)`` and ``\sigma(r)`` are coupled so that ``r^2\mu(r)\sigma(r)`` is constant. The isolated-conductor boundary condition ``\bar H(r_1)=0`` excludes current in the hollow interior and external proximity fields.

**Reference.** [BrandaoFaria2011](@cite).  Brandão Faria, “Skin Effect in Inhomogeneous Euler–Cauchy Tubular Conductors,” *Progress In Electromagnetics Research M* 18 (2011), equations (7)–(19), pp. 92–94.

**Transcription source.** Original publication page images. The material exponents, radical sign, root-quadrant selection, order and signs of the numerator/denominator in (19), and solid limit were visually checked on printed pp. 93–94.

## Source transcription

Before (19), the source prints

```math
\bar H(r)=\frac{\bar I}{2\pi r}
\frac{(r_1/r)^{m_1}-(r_1/r)^{m_2}}
{(r_1/r_2)^{m_1}-(r_1/r_2)^{m_2}},
\qquad\text{(17)}
```

```math
\bar J(r)=\frac{\bar I}{2\pi r^2}
\frac{m_2(r_1/r)^{m_2}-m_1(r_1/r)^{m_1}}
{(r_1/r_2)^{m_1}-(r_1/r_2)^{m_2}},
\qquad\text{(18)}
```

under ``\bar H(r_1)=0`` and ``\bar H(r_2)=\bar I/(2\pi r_2)``. These dependencies make the impedance boundary ratio evaluable.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\bar Z`` | unchanged | p.u.l. internal impedance | ``\Omega/\mathrm m`` |
| ``r_1,r_2`` | unchanged | inner and outer radii | m |
| ``\mu_2,\sigma_2`` | unchanged | outer-boundary material values | H/m; S/m |
| ``p`` | unchanged | permeability power-law exponent | dimensionless real number |
| ``m_1,m_2`` | unchanged | characteristic exponents | dimensionless complex roots |
| ``\bar k_2`` | unchanged | outer-boundary diffusion wavenumber | ``\mathrm m^{-1}`` |

No source variable was renamed.

## Evidence and approximation sources

Equations (4), (7), and (8) show the operation: variable-coefficient magnetic-vector-potential diffusion is reduced to a constant-coefficient Euler–Cauchy equation. Equations (11)–(13) are the necessary constitutive and branch definitions, and (15)–(19) apply the isolated-tube boundary conditions.

## Limitations and discrepancies

- Equation (13) prints set-membership notation for the signs and quadrants of the roots; the record states the accompanying prose instead of trying to normalize that typography.
- The paper's statement that conductor displacement current is negligible “up to the optical range” is an author assertion under ``\sigma\gg\omega\epsilon``, not converted here into a universal frequency bound.
