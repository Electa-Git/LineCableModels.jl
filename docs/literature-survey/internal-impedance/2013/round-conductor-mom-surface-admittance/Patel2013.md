# Patel–Gustavsen–Triverio round-conductor MoM–SO series impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Solid circular conductors of radii ``a_p``; insulation and hollow walls are not represented in this 2013 operator. |
| Calculated quantities | Frequency-dependent resistance and inductance matrices for systems of round conductors, including skin and proximity effects |
| Earth structure | Homogeneous lossless exterior only. |
| Model and approximation | Boundary electric field and equivalent surface current are truncated Fourier series. The operator eigenvalues are analytic for a solid cylinder; conductor interactions enter through the analytically integrated logarithmic Green matrix ``G``. No symmetric-current approximation is imposed. |
| Main source | U. R. Patel, B. Gustavsen, and P. Triverio (2013) |
| Citation key(s) | `:Patel2013` |
| Evidence status | Original publication page images checked |

**Description.** Boundary-only method that replaces each solid round conductor with the surrounding medium plus an equivalent longitudinal surface current, represents the boundary fields in Fourier modes, and derives the complete conductor-system series-impedance matrix.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fields are assumed longitudinally invariant; end effects are neglected. | Stated — §II-A, p. 2475. |
| Air propagation constant ``γ_air`` | Not an earth-return model; surrounding lossless medium wavenumber is ``k_{out}=\omega\sqrt{\mu_0\epsilon_{out}}``. | Stated — (15), p. 2476. |
| Earth propagation constant ``γ_earth`` | Not applicable in this paper's homogeneous lossless exterior derivation. | Stated — Fig. 1 and §II-A. |
| Earth permittivity and displacement current | Exterior displacement is retained through ``k_{out}``; no conductive earth is modeled in the published derivation. | Stated — (15). |
| Range of validity | Accuracy is controlled by conductor-specific Fourier truncation ``N_p``; examples report ``N_p=2`` or 3 as accurate, not a universal bound. | Stated — below (11), p. 2476 and §IV. |
| Earth permeability ``μ_earth`` | Not applicable; exterior permeability is ``\mu_0``. Magnetic conductors with ``\mu\ne\mu_0`` are explicitly allowed. | Stated — (5), p. 2475. |
| Arrangement | Arbitrary placement of ``P`` parallel round solid conductors. | Stated — §II-A and Fig. 1. |
| Earth structure | Homogeneous lossless exterior only. | Stated — §II-A. |
| Conductor and insulation geometry | Solid circular conductors of radii ``a_p``; insulation and hollow walls are not represented in this 2013 operator. | Stated — (6), Fig. 1. |
| Constitutive and field assumptions | Linear homogeneous isotropic conductors, harmonic fields, two-dimensional cross-section, longitudinal invariance; skin and proximity emerge from boundary harmonics. | Stated — §§II–III. |
| Conventions | ``e^{j\omega t}`` sign is implied by (1), (3), (14); ``R+j\omega L`` is a partial p.u.l. matrix in conductor coordinates. | Stated — (1), p. 2475. |

**Expression.** For conductor ``p``, the Fourier surface-admittance eigenvalue is

```math
J_n^{(p)}=E_n^{(p)}\frac{2\pi}{j\omega}
\left[
\frac{k a_p\mathcal J'_{|n|}(ka_p)}{\mu\mathcal J_{|n|}(ka_p)}
-\frac{k_{out}a_p\mathcal J'_{|n|}(k_{out}a_p)}{\mu_0\mathcal J_{|n|}(k_{out}a_p)}
\right],
\tag{13}
```

with ``k=\sqrt{\omega\mu(\omega\epsilon-j\sigma)}`` and ``k_{out}=\omega\sqrt{\mu_0\epsilon_{out}}``. After Fourier/Galerkin assembly,

```math
\mathbf R(\omega)+j\omega\mathbf L(\omega)=
\left[
\mathbf U^T(\mathbf 1-j\omega\mu_0\mathbf Y_s\mathbf G)^{-1}
\mathbf Y_s\mathbf U
\right]^{-1}.
\tag{31}
```

```math
\mathbf R=\Re\{\cdots\},\qquad
\mathbf L=\omega^{-1}\Im\{\cdots\}.
\tag{32,33}
```

**Approximation.** Boundary electric field and equivalent surface current are truncated Fourier series. The operator eigenvalues are analytic for a solid cylinder; conductor interactions enter through the analytically integrated logarithmic Green matrix ``G``. No symmetric-current approximation is imposed.

**Limitations.** The derived operator is for solid round conductors and a homogeneous lossless exterior. The output combines conductor loss with exterior magnetic field and is therefore a partial series matrix, not a scalar isolated-conductor internal impedance. Earth return, hollow conductors and insulation require later extensions.

**Reference.** [Patel2013](@cite), equations (1)–(33), printed pp. 2475–2477 (PDF pages 2–4).

**Transcription source.** Original IEEE page images. Both Bessel ratios, permeability factors, signs, inverse ordering and real/imaginary extraction were visually verified.

## Source transcription

The continuous equivalence operator is

```math
J_s=H_t-\widetilde H_t
=\frac1{j\omega}\left[
\frac1\mu\frac{\partial E_z}{\partial n}
-\frac1{\mu_0}\frac{\partial\widetilde E_z}{\partial n}
\right].
\tag{5}
```

The discretization satisfies ``\mathbf J=\mathbf Y_s\mathbf E`` and ``\mathbf I=\mathbf U^T\mathbf J``.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``J_s,E_z`` | unchanged | equivalent surface current and longitudinal boundary field | ``\mathrm{A/m}``, ``\mathrm{V/m}`` |
| ``\mathbf Y_s`` | unchanged | discrete surface-admittance operator | Fourier-basis operator |
| ``\mathbf G`` | unchanged | logarithmic exterior Green interaction matrix | source normalized |
| ``\mathbf U`` | unchanged | selector from constant current harmonics to conductor currents | dimensionless |
| ``\mathbf R,\mathbf L`` | unchanged | partial p.u.l. resistance/inductance matrices | ``\Omega/\mathrm m``, ``\mathrm H/\mathrm m`` |

No notation was renamed.

## Evidence and approximation sources

Equation (5) extends an earlier rectangular-conductor surface operator to magnetic round conductors; equations (13) and (31)–(33) are the paper's explicit round-conductor/Fourier system formulation.

## Limitations and discrepancies

- This 2013 result must not be credited with the hollow-conductor or layered-earth extensions published later.
