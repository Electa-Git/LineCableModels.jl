# Patel–Gustavsen–Triverio hollow-conductor MoM–SO impedance

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Solid cylinders and annuli with inner and outer Fourier boundaries; infinitesimally thin insulation assumed in the ground correction. |
| Calculated quantities | P.u.l. series resistance and inductance matrices for arbitrary systems of solid and hollow round conductors with skin and proximity effects |
| Earth structure | Homogeneous lossless operator exterior; homogeneous-ground contribution may be appended approximately. |
| Model and approximation | The boundary fields/currents are truncated Fourier series. Equations (37)–(38) are the direct conductor/exterior result; (39)–(40) adds a proximity correction to a separately computed conventional ground-return matrix and is explicitly approximate. |
| Main source | U. R. Patel, B. Gustavsen, and P. Triverio (2014) |
| Citation key(s) | `:Patel2014` |
| Evidence status | Original publication page images checked |

**Description.** Extension of the round-conductor surface-admittance method to annular conductors, with independent inner/outer boundary harmonics, enabling proximity-aware series matrices for sheaths, armors and pipes.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fields are longitudinally invariant; end effects are excluded. | Stated — §II. |
| Air propagation constant ``γ_air`` | Exterior/hole medium is lossless with its own wavenumber; not an exact earth-return operator. | Stated — §§II–IV. |
| Earth propagation constant ``γ_earth`` | Ground is appended through a conventional calculation, not solved in the hollow surface operator. | Stated — §V. |
| Earth permittivity and displacement current | Exterior lossless displacement is retained in the operator; the ground correction uses the conventional earth model selected by the user. | Stated — §§II, V. |
| Range of validity | Fourier orders ``N_p=3`` or 4 are reported as typical; simplified ground addition requires earth penetration depth much larger than conductor spacing. | Stated — §§II, V. |
| Earth permeability ``μ_earth`` | Not part of the conductor operator; conductor permeability and exterior ``\mu`` are retained. | Stated — §II. |
| Arrangement | Arbitrary parallel systems of solid and hollow round conductors. | Stated — abstract, §II. |
| Earth structure | Homogeneous lossless operator exterior; homogeneous-ground contribution may be appended approximately. | Stated — §§II, V. |
| Conductor and insulation geometry | Solid cylinders and annuli with inner and outer Fourier boundaries; infinitesimally thin insulation assumed in the ground correction. | Stated — §§II–III, V. |
| Constitutive and field assumptions | Linear homogeneous isotropic conductors, 2-D harmonic fields, equivalence theorem, truncated Fourier/Galerkin discretization. | Stated — §§II–IV. |
| Conventions | ``R+j\omega L`` p.u.l.; current selector includes opposite-signed constant harmonics at the inner and outer annular boundaries. | Stated — before (37). |

**Expression.** With the hollow-conductor surface-admittance operator ``\mathbf Y_s`` and exterior interaction matrix ``\mathbf G``,

```math
\mathbf R(\omega)=\Re\left\{\left[\mathbf U^T
(\mathbf 1-j\omega\mu_0\mathbf Y_s\mathbf G)^{-1}\mathbf Y_s\mathbf U
\right]^{-1}\right\},
\qquad\text{(37)}
```

```math
\mathbf L(\omega)=\omega^{-1}\Im\left\{\left[\mathbf U^T
(\mathbf 1-j\omega\mu_0\mathbf Y_s\mathbf G)^{-1}\mathbf Y_s\mathbf U
\right]^{-1}\right\}.
\qquad\text{(38)}
```

For the approximate ground inclusion,

```math
\mathbf Z=(\mathbf Z_c+\mathbf Z_g)+\Delta\mathbf Z_{prox},\qquad
\Delta\mathbf Z_{prox}=\mathbf Z_{MoM-SO}(N_p>0)-\mathbf Z_{MoM-SO}(N_p=0).
\qquad\text{(39,40)}
```

**Approximation.** The boundary fields/currents are truncated Fourier series. Equations (37)–(38) are the direct conductor/exterior result; (39)–(40) adds a proximity correction to a separately computed conventional ground-return matrix and is explicitly approximate.

**Limitations.** Circular parallel conductors only. The ground decomposition is valid only when earth skin depth is much larger than conductor spacing and insulation thickness is negligible relative to that depth. It is not the later exact multilayer Green-function extension.

**Reference.** [Patel2014](@cite), equations (1)–(44), especially (37)–(40), printed pp. 2105–2106.

**Transcription source.** Original IEEE page images. Matrix nesting, inverse order, real/imaginary operations, ``N_p`` subtraction and the ground-validity statement were visually verified.

## Source transcription

The annular ``\mathbf Y_s`` derives from the two-boundary Bessel system (11)–(28). Its inner and outer constant harmonics both enter the conductor-current selector; this structural difference from the 2013 solid-cylinder operator is essential.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``\mathbf Y_s`` | unchanged | block surface-admittance operator for solid/annular conductors | boundary operator |
| ``\mathbf G`` | unchanged | lossless exterior Green interaction matrix | source normalized |
| ``N_p`` | unchanged | boundary harmonic order of conductor ``p`` | integer |
| ``\Delta\mathbf Z_{prox}`` | unchanged | proximity-only impedance correction | ``\Omega/\mathrm m`` |

No notation was renamed.

## Evidence and approximation sources

The hollow-conductor operator is the paper's new equation-level contribution. Its ground term is explicitly a correction construction and is not described here as a new earth kernel.

## Limitations and discrepancies

- Ground return in (39) depends on the separately selected conventional model, so (39) is not independently reproducible without that source.
