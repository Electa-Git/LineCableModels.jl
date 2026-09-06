# Ametani–Miyamoto–Mahseredjian air-referenced displacement-current integral

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel round overhead conductors; mutual uses direct/image distances; self uses the conductor radius in ``d_{ii}``. |
| Calculated quantities | Self and mutual overhead earth-return impedance retaining earth displacement current relative to the air bulk constant |
| Earth structure | Homogeneous earth half-space below air, obtained as the homogeneous limit of (9). |
| Model and approximation | No finite closed-form approximation is applied in this record; it is the source's homogeneous spectral integral. It does impose the longitudinal/air-reference prescription embodied in ``m_1^2-m_0^2``. |
| Main source | A. Ametani, Y. Miyamoto, and J. Mahseredjian (2014), identifying the homogeneous limit of a stratified-earth formulation and its relation to Wise/Kikuchi |
| Citation key(s) | `:Ametani2014` |
| Evidence status | Original publication page images checked |

**Description.** Carson-form overhead integral in which the earth vertical constant is formed as the difference between earth and air bulk wavenumber squares, retaining conductivity, earth permittivity and air displacement current.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Imposed as the air bulk propagation constant through the subtraction ``m_1^2-m_0^2``. | Stated — (11)–(13), p. 937. |
| Air propagation constant ``γ_air`` | ``m_0^2=j\omega\mu_0(j\omega\epsilon_0)`` for ``\sigma_0=0``. | Stated — (13), p. 937. |
| Earth propagation constant ``γ_earth`` | ``m_1^2=j\omega\mu_e(\sigma_e+j\omega\epsilon_e)``. | Stated — (13). |
| Earth permittivity and displacement current | Retained explicitly in ``m_1``; the spectral root uses the air-referenced difference. | Stated — (11)–(13). |
| Range of validity | The source calls (12) the most general form within this overhead integral reduction; numerical comparisons cover its selected 100 Hz–1 MHz cases but give no universal bound. | Stated — §§2.7–4. |
| Earth permeability ``μ_earth`` | Equation (11) retains ``\mu_e``; the compact ``A=1/(s+a_1)`` in (12) assumes ``\mu_e=\mu_0``. | Stated — text before (12). |
| Arrangement | Overhead multiconductor; self and mutual. | Stated — Fig. 1 and (1). |
| Earth structure | Homogeneous earth half-space below air, obtained as the homogeneous limit of (9). | Stated — text before (11). |
| Conductor and insulation geometry | Parallel round overhead conductors; mutual uses direct/image distances; self uses the conductor radius in ``d_{ii}``. | Stated — Fig. 1 and (2). |
| Constitutive and field assumptions | Linear homogeneous isotropic media; transmission-line/Carson-form reduction rather than full TM/TE modal solution. | Stated — §§1–2. |
| Conventions | ``j=\sqrt{-1}``; ``s`` is the spectral variable; ``P_0=\ln(D_{ij}/d_{ij})``; p.u.l. impedance. | Stated — (1)–(4). |

**Expression.** With ``\mu_e=\mu_0`` and ``\sigma_0=0``, the homogeneous spectral factor is

```math
A=\frac{1}{s+a_1},\qquad
a_1=\sqrt{s^2+m_1^2-m_0^2},
\tag{12}
```

```math
m_1^2=j\omega\mu_e(\sigma_e+j\omega\epsilon_e),
\qquad
m_0^2=j\omega\mu_0(j\omega\epsilon_0).
\tag{13}
```

It enters

```math
Z_{ij}=j\omega\frac{\mu_0}{2\pi}[P_0+(Q-jR)],
\quad P_0=\ln\frac{D_{ij}}{d_{ij}},
\tag{1,2}
```

```math
Q-jR=2\int_0^\infty A,e^{-(h_i+h_j)s}\cos(ys)\,ds.
\tag{3,4}
```

**Approximation.** No finite closed-form approximation is applied in this record; it is the source's homogeneous spectral integral. It does impose the longitudinal/air-reference prescription embodied in ``m_1^2-m_0^2``.

**Limitations.** Homogeneous planar earth, infinite parallel thin conductors and the stated transmission-line reduction. The paper separately notes that a complete Kikuchi formulation also covers TM/TE modes; this integral does not thereby become full wave. Root branches are not explicitly printed.

**Reference.** [Ametani2014](@cite), equations (1)–(4), (9)–(14), printed pp. 936–937 (PDF pages 1–2).

**Transcription source.** Original IEEJ page images. The air subtraction, bulk constants, spectral prefactor, geometric logarithm and Fourier kernel were visually verified.

## Source transcription

For arbitrary permeability the preceding homogeneous form is

```math
A=A_1=\left[s+\frac{\mu_0}{\mu_1}\sqrt{s^2+m_1^2-m_0^2}\right]^{-1}.
\tag{11}
```

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{ij}`` | unchanged | earth-return self/mutual impedance | ``\Omega/\mathrm m`` |
| ``s`` | unchanged | transverse spectral variable | ``\mathrm m^{-1}`` |
| ``m_0,m_1`` | unchanged | air and earth intrinsic bulk constants | ``\mathrm m^{-1}`` |
| ``a_1`` | unchanged | air-referenced earth vertical root | ``\mathrm m^{-1}`` |
| ``d_{ij},D_{ij}`` | unchanged | direct and image distances | m |

No notation was renamed.

## Evidence and approximation sources

The source derives (11) as the homogeneous limit of the stratified expression (9), then identifies (12) with Wise and the corresponding Kikuchi limit. This record follows the printed derivation and does not assign priority beyond that attribution.

## Limitations and discrepancies

- Equation (14) visibly omits the ``j\omega`` multiplying ``\epsilon_0(\epsilon_r-1)`` inside its brace, although equations (8), (12)–(13), and (18) retain it. The corpus uses the unsimplified authoritative (12)–(13) and records, but does not repair, the (14) conflict.
