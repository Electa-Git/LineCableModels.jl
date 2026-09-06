# Pawlik–Woodhouse–Summers full-spectrum mixed mutual admittance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite parallel thin conductors on opposite sides of a planar interface. |
| Calculated quantities | Source-labelled full-spectrum cross-boundary mutual admittance ``Y_{12}^{jp}`` |
| Earth structure | Two homogeneous half-spaces with independent conductivity, permittivity, and permeability. |
| Model and approximation | Integral representation within a full TM/TE infinite-thin-wire model retaining ``Γ``; the source-defined scalar inverse is preserved and is not generalized into entrywise inversion of a potential matrix. |
| Main source | B. Pawlik, D. Woodhouse, and T. J. Summers (2018) |
| Citation key(s) | `:Pawlik2018` |
| Evidence status | Original publication page images checked; source normalization and homogeneous/infinite-wire restrictions explicit |

**Description.** Full-spectrum cross-boundary quantity explicitly labelled by the source as per-unit-length mutual admittance for an ``n``-wire system spanning two homogeneous half-spaces.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Retained through ``I=I_Ae^{j\omega t-\Gamma z}`` and the roots ``\gamma_i^2=\Gamma^2+k_i^2``; system modes are obtained from the simultaneous equations. | Stated — p. 267, (9), (21), and discussion after (83). |
| Air propagation constant ``γ_air`` | Medium-1 ``\gamma_1^2=\Gamma^2+k_1^2``; air is a later specialization, not imposed in (90)–(91). | Stated — (9), (39), §V.A. |
| Earth propagation constant ``γ_earth`` | Medium-2 ``\gamma_2^2=\Gamma^2+k_2^2``. | Stated — (21), (40). |
| Earth permittivity and displacement current | Retained in both media through ``k_i=\omega\sqrt{\mu_i(\epsilon_i-j\sigma_i/\omega)}``; the prefactor in (90) uses ``\sigma_1+j\omega\epsilon_1``. | Stated — nomenclature and (90). |
| Range of validity | Infinite thin wires; physical use requires current to decay before line ends. The source warns that homogeneous conductivity strongly limits power-system applications. No universal frequency bound is supplied. | Stated — prose after (91), p. 273. |
| Earth permeability ``μ_earth`` | Independent ``\mu_2`` retained. | Stated — abstract and (91). |
| Arrangement | Mutual cross-boundary term with ``j`` in medium 1 and ``p`` in medium 2; reciprocal placement is stated by index/material interchange. | Stated — before (88) and after (91). |
| Earth structure | Two homogeneous half-spaces separated by a plane interface. | Stated — Fig. 1 and §II.A. |
| Conductor and insulation geometry | Infinitely long thin circular wires, with finite conductor conductance allowed in the full system; no finite insulation shell enters (90)–(91). | Stated — abstract and §II.A. |
| Constitutive and field assumptions | Linear isotropic homogeneous media; coupled TM/TE full-spectrum solution; interface boundary conditions. | Stated — §§II–III. |
| Conventions | ``e^{j\omega t-\Gamma z}``; ``h_j,h_p`` positive into each medium; ``d_{jp}`` horizontal; decaying transverse-root selection. | Stated — p. 267 and (90)–(91). |

**Expression.** The source labels the following scalar quantity as cross-boundary mutual admittance:

```math
Y_{12}^{jp}=\pi(\sigma_1+j\omega\epsilon_1)
\left[N_{12}^{jp}-jM_{12}^{jp}\right]^{-1},
\tag{90}
```

```math
N_{12}^{jp}-jM_{12}^{jp}
=\int_0^\infty
\frac{\left[u_1+(\mu_2/\mu_1)u_2\right]
e^{-(u_1h_j+u_2h_p)}}
{\left[u_1+(\mu_1/\mu_2)u_2\right]
\left[(k_2^2/k_1^2)u_1+(\mu_2/\mu_1)u_2\right]}
\cos(\lambda d_{jp})\,d\lambda.
\tag{91}
```

Here ``k_i=\omega\sqrt{\mu_i(\epsilon_i-j\sigma_i/\omega)}``, ``\gamma_i^2=\Gamma^2+k_i^2``, and ``u_i=\sqrt{\lambda^2-\gamma_i^2}``.

**Approximation.** Not an analytical approximation within the paper's full-spectrum infinite-thin-wire and homogeneous-half-space model. The scalar inverse is part of the source's printed definition; it is documented, not generalized into entrywise inversion of a separately assembled potential matrix.

**Limitations.** The paper prints a scalar inverse for the cross-pair coefficient and calls the result mutual admittance. It does not restate this pair as a Maxwell-potential-matrix entry or give the ``Y=j\omega P^{-1}`` assembly used by Martins-Britto 2024. Those representations are reported separately rather than silently equated element by element.

**Reference.** [Pawlik2018](@cite), equations (90)–(91), printed p. 273.

**Transcription source.** Original IEEE page image. The prefactor medium index, both permeability ratios, ``k_2^2/k_1^2`` factor, inverse exponent, and cosine term were visually verified.

## Source transcription

The source introduces (90)–(91) immediately after its cross-boundary impedance (88)–(89), under the heading “Mutual Impedance and Admittance.” It then states that the reversed medium/conductor placement follows by interchanging the relevant subscripts.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Y_{12}^{jp}`` | unchanged | source-labelled mixed mutual admittance | ``\mathrm S/\mathrm m`` |
| ``N_{12}^{jp}-jM_{12}^{jp}`` | unchanged | source spectral coefficient inverted in (90) | source normalization |
| ``h_j,h_p,d_{jp}`` | unchanged | distances into the two media and horizontal separation | m |
| ``k_i,\gamma_i,u_i`` | unchanged | material, longitudinally modified, and transverse wave numbers | ``\mathrm m^{-1}`` |

No notation was renamed.

## Evidence and approximation sources

This is direct evidence that mixed external admittance was not empty in the accessible collection. Martins-Britto et al. 2024 state that their equation (6), after matrix conversion (7), agrees with Pawlik's mixed admittance. Martins-Britto's potential-coefficient representation remains a distinct record because it explicitly supports arbitrary multiconductor matrix assembly.

## Limitations and discrepancies

- Do not reinterpret (90) as permission to replace an entry of ``P^{-1}`` by ``1/P_{ij}``; it is a source-defined scalar relation in Pawlik's formulation.
- No separate insulation-potential matrix is included in (90)–(91).
- The homogeneous two-half-space and infinite-wire/end-effect restrictions remain.
