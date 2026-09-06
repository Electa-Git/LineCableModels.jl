# Ametani–Yoneda–Baba–Nagaoka power-frequency mixed approximation

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Parallel horizontal wires, horizontal separation ``y``, signed source coordinates ``h_1,h_2``; ``d^2=(h_1-h_2)^2+y^2`` in section II. |
| Calculated quantities | Low-frequency mutual earth-return impedance between overhead and buried conductors |
| Earth structure | Homogeneous half-space below air. |
| Model and approximation | The authors invoke a large complex penetration-depth magnitude in (30) to reduce (27) and the Lucca expression (28) to (31)–(32). The paper supplies no explicit truncation remainder or dimensionless geometry ratio in this subsection. The author-stated low-frequency result is retained even though its printed factors are suspect. |
| Main source | Ametani, Yoneda, Baba, and Nagaoka (2009), who identify this reduction with a previously known overhead approximation rather than claim priority for that older expression |
| Citation key(s) | `:Ametani2009` |
| Evidence status | Original PDF equations checked against page images; published mathematical/convention discrepancies unresolved |

**Description.** Published power-frequency reduction of the mixed overhead/underground earth-return impedance approximation, expressed in milliohms per kilometer.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | No independent imposed constant is present; inherits the source's parent TEM assumption. | Stated — sections II–IV; equation-implied — (31), p. 865. |
| Air propagation constant ``γ_air`` | Air conductivity and displacement current omitted in the parent, yielding ``m_1=0``. | Equation-implied — (4), (6), p. 861; inherited through (27) to (31). |
| Earth propagation constant ``γ_earth`` | ``m=\sqrt{j\omega\mu_0/\rho_e}=|m|e^{j\pi/4}``; ``h_e=1/m``. | Stated — (22), (27), p. 864. |
| Earth permittivity and displacement current | Omitted; conduction-only parent. | Stated — (4), p. 861; (22), p. 864. |
| Range of validity | “Power frequency region” and ``|h_e|\gg1`` as printed in (30); no reference length is supplied to make this inequality dimensionless. | Stated — section IV-C, p. 865. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0``. | Equation-implied — (22), (30). |
| Arrangement | Mixed mutual interaction; the authors also identify the printed reduction with an existing overhead expression. | Stated — section IV-C, p. 865. |
| Earth structure | Homogeneous half-space below air. | Stated — Fig. 1 and parent formulation. |
| Conductor and insulation geometry | Parallel horizontal wires, horizontal separation ``y``, signed source coordinates ``h_1,h_2``; ``d^2=(h_1-h_2)^2+y^2`` in section II. | Stated — definitions following (3), p. 861; (31), p. 865. |
| Constitutive and field assumptions | Scalar earth resistivity, nonmagnetic media, source TEM approximation. No conductor internal or insulation contribution in this output. | Stated/equation-implied — sections II–IV. |
| Conventions | ``f`` is frequency, ``\omega=2\pi f``, ``\rho_e`` earth resistivity. Equation (31) explicitly uses milliohms per kilometer. The final paragraph requests negative ``h_2`` for an underground conductor; the preceding section's positive-depth usage conflicts with it. Time factor is not separately stated. | Stated — (30)–(32) and final paragraph of section IV-C, p. 865. |

**Expression.** Equations (30)–(32) exactly as printed.

```math
|h_e|=\sqrt{\frac{\rho_e}{\omega\mu_0}}
\simeq\sqrt{\frac{\rho_e}{8f}}\times10^3\gg1,
\tag{30}

Z_m\simeq f+j\left\{8.253+0.628\ln\left(\frac{\rho_e}{fd^2}\right)\right\}
\quad\text{(in milliohms per kilometer)},
\tag{31}

Z_L\simeq Z_m+\frac{1}{12}\simeq Z_m.
\tag{32}
```

Here ``Z_L`` denotes the Lucca comparison formula in (28), ``Z_m`` the authors' (27), and ``d^2=(h_1-h_2)^2+y^2`` is the geometry defined on p. 861. The lack of an ``f`` multiplier on the imaginary bracket in (31), and the bare ``1/12`` addition in (32), are retained without correction.

**Approximation.** The authors invoke a large complex penetration-depth magnitude in (30) to reduce (27) and the Lucca expression (28) to (31)–(32). The paper supplies no explicit truncation remainder or dimensionless geometry ratio in this subsection. The author-stated low-frequency result is retained even though its printed factors are suspect.

**Limitations.** This is a verified transcription of a problematic printed reduction, not a resolved physical formula. The missing scale in ``|h_e|\gg1``, factor questions in (31)–(32), and changing burial-depth convention remain open.

**Reference.** [Ametani2009](@cite), section IV-C, equations (30)–(32), printed p. 865. Actual inspected byline: Akihiro Ametani, Tetsuzo Yoneda, Yoshihiro Baba, Naoto Nagaoka.

**Transcription source.** Original PDF page image, p. 865; the constants 8.253, 0.628, ``1/12``, the standalone ``f``, parentheses, logarithm argument, and printed output units were visually verified. The geometric definition was checked on p. 861.

## Source transcription

The formula section retains the source order and notation of (30)–(32). The parent is

```math
Z_m=j\omega\left(\frac{\mu_0}{2\pi}\right)
\exp\left(-\frac{h_2}{h_e}\right)\ln\left(\frac{S}{D}\right),
\quad h_e=\frac1m,\quad m=\sqrt{j\omega\mu_0/\rho_e},
\tag{27}

S=\sqrt{H^2+y^2},\quad D=\sqrt{(h_1+h_2)^2+y^2},\quad
H=h_1+h_2+2h_e.
```

The comparison expression is explicitly attributed to Lucca [21] by the inspected source:

```math
Z_L=j\omega\left(\frac{\mu_0}{2\pi}\right)
\left[\ln\left(\frac SD\right)
-\left(\frac23\right)\left(\frac{h_e}{S^2}\right)^3H(H^2-3y^2)\right].
\tag{28, secondary witness only}
```

This inclusion identifies the parent of (32); it does not claim that Lucca's original was verified here. Detailed parent evidence is in the [exponential-image record](../homogeneous-earth-mixed-exponential-image/Ametani2009.md).

## Notation map

Notation is unchanged.

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_m,Z_L`` | unchanged | Proposed and Lucca mutual impedance approximations | mΩ/km in (31)–(32); per-unit-length SI parent in (27)–(28) |
| ``f,\omega,j`` | unchanged | Frequency, angular frequency, imaginary unit | Hz, rad/s, dimensionless |
| ``\rho_e,\mu_0`` | unchanged | Earth resistivity, vacuum permeability | ``\Omega\,\mathrm m``, H/m |
| ``m,h_e`` | unchanged | Conductive-earth propagation constant and inverse | inverse m, m |
| ``h_1,h_2,y`` | unchanged | Source heights/depths and horizontal separation | m; signed-depth ambiguity retained |
| ``d`` | unchanged | Distance from the section II definition | ``d^2=(h_1-h_2)^2+y^2`` |
| ``H,S,D`` | unchanged | Parent image geometry as defined after (27) | m |

## Evidence and approximation sources

The source explicitly calls (31) identical to a known overhead formula and says the Carson series implemented in EMTP can be used for the mixed arrangement at low frequency with negative buried ``h_2``. This is a distinct published reduction/witness, not a claim of new historical priority. The existing LCM implementation maps only the parent mixed expression; no new ID is assigned here.

## Limitations and discrepancies

- Suspected published defect: the imaginary term in (31) lacks a frequency multiplier present in the overall ``j\omega`` parent. No replacement expression is supplied.
- Unresolved published unit treatment: (32) adds a dimensionless-looking ``1/12`` directly to the impedance; the paper supplies no unit or frequency factor for this addition.
- The source's ``|h_e|\gg1`` is an inequality involving a dimensional length without an explicit normalization.
- ``d`` in (31) and ``D`` in (27) are distinct source symbols with different definitions; neither is silently replaced by the other.
- The copied existing BibTeX entry has an incorrect/incomplete author list; the inspected four-author byline is recorded above. This conflict is shared with the parent record.
