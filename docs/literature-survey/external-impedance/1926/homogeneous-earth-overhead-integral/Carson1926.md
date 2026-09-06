# Carson homogeneous-earth overhead correction integral

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Infinitely long parallel wires; self expression includes wire radius ``a`` only in the perfect-ground base term, while the finite-ground correction uses height. |
| Calculated quantities | Per-unit-length self and mutual finite-conductivity ground-return corrections for overhead wires |
| Earth structure | Plane homogeneous semi-infinite ground ``y\leq0`` beneath a nonconducting dielectric ``y>0``. |
| Model and approximation | Integral representation exact within Carson's stated reduced field model. The reductions are the very-small-``\Gamma`` assumption and neglect of transverse ground electric-field components; earth displacement current and an independent permeability are absent. Carson's later series evaluations of ``J`` are distinct evaluator forms and are not substituted here. |
| Main source | John R. Carson (1926) |
| Citation key(s) | `:Carson1926` |
| Evidence status | PDF page images checked |

**Description.** Carson's per-unit-length correction to the self and mutual series impedance of parallel overhead wires caused by finite conductivity of a homogeneous semi-infinite ground, in the source's electromagnetic c.g.s. notation.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Fields initially contain ``\exp(-\Gamma z+i\omega t)`` and ``\Gamma`` is assumed very small; it is absent from the final correction integral. | Stated — printed p. 539; equations (27)–(31), p. 545. |
| Air propagation constant ``γ_air`` | Not defined; the upper dielectric is assigned zero conductivity and no air propagation term occurs in the correction. | Stated — printed p. 539. |
| Earth propagation constant ``γ_earth`` | No independently named ``\gamma_{earth}``; conduction enters through ``\alpha=4\pi\lambda\omega`` and ``\sqrt{\mu^2+i}`` after normalization. | Stated — printed p. 540 and (27)–(29), p. 545. |
| Earth permittivity and displacement current | Earth permittivity does not occur; the ground is represented by conductivity ``\lambda``. | Equation-implied — definition of ``\alpha`` and field equation (1), p. 540. |
| Range of validity | ``\Gamma`` is assumed very small in electromagnetic c.g.s. units; the ground transverse field components ``E_x,E_y`` are assumed negligible compared with axial ``E_z``. No universal frequency bound is supplied. | Stated — printed pp. 539–540. |
| Earth permeability ``μ_earth`` | No independent earth-permeability parameter is retained in the source expression. | Equation-implied — ``\alpha=4\pi\lambda\omega`` in the declared c.g.s. formulation, p. 540. |
| Arrangement | Overhead; self term at height ``h`` and mutual term for parallel wires at heights ``h_1,h_2`` separated horizontally by ``x``. | Stated — printed pp. 539, 544–545. |
| Earth structure | Plane homogeneous semi-infinite ground ``y\leq0`` beneath a nonconducting dielectric ``y>0``. | Stated — printed p. 539. |
| Conductor and insulation geometry | Infinitely long parallel wires; self expression includes wire radius ``a`` only in the perfect-ground base term, while the finite-ground correction uses height. | Stated — (23)–(28), printed pp. 544–545. |
| Constitutive and field assumptions | Homogeneous conductive ground; axial ground electric field retained, transverse components neglected; linear harmonic fields. | Stated — printed pp. 539–540. |
| Conventions | Source coordinates place the wire parallel to ``z`` and ground below ``y=0``; common factor ``\exp(-\Gamma z+i\omega t)``; ``i=\sqrt{-1}``; all displayed formulas use electromagnetic c.g.s. units. | Stated — printed pp. 539–540. |

**Expression.** Self and mutual finite-conductivity ground-return correction integrals, equations (27)–(31).

```math
Z=Z^0+Z',\qquad Z_{12}=Z_{12}^{0}+Z'_{12},
\qquad\text{(25--26)}

Z'=4\omega\int_0^\infty
\left(\sqrt{\mu^2+i}-\mu\right)e^{-2h'\mu}\,d\mu,
\qquad\text{(27)}

Z'_{12}=4\omega\int_0^\infty
\left(\sqrt{\mu^2+i}-\mu\right)
e^{-(h'_1+h'_2)\mu}\cos(x'\mu)\,d\mu,
\qquad\text{(28)}

J(p,q)=\int_0^\infty
\left(\sqrt{\mu^2+i}-\mu\right)e^{-p\mu}\cos(q\mu)\,d\mu,
\qquad\text{(29)}

Z'=4\omega J(2h',0),\qquad
Z'_{12}=4\omega J(h'_1+h'_2,x'),
\qquad\text{(30--31)}

\alpha=4\pi\lambda\omega,\qquad
h'=h\sqrt{\alpha},\quad h'_1=h_1\sqrt{\alpha},\quad
h'_2=h_2\sqrt{\alpha},\quad x'=x\sqrt{\alpha}.
```

``Z^0`` and ``Z_{12}^0`` are the corresponding self and mutual impedances for perfectly conducting ground. ``\lambda`` is ground conductivity in electromagnetic c.g.s. units; ``\mu`` is the dimensionless integration variable and is not permeability. The source's correction is reported in its c.g.s. impedance-per-length convention.

**Approximation.** Integral representation exact within Carson's stated reduced field model. The reductions are the very-small-``\Gamma`` assumption and neglect of transverse ground electric-field components; earth displacement current and an independent permeability are absent. Carson's later series evaluations of ``J`` are distinct evaluator forms and are not substituted here.

**Limitations.** Homogeneous plane ground and infinite parallel overhead wires only. The displayed ``Z'`` terms are finite-ground corrections, not total line impedances. The source uses electromagnetic c.g.s. units; this record performs no hidden SI rescaling. The square-root branch is not explicitly stated beside (27)–(29).

**Reference.** [Carson1926](@cite), equations (25)–(31), printed p. 545 (PDF page 7), with assumptions and definitions on printed pp. 539–540.

**Transcription source.** Original publication scan. Integral limits, the two separate self/mutual exponentials, cosine, normalization, and all factors were checked visually against printed p. 545. The scan has no usable text layer; no OCR reconstruction was treated as authority.

## Source transcription

The formula section retains Carson's source notation. The immediately preceding total impedances are

```math
Z=z+i2\omega\log(\rho''/a)+Z',
\qquad\text{(23)}

Z_{12}=i2\omega\log(\rho''/\rho')+Z'_{12},
\qquad\text{(24)}
```

with ``\rho''=\sqrt{(h_1+h_2)^2+x^2}`` and ``\rho'=\sqrt{(h_1-h_2)^2+x^2}``. They show explicitly that (27)–(28) are corrections rather than complete conductor impedances.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z',Z'_{12}`` | unchanged | Self/mutual finite-ground corrections | per unit length, source c.g.s. system |
| ``Z^0,Z^0_{12}`` | unchanged | Perfectly conducting-ground base impedances | per unit length |
| ``J(p,q)`` | unchanged | Dimensionless normalized Carson integral | source definition (29) |
| ``\mu`` | unchanged | Integration variable | dimensionless; not permeability |
| ``\lambda`` | unchanged | Ground conductivity | electromagnetic c.g.s. units |
| ``\alpha`` | unchanged | Normalization ``4\pi\lambda\omega`` | inverse-length squared in source normalization |
| ``h,h_1,h_2`` | unchanged | Wire heights | source length units |
| ``x`` | unchanged | Horizontal separation of wire planes | source length units |
| primed geometry | unchanged | Geometry multiplied by ``\sqrt\alpha`` | dimensionless |
| ``i`` | unchanged | Imaginary unit | ``i^2=-1``; time factor ``e^{i\omega t}`` |
| ``\Gamma`` | unchanged | Longitudinal propagation constant | common factor ``e^{-\Gamma z}`` |

## Evidence and approximation sources

Problem geometry, half-spaces, and longitudinal/time convention: p. 539. Conductivity normalization and c.g.s. declaration: p. 540. Total/self/mutual decomposition and normalized correction integrals: (23)–(31), pp. 544–545.

## Limitations and discrepancies

- All available copies are scans with no useful mathematical text layer; visual transcription was mandatory.
- This record does not convert Carson's c.g.s. coefficient ``4\omega`` into an SI kernel.

