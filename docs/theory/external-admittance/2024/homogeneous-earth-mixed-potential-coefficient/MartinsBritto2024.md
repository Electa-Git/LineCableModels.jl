# Martins-Britto–Papadopoulos–Chrysochos mixed potential coefficient

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Infinite parallel thin conductors, one in air and one in homogeneous soil. |
| Calculated quantities | Reciprocal mixed mutual Maxwell potential coefficient and ``Y_{tot}=j\omega P_{tot}^{-1}`` after full matrix assembly |
| Earth structure | Homogeneous soil half-space below air. |
| Model and approximation | Integral representation within the paper's quasi-TEM model; longitudinal propagation, both permeabilities, conductivity, and displacement current are retained. Numerical examples separately set ``k_x=0``. |
| Main source | A. G. Martins-Britto, T. A. Papadopoulos, and A. I. Chrysochos (2024) |
| Citation key(s) | `:MartinsBritto2024` |
| Evidence status | Original publication page images checked; matrix assembly explicit and no entrywise reciprocal introduced |

**Description.** Reciprocal mixed air/soil per-unit-length Maxwell potential coefficient for two infinite parallel thin conductors, with the source's required assembly of all potential entries and internal-insulation coefficients before matrix inversion to line admittance.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | Retained as ``\gamma_x=jk_x``; numerical examples set ``k_x=0``. | Stated — p. 984 and §III, p. 985. |
| Air propagation constant ``γ_air`` | ``\gamma_0=\sqrt{j\omega\mu_0(\sigma_0+j\omega\varepsilon_0)}``; retained in ``a_0`` and (6). | Stated — p. 984. |
| Earth propagation constant ``γ_earth`` | ``\gamma_1=\sqrt{j\omega\mu_1(\sigma_1+j\omega\varepsilon_1)}``; retained in ``a_1`` and (6). | Stated — p. 984. |
| Earth permittivity and displacement current | Retained in ``\gamma_1`` and hence the mixed potential kernel. | Stated/equation-implied — definition p. 984 and (6). |
| Range of validity | Quasi-TEM derivation; Appendix A reports prior support up to 10 MHz for dominant transmission-line modes. Validation here spans 1 kHz–1 MHz. These are stated/tested ranges, not universal guarantees. | Stated — Appendix A p. 990 and §III.A p. 985. |
| Earth permeability ``μ_earth`` | Independent ``\mu_1`` retained. | Stated — (6). |
| Arrangement | Mixed mutual air/soil coefficient with reciprocity ``P_{eij}^{01}=P_{eij}^{10}``; multiconductor assembly supported. | Stated — (6)–(7). |
| Earth structure | Homogeneous soil half-space below homogeneous air. | Stated — Fig. 1 and §II.A. |
| Conductor and insulation geometry | Infinite parallel thin conductors, ``h_i>0``, ``h_j<0``, lateral spacing ``y_{ij}``; coated-conductor insulation enters only through the separate internal potential matrix. | Stated — pp. 984–985. |
| Constitutive and field assumptions | Linear homogeneous media; quasi-TEM Hertzian-vector solution under Lorenz gauge. The full potential matrix includes same-medium external entries and internal insulation terms. | Stated — text before (5), Appendix A, and prose around (7). |
| Conventions | Vertical ``z`` upward, interface ``z=0``; source/observation superscripts ``0`` air and ``1`` soil; ``Y_{tot}=j\omega P_{tot}^{-1}``. | Stated — Fig. 1, pp. 984–985. |

**Expression.** The mixed coefficient is

```math
P_{eij}^{01}=P_{eij}^{10}
=-\frac{\omega^2\mu_0\mu_1}{\pi}
\int_0^\infty
\frac{(a_0\mu_0+a_1\mu_1)e^{-a_0h_i+a_1h_j}}
{(a_1\mu_1\gamma_0^2+a_0\mu_0\gamma_1^2)
 (a_0\mu_1+a_1\mu_0)}
\cos(\lambda y_{ij})\,d\lambda.
\qquad\text{(6)}
```

After assembling same-medium terms (2a), (4a), the mixed terms (6), and the internal potential matrix, the source requires

```math
\mathbf Y_{tot}=j\omega\mathbf P_{tot}^{-1}.
\qquad\text{(7)}
```

Here ``\gamma_k=\sqrt{j\omega\mu_k(\sigma_k+j\omega\varepsilon_k)}``, ``a_k=\sqrt{\lambda^2+\gamma_k^2+k_x^2}``, and ``\gamma_x=jk_x``.

**Approximation.** Not an analytical approximation after the source's quasi-TEM, infinite-thin-conductor, homogeneous-medium assumptions. Setting ``k_x=0`` is a separate numerical prescription used in the examples, not applied to the displayed general coefficient.

**Limitations.** Equation (6) is a potential-matrix entry, not an admittance entry. The source requires inversion of assembled ``\mathbf P_{tot}``; ``(P^{-1})_{ij}`` is not replaced by ``1/P_{ij}``, and the earth-only mixed coefficient is not independently inverted.

**Reference.** [MartinsBritto2024](@cite), equations (6)–(7), printed p. 985; definitions p. 984; Appendix derivation pp. 990–991.

**Transcription source.** Original IEEE page images. The leading minus sign, ``\omega^2\mu_0\mu_1/\pi`` factor, both denominator products, signed exponent, reciprocity, and matrix inversion were visually verified.

## Source transcription

Appendix B defines ``P_{eij}^{01}/(j\omega)`` by longitudinal integration of the source's scalar-potential-related function in (B.5), expands it in (B.7), and obtains (6) after (B.14) and the variable transform. Section II.B states that (6) combined with (7) agrees with Pawlik's mixed admittance, while emphasizing that the potential coefficient supports generalized multiconductor assembly.

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``P_{eij}^{01},P_{eij}^{10}`` | unchanged | reciprocal mixed mutual potential coefficients | p.u.l. potential-matrix normalization |
| ``\mathbf P_{tot}`` | unchanged | assembled external plus internal potential matrix | invert as a matrix |
| ``\mathbf Y_{tot}`` | unchanged | total p.u.l. shunt-admittance matrix | ``j\omega\mathbf P_{tot}^{-1}`` |
| ``h_i,h_j,y_{ij}`` | unchanged | signed vertical coordinates and horizontal spacing | m |
| ``a_k,\gamma_k,k_x`` | unchanged | transverse, material, and longitudinal spectral quantities | ``\mathrm m^{-1}`` |

No notation was renamed.

## Evidence and approximation sources

This record fills the mixed potential-coefficient coverage directly from an original publication. It is distinct from Pawlik's source-labelled scalar mutual admittance because it exposes the coefficient to be inserted into a complete ``P`` matrix and prints the matrix-to-admittance operation.

## Limitations and discrepancies

- Homogeneous two-medium geometry only; no layered mixed-potential recursion is supplied.
- Numerical examples set ``k_x=0`` and do not establish a universal truncation/frequency theorem.
- Agreement with Pawlik is source-stated; the corpus does not algebraically collapse the two differently normalized formulations.

## Numerical interpretation

Numerical evaluation retains independent half-space permeabilities and the prescribed longitudinal input. Common material scales cancel in the magnetic and electric denominator factors before quadrature. The same factored expression serves same-medium and mixed pairs, with their different propagation distances retained. Tests include reversed pairs and negative-frequency conjugacy. The returned coefficients enter the complete potential matrix before inversion.
