# Ametani–Miyamoto–Nagaoka bonded two-layer cylindrical conductor

## Identification and source

| Field | Value |
| --- | --- |
| Family | Internal impedance |
| Geometry | Medium 1 occupies ``a\le r\le b`` and medium 2 occupies ``b'\le r\le c``; for bonded media ``b=b'``. The first medium may be solid (``a=0``) or hollow. Both ends are short-circuited for the outer-terminal reduction. |
| Calculated quantities | Inner-surface, outer-surface, and transfer impedance of two electrically bonded concentric conductive media; stable assembly from single-layer surface terms and direct Maxwell/Bessel representation |
| Earth structure | Not applicable. |
| Model and approximation | Appendix (A.11) is the direct Maxwell/Bessel representation with ``\Gamma`` and intrinsic admittance retained. The main-text construction (9)–(10) is algebraically equivalent under the paper's stated good-conductor and penetration-depth assumptions: ``j\omega\mu\sigma\gg\omega^2\mu\varepsilon+\Gamma^2`` and no radial current between media. It is a model reduction, not merely a notation change. |
| Main source | A. Ametani, Y. Miyamoto, and N. Nagaoka (2004); appendix says the direct Maxwell derivation summarizes N. Amekawa's 2001 thesis |
| Citation key(s) | `:Ametani2004` |
| Evidence status | PDF page images checked |

**Description.** Inner, outer, and transfer surface impedances of two concentric conductive cylindrical media in electrical contact, developed for metallic conductor/semiconducting-screen combinations and also applicable to a semiconductor on either surface of a cable sheath.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | The appendix retains ``\Gamma`` in ``m^2=j\omega\mu\eta-\Gamma^2``. The reduced main-text surface terms use the good-conductor condition (A.12), which omits ``\Gamma^2``. | Stated/equation-implied — (A.2)–(A.3), (A.12), printed pp. 1530–1531. |
| Air propagation constant ``γ_air`` | Not applicable. | The formulation is confined to conductive cylindrical layers. |
| Earth propagation constant ``γ_earth`` | Not applicable. | The formulation is earth-independent. |
| Earth permittivity and displacement current | Not applicable. | The formulation is earth-independent. |
| Range of validity | For the independent-layer circuit construction the source requires penetration depth ``d`` greater than semiconductor thickness ``c-b'`` and says the condition is almost always satisfied below 1 GHz for the cited 1–5 mm layers. The direct appendix expression precedes the good-conductor reduction. | Stated — (1) and following text, printed p. 1524; (A.12), printed p. 1531. |
| Earth permeability ``μ_earth`` | Not applicable. | The formulation is earth-independent. |
| Arrangement | Not applicable. | This is a local cross-sectional conductor surface relation. |
| Earth structure | Not applicable. | The formulation is earth-independent. |
| Conductor and insulation geometry | Medium 1 occupies ``a\le r\le b`` and medium 2 occupies ``b'\le r\le c``; for bonded media ``b=b'``. The first medium may be solid (``a=0``) or hollow. Both ends are short-circuited for the outer-terminal reduction. | Stated — Fig. 1 and (5)–(10), printed pp. 1523–1524. |
| Constitutive and field assumptions | Linear isotropic media with ``\eta=\sigma+j\omega\varepsilon`` in the appendix; sinusoidal, axially invariant, circularly symmetric fields. The main circuit construction additionally assumes no radial current across the interface and treats each material as an isolated Schelkunoff cylindrical conductor before parallel assembly. | Stated — (1), Fig. 10, and (A.1)–(A.3), printed pp. 1524, 1530. |
| Conventions | ``d/dt=j\omega`` and fields vary as ``e^{j\omega t}``; infinitely long conductor on the ``z`` axis; inner and outer surface currents/orientations are those of (A.8)–(A.10); impedances are per unit length in ``\Omega/\mathrm m``. | Stated — appendix opening and (A.8)–(A.10), printed pp. 1530–1531. |

**Expression.** Source equations (8)–(10), assembled from the component surface impedances in (3); the direct Maxwell witness is (A.11)–(A.13).

```math
Z_{out}=\frac{V}{I}
=\frac{Z_{11}Z_{22}-Z_{12}^{2}}
       {Z_{11}+Z_{22}-2Z_{12}}
=z_{2o}-\frac{z_{2m}^{2}}{z_{1o}+z_{2i}},
\qquad\text{(8--9)}

Z_{in}=z_{1i}-\frac{z_{1m}^{2}}{z_{1o}+z_{2i}},
\qquad
Z_m=\frac{z_{1m}z_{2m}}{z_{1o}+z_{2i}},
\qquad\text{(10)}

Z_{in}=\frac{m_1\rho_1}{2\pi aD}
        \left(m_1\rho_1FQ+m_2\rho_2EP\right),

Z_{out}=\frac{m_2\rho_2}{2\pi cD}
         \left(m_1\rho_1GR+m_2\rho_2HS\right),

Z_m=\frac{\rho_1\rho_2}{2\pi abcD},
\qquad\text{(A.11)}

m_n^2=j\omega\mu_n\eta_n-\Gamma^2,
\qquad \eta_n=\sigma_n+j\omega\varepsilon_n,
\qquad\text{(A.2--A.3)}

m_n^2\simeq j\omega\mu_n\sigma_n=\frac{j\omega\mu_n}{\rho_n}
\quad\text{under (A.12).}
```

The complete source definitions of ``D,E,F,G,H,P,Q,R,S`` are retained below. In the circuit form, every ``z_{ni}``, ``z_{no}``, and ``z_{nm}`` is the corresponding single-layer inner, outer, or transfer surface impedance. The source obtains (9) after setting ``b=b'`` and therefore ``z_{12}=0``.

**Approximation.** Appendix (A.11) is the direct Maxwell/Bessel representation with ``\Gamma`` and intrinsic admittance retained. The main-text construction (9)–(10) is algebraically equivalent under the paper's stated good-conductor and penetration-depth assumptions: ``j\omega\mu\sigma\gg\omega^2\mu\varepsilon+\Gamma^2`` and no radial current between media. It is a model reduction, not merely a notation change.

**Limitations.** Perfect concentric circular geometry and axial invariance exclude sector, strand, and proximity effects. The no-radial-current assumption is tied to ``d>c-b'``. The paper states a practical below-1-GHz observation for its cited semiconductor dimensions, not a universal frequency guarantee. The stable circuit expression requires independently defined single-layer surface terms.

**Reference.** [Ametani2004](@cite).  A. Ametani, Y. Miyamoto, and N. Nagaoka, “Semiconducting Layer Impedance and its Effect on Cable Wave-Propagation and Transient Characteristics,” *IEEE Transactions on Power Delivery* 19(4), 1523–1531 (2004), DOI `10.1109/TPWRD.2003.822502`; equations (8)–(10), printed p. 1524; appendix (A.11)–(A.13), printed p. 1531.

**Transcription source.** Original publication PDF. The matrix reduction, all three terminal surface terms, appendix material propagation definition, and the complete Bessel auxiliaries were checked against rendered page images. The Markdown drops most symbols in (1)–(10) and (A.1)–(A.12) and was not used as mathematical authority.

## Source transcription

The source's appendix auxiliaries are

```math
D=m_1\rho_1FG+m_2\rho_2EH,

E=I_0(x_3)K_1(x_4)+I_1(x_4)K_0(x_3),
\qquad
F=I_1(x_4)K_1(x_3)-I_1(x_3)K_1(x_4),

G=I_0(x_2)K_1(x_1)+I_1(x_1)K_0(x_2),
\qquad
H=I_1(x_2)K_1(x_1)-I_1(x_1)K_1(x_2),

P=I_0(x_1)K_1(x_2)+I_1(x_2)K_0(x_1),
\qquad
Q=I_0(x_2)K_0(x_1)-I_0(x_1)K_0(x_2),

R=I_0(x_4)K_1(x_3)+I_1(x_3)K_0(x_4),
\qquad
S=I_0(x_4)K_0(x_3)-I_0(x_3)K_0(x_4),
\qquad\text{(A.13)}

x_1=m_1a,\quad x_2=m_1b,\quad x_3=m_2b',\quad x_4=m_2c.
```

The boundary relations from which (A.11) is solved are continuity ``E_{z1}(x_2)=E_{z2}(x_3)`` and ``H_{\theta1}(x_2)=H_{\theta2}(x_3)`` in (A.7), with the terminal relation in (A.10).

## Notation map

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_{in},Z_{out},Z_m`` | unchanged | Inner, outer, and transfer surface impedances of the bonded composite | ``\Omega/\mathrm m`` |
| ``z_{ni},z_{no},z_{nm}`` | unchanged | Single-layer inner, outer, and transfer surface terms for medium ``n`` | ``\Omega/\mathrm m`` |
| ``a,b,b',c`` | unchanged | Successive cylindrical radii | m; bonded case ``b=b'`` |
| ``\rho_n,\mu_n,\varepsilon_n,\sigma_n`` | unchanged | Resistivity, permeability, permittivity, conductivity of medium ``n`` | SI |
| ``\eta_n`` | unchanged | Intrinsic admittance ``\sigma_n+j\omega\varepsilon_n`` | ``\mathrm{S/m}`` |
| ``m_n,x_1,\ldots,x_4`` | unchanged | Radial propagation and Bessel arguments | ``m_n`` in ``\mathrm m^{-1}``; ``x_k`` dimensionless |
| ``I_n,K_n`` | unchanged | Modified Bessel functions of order ``n`` | dimensionless |
| ``\Gamma`` | unchanged | Impressed longitudinal propagation constant | ``\mathrm m^{-1}`` |

## Evidence and approximation sources

- Physical geometry and material roles: Fig. 1, printed p. 1523.
- Penetration-depth condition and isolated-layer construction: (1) and following text, printed p. 1524.
- Component terms and terminal reduction: (2)–(10), printed p. 1524.
- Exact Maxwell field assumptions and propagation factor: (A.1)–(A.5), printed p. 1530.
- Boundaries, direct terminal solution, good-conductor reduction, and auxiliaries: (A.7)–(A.13), printed p. 1531.

## Limitations and discrepancies

- The Markdown conversion preserves fragments of (3), (12), and (A.13) but drops the headline (8)–(10) and direct (A.11) equations. This record follows the PDF images.
- The source calls the appendix expression “accurate”; this record describes it as the direct Maxwell/Bessel representation under the printed circular, harmonic, and axial assumptions rather than assigning a universal accuracy label.
