# Buried-wire ground admittance in Theethayi et al.

## Identification and source

| Field | Value |
| --- | --- |
| Family | External admittance |
| Geometry | Bare circular radius ``a`` or concentric insulated wire with outer radius ``b`` at depth ``d``. Ground impedance selects ``R_{ab}=a`` or ``b`` respectively. |
| Calculated quantities | Scalar ground admittance of a bare wire and of an insulated wire; source-supplied insulation series assembly |
| Earth structure | Homogeneous conducting dielectric half-space below air. |
| Model and approximation | A source-attributed scalar relation in the selected TL description. When ``Z_g`` is approximated by (9), its empirical approximation carries into ``Y_g``. |
| Main source | E. F. Vance (1978), as reproduced by later cable-parameter publications |
| Citation key(s) | Attributed source: `:Vance1978`; equation witnesses: `:Theethayi2007`, `:Zhang2017`, `:Guneri2018` |
| Evidence status | The relation is independently printed in three accessible publications |

**Description.** Scalar ground admittance obtained from the corresponding ground impedance and bulk-earth propagation constant for a bare or concentrically insulated horizontal buried wire.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | No imposed longitudinal constant in (10); the line's voltage/current propagation is governed by the source's TL equations. It must not be identified with bulk ``\gamma_g``. | Stated/equation-implied — section 2, p. 752; (1)–(2), p. 753; (10), p. 755. |
| Air propagation constant ``γ_air`` | Air ``\mu_0,\varepsilon_0`` are shown in Fig. 1; no air propagation constant appears in the scalar admittance relation. | Stated/equation-implied — Fig. 1 and (10). |
| Earth propagation constant ``γ_earth`` | ``\gamma_g=\sqrt{j\omega\mu_0(\sigma_g+j\omega\varepsilon_g)}``. | Stated — definition following (5), p. 753. |
| Earth permittivity and displacement current | Retained in ``\gamma_g``; ``\varepsilon_g=\varepsilon_{rg}\varepsilon_0``. | Stated — section 2, p. 752; definition after (5). |
| Range of validity | Source TEM/TL scope and the restrictions of the selected ground-impedance model. The surrounding discussion gives ``d\sqrt{\varepsilon_0\mu_0\omega^2}\ll1`` for the quasi-static/TL treatment; it is not a separate exactness proof of (10). | Stated — p. 754 opening; (10), p. 755. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0``. | Stated — section 2 and Fig. 1. |
| Arrangement | Underground scalar bare and insulated-wire cases. No explicit mutual-admittance matrix inversion is printed with (10). | Stated — Figs. 1–2 and (10a)–(10b). |
| Earth structure | Homogeneous conducting dielectric half-space below air. | Stated — Fig. 1, p. 752. |
| Conductor and insulation geometry | Bare circular radius ``a`` or concentric insulated wire with outer radius ``b`` at depth ``d``. Ground impedance selects ``R_{ab}=a`` or ``b`` respectively. | Stated — section 2 and section 3, pp. 752–753. |
| Constitutive and field assumptions | Linear scalar earth properties and TEM/TL model. Insulation capacitance is lossless ``C=2\pi\varepsilon_{in}/\ln(b/a)`` in the printed assembly. | Stated/equation-implied — section 2 and (3b). |
| Conventions | Per-length ``Y_{gb},Y_{gi}`` in S/m, ``Z_{gb},Z_{gi}`` in Ω/m. Depth positive downward, longitudinal coordinate ``x``. The paper prints positive right-hand sides in (1)–(2); that sign is retained below. Separate time-factor convention not stated. | Stated/equation-implied — Figs. 1–2, (1)–(3), and (10). |

**Expression.** The secondary witness's scalar ground relations, equations (10a)–(10b).

```math
Y_{gb}=\frac{\gamma_g^2}{Z_{gb}},
\qquad\text{(10a)}

Y_{gi}=\frac{\gamma_g^2}{Z_{gi}},
\qquad\text{(10b)}

\gamma_g=\sqrt{j\omega\mu_0(\sigma_g+j\omega\varepsilon_g)}.
```

``b`` and ``i`` in the subscripts distinguish bare and insulated cases. ``Z_{gb}`` and ``Z_{gi}`` are the corresponding earth impedance components, not the full conductor-plus-insulation series impedance. For the logarithmic-exponential option presented in the same source,

```math
Z_g^{\mathrm{LOGEXP}}=\frac{j\omega\mu_0}{2\pi}
\left\{\ln\left(\frac{1+\gamma_gR_{ab}}{\gamma_gR_{ab}}\right)
+\left[\frac{2e^{-2d|\gamma_g|}}{4+\gamma_g^2R_{ab}^2}\right]\right\},
\quad R_{ab}=a\ \text{(bare)},\quad R_{ab}=b\ \text{(insulated)}.
\qquad\text{(9)}
```

The paper also considers other ground-impedance models; (10) is not printed as exclusive to (9). ``d`` is burial depth, ``a`` conductor radius, and ``b`` outer insulation radius. For the insulated wire the source-prescribed assembly is retained in its actual line-equation form:

```math
\frac{dI(x,j\omega)}{dx}
=j\omega\left(\frac{CY_{gi}}{j\omega C+Y_{gi}}\right)V(x,j\omega),
\qquad\text{(2b)}

C=\frac{2\pi\varepsilon_{in}}{\ln(b/a)},\qquad
\varepsilon_{in}=\varepsilon_{rin}\varepsilon_0.
\qquad\text{(3b and section 2 definition)}
```

**Approximation.** A source-attributed scalar relation in the selected TL description. When ``Z_g`` is approximated by (9), its empirical approximation carries into ``Y_g``. The inspected paper supplies no separate derivation or asymptotic order for the Vance-attributed relation (10); original-source verification is unresolved.

**Limitations.** This record does not turn the scalar quotient into elementwise multiconductor admittances, a matrix inverse, or potential coefficients. None of those operations is printed in (10). The ground admittance and full insulated-wire shunt contribution have different meanings, as shown in (2b).

**Reference.** [Vance1978](@cite).  [Theethayi2007](@cite), equations (10a)–(10b), p. 755; assembly (2b)–(3b), p. 753. The inspected source cites E. F. Vance, *Coupling to Shielded Cables* (Wiley, 1978), reference [7], for (10); the equation witness does not supply a book page locator.

**Transcription source.** Secondary-source PDF-image verification only, pp. 752–755. The numerator ``\gamma_g^2``, scalar denominators, insulation assembly, and selected empirical impedance dependency were checked visually. The original Vance book was not inspected.

## Source transcription

The formula section preserves (10a)–(10b), with (9) identified as an optional impedance model, and (2b)–(3b) as source-prescribed assembly. For the bare wire the corresponding printed equation is

```math
\frac{dI(x,j\omega)}{dx}=Y_{gb}\,V(x,j\omega).
\qquad\text{(1b)}
```

No source equation declares ``\mathbf Y_g=\gamma_g^2\mathbf Z_g^{-1}`` or ``P_{ij}=j\omega Z_{ij}/\gamma_g^2`` in the inspected pages; those existing LCM expressions are therefore not transcribed as if printed by this paper.

## Notation map

Notation is unchanged.

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Y_{gb},Y_{gi}`` | unchanged | Bare/insulated ground admittance components | S/m |
| ``Z_{gb},Z_{gi},Z_g^{\mathrm{LOGEXP}}`` | unchanged | Corresponding ground impedance components/model | Ω/m |
| ``\gamma_g`` | unchanged | Bulk-earth propagation constant | inverse m |
| ``\sigma_g,\varepsilon_g,\mu_0`` | unchanged | Earth conductivity, permittivity, permeability | S/m, F/m, H/m |
| ``a,b,R_{ab},d`` | unchanged | Conductor/insulation radii, selected radius, burial depth | m |
| ``C,\varepsilon_{in},\varepsilon_{rin},\varepsilon_0`` | unchanged | Insulation capacitance per length, absolute/relative insulation permittivity, vacuum permittivity | F/m, F/m, dimensionless, F/m |
| ``I,V,x,j,\omega`` | unchanged | Line current, voltage, longitudinal position, imaginary unit, angular frequency | A, V, m, dimensionless, rad/s |

## Evidence and approximation sources

Section 2 distinguishes ground admittance from the insulation capacitance and gives the series connection. Equation (10) explicitly cites Vance [7]. Equation (9) explicitly credits Theethayi's 2005 thesis; its later-author witness and depth approximation are documented in the [impedance record](../../../external-impedance/2007/homogeneous-earth-underground-logarithmic-exponential/Theethayi2007.md).

## Limitations and discrepancies

- Original Vance equation locator is missing; verification against the 2007 page is not verification against the original book.
- The 2005 thesis also attributes ``Y_{gbbi}=\gamma_g^2/Z_g`` to Vance and calls the relation approximate: thesis (7.11), printed p. 133. This is a secondary witness, not the Vance original.
- The scalar source statement is narrower than the existing LCM matrix assembly. This is a documentation distinction, not a proposed implementation change.
- The source's line-equation signs are retained as printed despite the familiar alternative convention.
- The branch of ``\gamma_g`` and complex logarithm is not separately stated in the inspected definitions.
