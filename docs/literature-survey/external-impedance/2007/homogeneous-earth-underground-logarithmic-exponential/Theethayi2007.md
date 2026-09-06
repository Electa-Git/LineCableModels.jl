# Theethayi logarithmic-exponential buried-wire impedance

## Identity and source

| Field | Value |
| --- | --- |
| Family | External impedance |
| Geometry | Circular bare conductor radius ``a`` or concentric insulated wire with outer radius ``b``; ``R_{ab}=a`` for bare, ``R_{ab}=b`` for insulated. Burial depth is ``d``. |
| Calculated quantities | Buried bare/insulated-wire earth impedance; source-prescribed parallel-wire mutual substitution |
| Earth structure | Homogeneous conducting dielectric half-space below air; empirical correction represents burial-depth/interface effects. |
| Model and approximation | Empirical modification of the infinite-earth logarithmic approximation attributed to Petrache et al. (2005) in (8). A depth-dependent exponential term is added in (9); the inspected paper does not derive an expansion parameter, retained order, or discarded remainder for it. Its Sunde/Wait comparisons are presented as numerical evidence within the stated geometries. |
| Main source | Nelson Theethayi's 2005 doctoral thesis, explicitly cited as reference [14] by the inspected 2007 paper; the record year identifies this later published witness |
| Citation key(s) | Primary publication: `:Theethayi2007`; thesis source: `:Theethayi2005` |
| Evidence status | Author's later explicit restatement and original thesis equation (7.10) checked against page images; the relevant thesis sections contain no additional distinct formula in the requested families. |

**Description.** Logarithmic-exponential approximation to the earth-return impedance of a horizontal bare or insulated buried wire in a homogeneous half-space, with an empirical depth-dependent addition to an infinite-earth logarithmic expression.

**Assumptions.**

| Field | Treatment | Evidence |
| --- | --- | --- |
| Impressed longitudinal propagation constant ``Γ`` | No independent imposed longitudinal constant appears in (9). The source adopts a TEM/transmission-line description; a prescribed longitudinal phase constant for deriving (9) is not stated. | Stated — section 2, p. 752; equation-implied — (9), p. 754. |
| Air propagation constant ``γ_air`` | Air is drawn with ``\mu_0,\epsilon_0``. No air bulk propagation constant is retained in (9); air electrical size enters the surrounding TL-validity discussion. | Stated/equation-implied — Fig. 1, p. 752; (9) and opening of p. 754. |
| Earth propagation constant ``γ_earth`` | ``\gamma_g=\sqrt{j\omega\mu_0(\sigma_g+j\omega\epsilon_g)}``; no separate square-root branch stated. | Stated — definition after (5), p. 753. |
| Earth permittivity and displacement current | Retained through ``j\omega\epsilon_g``; ``\epsilon_g=\epsilon_{rg}\epsilon_0``. Conductivity and permittivity are separate scalar inputs, without a dielectric-loss model in this expression. | Stated — section 2, p. 752; definition after (5), p. 753. |
| Range of validity | Source TL condition ``d\sqrt{\epsilon_0\mu_0\omega^2}\ll1``, with approximately 5 MHz discussed for depths 0.5–1 m. Comparisons extend to 10 MHz for radius 2 cm, depth 0.5 m, ``\epsilon_{rg}=10``, conductivities 0.1, 1, and 10 mS/m; these are test cases, not a universal accuracy guarantee. | Stated — opening and final paragraphs of p. 754; Figs. 3a–3b. |
| Earth permeability ``μ_earth`` | Fixed ``\mu_0``. | Stated — Fig. 1 and section 2, p. 752. |
| Arrangement | Underground: self earth impedance of bare and insulated wires. For parallel-wire mutual impedance, the source says to replace ``R_{ab}`` by horizontal distance and ``d`` by average depth. It does not print a separate unequal-depth correction. | Stated — opening of section 3, p. 753. |
| Earth structure | Homogeneous conducting dielectric half-space below air; empirical correction represents burial-depth/interface effects. | Stated — Fig. 1, p. 752; paragraph introducing (9), p. 754. |
| Conductor and insulation geometry | Circular bare conductor radius ``a`` or concentric insulated wire with outer radius ``b``; ``R_{ab}=a`` for bare, ``R_{ab}=b`` for insulated. Burial depth is ``d``. | Stated — section 2 and Fig. 1, p. 752; section 3, p. 753. |
| Constitutive and field assumptions | Linear scalar ``\sigma_g,\epsilon_g,\mu_0``; source TEM/transmission-line model. Conductor internal impedance is omitted from the analyzed total; insulation inductive contribution remains separate. | Stated — section 2, p. 752. |
| Conventions | ``d`` is positive burial depth; ``x`` is longitudinal position in the line equations. Frequency-domain equations use ``j\omega`` but no separate time-factor prescription is stated. ``Z_g`` is per unit length. The absolute value in ``e^{-2d|\gamma_g|}`` is part of the published formula. | Stated/equation-implied — Figs. 1–2 and (1)–(3), pp. 752–753; (9), p. 754. |

**Expression.** The author's later restatement of the logarithmic-exponential formula, equation (9).

```math
Z_g^{\mathrm{LOGEXP}}=\frac{j\omega\mu_0}{2\pi}
\left\{\ln\left(\frac{1+\gamma_gR_{ab}}{\gamma_gR_{ab}}\right)
+\left[\frac{2e^{-2d|\gamma_g|}}{4+\gamma_g^2R_{ab}^2}\right]\right\},
\qquad\text{(9)}

\gamma_g=\sqrt{j\omega\mu_0(\sigma_g+j\omega\epsilon_g)},\qquad
\epsilon_g=\epsilon_{rg}\epsilon_0.
```

``Z_g`` has units ``\Omega/\mathrm m``; ``\gamma_g`` has units ``\mathrm m^{-1}``. ``R_{ab}=a`` for a bare conductor and ``R_{ab}=b`` for an insulated one, with ``a,b,d`` in meters. For mutual impedance, preserve the source's textual substitutions: horizontal wire separation for ``R_{ab}``, average burial depth for ``d``. This record does not supply a reconstructed direct-distance term for unequal depths.

**Approximation.** Empirical modification of the infinite-earth logarithmic approximation attributed to Petrache et al. (2005) in (8). A depth-dependent exponential term is added in (9); the inspected paper does not derive an expansion parameter, retained order, or discarded remainder for it. Its Sunde/Wait comparisons are presented as numerical evidence within the stated geometries.

**Limitations.** Only the earth impedance is computed. ``|\gamma_g|`` appears in the exponential, whereas the denominator retains the complex ``\gamma_g^2``. No explicit branch for the complex square root/logarithm or claimed extension to arbitrary stratification is supplied. The paper's simplified mutual substitution is retained as stated, with its unequal-depth applicability unresolved.

**Reference.** [Theethayi2005](@cite).  [Theethayi2007](@cite), (8)–(9), p. 754; definitions/substitutions, pp. 752–753. Original attribution given in reference [14]: N. Theethayi, *Electromagnetic Interference in Distributed Outdoor Electrical Systems, with an Emphasis on Lightning Interaction with Electrified Railway Network*, Uppsala University doctoral thesis (2005), ISBN 91-554-6301-0.

**Transcription source.** Author's later explicit restatement, verified against the 2007 PDF images. The [Uppsala thesis](https://www.diva-portal.org/smash/get/diva2%3A166746/FULLTEXT01.pdf) equation (7.10), printed p. 133, was also image-checked. The identity table retains the 2007 witness's form and assumptions.

## Source transcription

The source first prints the Petrache-attributed parent

```math
Z_g^{\mathrm{LOG}}=\frac{j\omega\mu_0}{2\pi}
\ln\left(\frac{1+\gamma_g\cdot R_{ab}}{\gamma_g\cdot R_{ab}}\right).
\qquad\text{(8, secondary parent witness)}
```

The immediately following equation (9) and definitions are reproduced without renaming in the formula section. For an insulated wire the source places its ground impedance ``Z_{gi}`` in series with ``j\omega L`` in (2a), where ``L=(\mu_0/2\pi)\ln(b/a)`` from (3a). That insulation expression is attributed to earlier references [5,10], not a new 2007 contribution.

The 2005 original prints two separate prefactors:

```math
Z_{gbbi}^{\mathrm{Log-Exp}}
=\frac{j\omega\cdot\mu_0}{2\pi}\cdot
\ln\left(\frac{1+\gamma_g\cdot R_{ab}}{\gamma_g\cdot R_{ab}}\right)
+\frac{j\omega\cdot\mu_0}{2\pi}\cdot
\frac{2e^{-2d\cdot|\gamma_g|}}{4+R_{ab}^2\cdot\gamma_g^2}.
\qquad\text{(thesis 7.10)}
```

The same page gives the unnumbered limit

```math
\omega\to\infty:\qquad
Z_{gbbi}^{\mathrm{Log-Exp}}\to
\frac{1}{2\pi\cdot R_{ab}}\sqrt{\frac{\mu_0}{\epsilon_g}}.
```

The thesis explicitly describes combining its Petrache logarithm (7.9) and Saad depth term (7.5), followed by sensitivity analyses. This resolves the empirical-parent attribution. The high-frequency limit is source-stated; it does not extend TL validity to infinite frequency. Thesis ``Z_{gbbi}`` denotes the ground impedance for the bare/insulated cases represented by ``Z_g`` in the 2007 witness.

## Notation map

Notation is unchanged.

| Source symbol | Display symbol | Physical meaning | Units/convention |
| --- | --- | --- | --- |
| ``Z_g^{\mathrm{LOGEXP}},Z_g^{\mathrm{LOG}}`` | unchanged | Empirical and parent ground impedance | ``\Omega/\mathrm m`` |
| ``Z_{gb},Z_{gi}`` | unchanged | Ground impedance for bare and insulated cases | ``\Omega/\mathrm m`` |
| Thesis ``Z_{gbbi}^{\mathrm{Log-Exp}}`` | unchanged in original witness; corresponds to 2007 ``Z_g^{\mathrm{LOGEXP}}`` | Original thesis's bare/insulated ground impedance | ``\Omega/\mathrm m``; both printed decompositions retained |
| ``\gamma_g`` | unchanged | Bulk-earth propagation constant | inverse m; not longitudinal line ``Γ`` |
| ``\sigma_g,\epsilon_g,\epsilon_{rg}`` | unchanged | Earth conductivity, absolute and relative permittivity | S/m, F/m, dimensionless |
| ``\mu_0,\epsilon_0`` | unchanged | Vacuum permeability and permittivity | H/m, F/m |
| ``a,b,R_{ab}`` | unchanged | Conductor radius, insulation outer radius, selected earth-boundary radius | m |
| ``d`` | unchanged | Burial depth; average depth for the stated mutual substitution | m, positive downward |
| ``j,\omega`` | unchanged | Imaginary unit and angular frequency | dimensionless, rad/s |
| ``L`` | unchanged | Insulation inductance per length | H/m |

## Evidence and approximation sources

The 2007 authors identify the extra term as empirical and credit Theethayi's 2005 thesis. This record therefore preserves 2007 as the inspected publication year without asserting a 2007 priority date. The source's (4) is a Saad approximation with conduction-only ``\gamma'_g``; that is not substituted into the displacement-current-retaining (9). The unmodified Saad original is recorded separately in this corpus.

## Limitations and discrepancies

- The PDF text extraction loses radicals and absolute-value bars in several comparison/validity formulas; the page image controls the displayed ``|\gamma_g|`` and the TL electrical-size condition.
- The proposed term uses the modulus only in its exponential. Replacing it with complex ``\gamma_g`` would produce a different expression.
- Actual existing LCM ID year 2007 identifies the later article, while that article credits the 2005 thesis. No ID was renamed.
- The mutual substitution contains no vertical separation difference; its scope should not be silently expanded from the source's short textual prescription.
