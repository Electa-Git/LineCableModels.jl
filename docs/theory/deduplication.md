# Formula equivalence and selection

The [implementation table](assimilation.tsv) distinguishes a new implementation,
a checked existing implementation, and an equivalent equation witness. Equivalent
records remain available for source comparison; they are not additional
physical formulas or additional engine registrations. Each equivalence below
has an independent numerical comparison in the literature tests.

## Equivalent equation witnesses

| Source record | Retained formula | Equality and limits |
| --- | --- | --- |
| [Wait1978](internal-impedance/1978/propagation-dependent-solid-core/Wait1978.md#identification-and-source) | `InternalImpedance.Formula(:Schelkunoff1934)` | Only the good-conductor solid-core reduction (3b): ``j k_w=\sqrt{j\omega\mu_w\sigma_w}`` gives the same scaled Bessel ratio. The propagation-dependent, displacement-retaining parent (3a) is not implemented by that selection. |
| [Ametani2014](external-impedance/2014/homogeneous-earth-overhead-air-referenced-displacement-integral/Ametani2014.md#identification-and-source) | `EarthImpedance.Formula(:Wise1934)` | The nonmagnetic restriction of equations (11)–(13), including compact equation (12), gives the same air-referenced overhead integral. The general magnetic extension is not identified with the unit-permeability Wise1934 selection. Direct/image distances and the thin-wire self prescription must be retained. |
| [Xue2018a](external-impedance/2018/quasi-tem-overhead-mom-so-green-function/Xue2018a.md#identification-and-source) | `EarthImpedance.Formula(:Wise1934)` | The phase-domain series term (17)–(20), with ``k_a^2-k_e^2=\gamma_e^2-\gamma_a^2`` and nonmagnetic earth, is the same overhead integral. This does not implement the paper's conductor MoM–SO solve or its full modal parent. |
| [Uribe2008](external-impedance/2008/pollaczek-mixed-algorithmic-evaluation/Uribe2008.md#identification-and-source) | `EarthImpedance.Formula(:Pollaczek1926)`, mixed pair | Positive-domain normalization and rationalization transform the same mixed integral. The engine evaluates the infinite integral adaptively; it does not reproduce the paper's empirical finite-tail truncation. |
| [Nguyen1998](external-impedance/1998/underground-earth-return-direct-quadrature/Nguyen1998.md#identification-and-source) | `EarthImpedance.Formula(:Pollaczek1926)`, underground pair | Substituting ``\lambda=\sqrt{\alpha}u`` and rationalizing the denominator gives the same underground integral. The source's finite trapezoidal sum is an alternative numerical evaluator, not a new physical kernel. |

The earlier equations retain their attribution: [Schelkunoff1934](@cite),
[Wise1934](@cite), and [Pollaczek1926](@cite). The complementary derivations are
[Wait1978](@cite), [Ametani2014](@cite), [Xue2018a](@cite),
[Uribe2008](@cite), and [Nguyen1998](@cite).

## Registry consolidation

The former `Gary1976` and description-only `DeriSemlyen1981` registrations
describe the homogeneous complex return plane. They now resolve to
`Dubanton1969`, matching the source record. [Deri1981](@cite) remains an
English equation-complete reference, with [Wait1969](@cite) as a complementary
image derivation. The independently selectable `Deri1981` registration
implements the different multilayer effective-depth construction.

The homogeneous self term uses ``2(h+p)/r``; the mutual term uses the
complex image distance. The finite wire radius is not inserted as a horizontal
self separation in the image numerator.

## Displacement-current image witnesses

The overhead series image in `EarthImpedance.Formula(:Pettersson1994)`
is shared by the following later records. Its air-referenced depth differs
from the conduction-only depth of Dubanton1969. The `Ametani2014`
impedance selector remains a compatibility alias with `lossless_air=true`;
its duplicate numerical registration has been removed.

| Source record | Retained formula | Equality and limits |
| --- | --- | --- |
| [Maaouni2001](external-impedance/2001/two-half-space-qtem-image/Maaouni2001.md#identification-and-source) | `EarthImpedance.Formula(:Pettersson1994)` | With lossless air and nonmagnetic media, source (12) gives the same shifted image distance after conversion from `e^-jωt`. Only the qTEM overhead image is covered, not the modal determinant (1). |
| [Ametani2014](external-impedance/2014/complex-depth-displacement-current-image/Ametani2014.md#identification-and-source) | `EarthImpedance.Formula(:Pettersson1994; lossless_air=true)` | Equations (17)–(18) use the same `γ_earth²−γ_air²` image depth. Radius enters the self direct distance, while the image correction uses zero lateral self separation. |

`EarthAdmittance.Formula(:Maaouni2001)` implements the distinct
[closed potential approximation](external-admittance/2001/two-half-space-potential/Maaouni2001.md#identification-and-source),
including the source helper (25) and expression (26). Its unapproximated
qTEM parent is Wise1948, but the analytical approximation is not identical
to that integral. Both transformed separation arguments must have positive
real parts. The source concerns overhead wires, not buried or mixed pairs.
The potential coefficients are assembled before matrix inversion.

## Buried-wire potential normalization

[Zhang2017](external-admittance/2017/buried-ground-admittance-asymptotic-extraction/Zhang2017.md#identification-and-source)
has its own `EarthAdmittance.Formula(:Zhang2017)` selection. Its separate
lossless axial references give `u_air=λ` and
`u_earth²=λ²+jωμ0σ`, whereas the interface weights retain displacement
current. This differs from the default Papadopoulos2010b prescription.

The right-hand side of source (11) is returned as a potential coefficient.
The default integral and the optional `evaluation=:asymptotic_tail`
are independently compared with (11)–(16) and (A.6)–(A.8), respectively.
The latter includes the dimensional threshold powers in its exponential
integrals; its finite interval uses adaptive quadrature.

The printed radius is the metal-core radius, not the jacket radius.
Standard matrix assembly samples the outer earth interface;
`source_radius=r` explicitly reproduces the source's core sampling.
Normal scalar evaluation rejects the negative ground conductance obtained
from this prescribed kernel in some high-frequency cases.

The separate logarithmic approximation (17) is
`EarthImpedance.Formula(:Petrache2005)`. Equation (18) uses the existing
Vance1978 scalar conversion with that impedance as its selected dependency.
No extra impedance or conversion registration is introduced.


## Distinct approximations within one publication

`EarthImpedance.Formula(:Saad1996)` selects the main Bessel/exponential
approximation. `EarthImpedance.Formula(:Saad1996; approximation=:small_argument)`
selects equations (31)–(32). These are different approximations with different
validity ranges, not duplicated bibliography entries. Both source tables are
included in the registration's documentation.

`EarthImpedance.Formula(:Uribe2008; approximation=:ccitt)` and
`approximation=:wedepohl` supply the two
[mixed CCITT](external-impedance/2008/mixed-ccitt-recommended-approximation/Uribe2008.md#identification-and-source)
and [mixed Wedepohl](external-impedance/2008/mixed-wedepohl-low-frequency-approximation/Uribe2008.md#identification-and-source)
restatements. Numerical use follows dimensional (6c) and (6e), with the
coordinates and Bessel Euler factor specified in the records. The normalized
(6d) and (6f) do not create extra registrations. The logarithmic calculation
is shared with Wedepohl1973, but the mixed height coefficient is distinct.
The selected same-medium coefficients remain the existing integrals.

`EarthImpedance.Formula(:Ametani2009; approximation=:power_frequency)`
implements the [leading image reduction](external-impedance/2009/homogeneous-earth-mixed-power-frequency/Ametani2009.md#identification-and-source).
The full frequency scale and coordinate distances follow parent (27).
The default `approximation=:image` retains that
[exponential image](external-impedance/2009/homogeneous-earth-mixed-exponential-image/Ametani2009.md#identification-and-source).
Both are compared with independent source expressions and complete
mixed two-wire matrix assembly. The reduction is not uniformly valid
outside its small propagation-distance limit.


`EarthAdmittance.Formula(:Wise1948)` retains the parent integral by
default. Its [analytical alternatives](external-admittance/1948/homogeneous-earth-overhead-potential-coefficient-approximation/Wise1948.md#identification-and-source)
are `approximation=:rational` for (9)–(10), `:coarse` for the first
unnumbered formula, and `:small_g` for the final expression on p. 371.
The two coarse forms require zero horizontal separation; the last also
requires `H*abs(beta)` much smaller than one. Their restrictions do not
constitute a uniform error bound.

All forms use the source normalization `C=2(M+jN)` and
`P=[log(D/d)+C]/(2πε_air)`. Source and independent rational-integral
checks include separation angles crossing the principal exponential-integral
cut. Self sampling uses zero horizontal displacement and the wire radius
in the direct distance only. Full two-wire potential matrices are checked
before and after inversion.


`EarthAdmittance.Formula(:DiLorenzo2023)` now supplies the
[seabed potential coefficient](external-admittance/2023/three-medium-seabed-potential-matrix/DiLorenzo2023.md#identification-and-source).
Its electric correction is normalized by the magnetic and electric
boundary determinants obtained from source (27). The printed numerator
(34) is preserved in the record. An independent solution of all eight
Hertz-potential boundary conditions confirms the normalized correction,
including independent layer permeabilities.

The zero-water-depth limit agrees with the permeability-general
MartinsBritto2024 homogeneous coefficient at zero longitudinal input.
For nonmagnetic media it also agrees with the Papadopoulos2010b reduction.
Removing the water–seabed contrast recovers a single half-space at the
original air interface. The two-cable, four-terminal assembly reproduces
source (36)–(40), with radial insulation terms supplied by the existing
local matrix construction.


## Equivalent-soil formulas

`EarthProps.EHEM.Formula(:MartinsBritto2020)` already supplies the real
equivalent-conductivity recursion. Each layer's diffusion constant retains its
permeability; the published numerical examples use nonmagnetic soil. The
result is composed with a homogeneous impedance formula, normally Carson1926.
It does not itself return an impedance matrix.

`EarthProps.EHEM.Formula(:Xue2021)` already supplies the air-referenced
equivalent-propagation recursion and material reconstruction. Compose it with
Wise1934 for the surveyed overhead impedance. For the separately printed
[EHEM potential formula](external-admittance/2021/n-layer-equivalent-propagation-overhead/Xue2021.md#identification-and-source),
select `EarthAdmittance.Formula(:Xue2021)`. It recovers the transverse
constant from the equivalent material before evaluating equation (8).
The separate exact four-layer field expressions are implemented as described below.

The potential denominator uses the transverse constant squared divided by
the air constant squared. Substituting the reconstructed bulk constant would
add one to that ratio and produce a different formula. The retained
underground Xue book approximation is a separate branch, not the overhead
EHEM paper's expression.

## Radial dielectric relations

| Source record | Retained formula | Equality and limits |
| --- | --- | --- |
| [Pawlik2020](insulation-admittance/2020/thin-wire-lossy-coating/Pawlik2020.md#identification-and-source) | `InsulationAdmittance.Formula(:Ametani2004)` and the radial assembler | Equation (13) is the same conducting-dielectric coaxial relation. Equation (11) follows by adding the coating and exterior potential coefficients before the scalar inverse. No elementwise matrix inverse is introduced. |
| [Ghosh2019](insulation-admittance/2019/multiple-semiconducting-screens-series-coaxial/Ghosh2019.md#identification-and-source) | Existing radial series assembler with the selected insulation and semicon relations | The reciprocal three-annulus network (35) is identical. This mapping covers the network algebra, not the printed single-layer normalization in (36)–(37). |
| [Ghosh2022](insulation-admittance/2022/n-semiconductor-screens-radial-network/Ghosh2022.md#identification-and-source) | Existing radial series assembler with the selected insulation and semicon relations | The arbitrary-screen grouping (17) is represented by conductor-owned ranges of physical dielectric annuli. Loop entries use the specified series groups; the full nodal matrix is assembled separately. This mapping does not assert equality of the printed single-layer coefficients. |

The physical SI annular coefficient is
`y = 2π(σ + jωε)/log(r_out/r_in)`.
The printed Ghosh expressions omit `2π` and combine a relative permittivity
with a term having absolute-permittivity units. Their network equations
remain usable with independently defined physical layer admittances;
the unnormalized layer expressions are not added as alternative physical laws.

The former `InsulationAdmittance.Formula(:Gustavsen2013)` selection is an
alias for the earlier `Ametani1980` lossless relation. Its complete potential
matrix and nodal admittance are tested against equations (21)–(24).
`Weeks1984` remains distinct: its selected loss tangent represents the
complete dielectric loss, rather than adding a second conduction term to it.

## Circular finite-wall pipe

`PipeImpedance.Formula(:DaSilva2006)` supplies the finite-wall Method 3
cavity and inner-surface coefficient. Its source table is copied exactly
from the [2006 conference-paper record](internal-impedance/2006/finite-pipe-without-core-proximity-auxiliary-model/DaSilva2006.md#identification-and-source).

The evaluator excludes the individual core-skin contribution already
provided by the coaxial conductor formula. Pipe outer-surface, transfer,
jacket, and exterior-earth contributions must be assembled separately.
Explicit circular common-pipe designs now assemble these contributions.
The pipe-return transformation cancels the outer-surface, transfer, jacket,
and exterior-earth contributions from the interior loop matrix as required.
One enclosing pipe supplies one exterior-earth representative, irrespective
of its number of internal coaxial units. The tested construction uses a
homogeneous nonmagnetic cavity and contiguous radial units on distinct axes.

`PipeImpedance.Formula(:Hoidalen2013)` selects the
[finite-wall low-frequency surface and cavity terms](internal-impedance/2013/finite-pipe-low-frequency-surface-and-loop-terms/Hoidalen2013.md#identification-and-source).
Its inner, outer, and transfer impedances satisfy the source connection
limit (22). For magnetic pipes the harmonic limit is taken from the full
finite-wall parent (9)–(11), retaining the outer radius. The printed
low-frequency (13) agrees for nonmagnetic pipes or infinite outer radius,
but omits that finite-thickness dependence otherwise.

`PipeImpedance.Formula(:Yang2001)` selects the different
[finite-surface/infinite-harmonic hybrid](internal-impedance/2001/finite-pipe-hybrid-thickness/Yang2001.md#identification-and-source).
Its full finite-wall surface terms do not make its positive-order
eddy-current terms finite-wall expressions. This distinction matters for
thin magnetic pipes at low frequency. Both new selections use the same
common-pipe matrix assembly as DaSilva2006.

| Source record | Retained formula | Equality and limits |
| --- | --- | --- |
| [Fortin2005](internal-impedance/2005/finite-pipe-eddy-current/Fortin2005.md#identification-and-source) | `PipeImpedance.Formula(:DaSilva2006)` | The four finite-wall boundary equations (6)–(7), solved independently for each angular order, give the same pipe harmonic response. The positive pipe-return inner-surface response is shared. No discretized core or additional core proximity is implied. |
| [DeSilva2019](internal-impedance/2019/finite-pipe-appendix-correction/DeSilva2019.md#identification-and-source) | `PipeImpedance.Formula(:Yang2001)` | The finite-wall inner surface (A6) and infinite-wall `K_(n-1)/K_n` harmonics (A2) repeat the hybrid. The assembler supplies (A3)–(A5); they are not additional pipe material laws. The earlier Markdown `I_(n-1)` in (A2) was a transcription error, corrected against the original page. |

## Mixed and seabed kernels

| Source record | Retained formula | Equality and limits |
| --- | --- | --- |
| [MartinsBritto2024](external-impedance/2024/homogeneous-earth-mixed-generalized/MartinsBritto2024.md#identification-and-source) | `EarthImpedance.Formula(:Pawlik2018)` | The mixed series kernels are identical after changing propagation variables. For Pawlik's axial attenuation constant `Γ_p`, the engine spectral input is `k_x=jΓ_p`. The existing `MartinsBritto2024` impedance selector remains an alias. |
| [Dawalibi1989](external-impedance/1989/homogeneous-earth-mixed-low-frequency/Dawalibi1989.md#identification-and-source) | `EarthImpedance.Formula(:Pawlik2018)` with explicit `k_x=0` | The permeability-weighted mixed integral is the zero-longitudinal-propagation restriction. Use `γ_i²=jωμ_i(σ_i+jωε_i)`; the manual's unsquared left side in its propagation definition is dimensionally inconsistent. No mixed admittance is attributed to this source. |
| [Pawlik2018](external-admittance/2018/two-half-space-mixed-full-spectrum/Pawlik2018.md#identification-and-source) | `EarthAdmittance.Formula(:MartinsBritto2024)` | Multiplying the source pair integral by `jω/[π(σ_1+jωε_1)]` gives the implemented mixed potential coefficient: `P_pair=jω/Y_Pawlik`. This normalization does not identify `Y_Pawlik` with an element of the final nodal admittance matrix. The complete matrix is assembled and inverted using the 2024 formulation. |
| [DiLorenzo2023](external-impedance/2023/three-medium-seabed-return-integral/DiLorenzo2023.md#identification-and-source) | `EarthImpedance.Formula(:Tsiamitros2008)` | Equations (29)–(30) are the three-medium, bottom-layer specialization with explicit `Γ=0`, full material admittivity, and physical air/sea/seabed thicknesses. Source depths map as `h_i=h_s+d_i`; only the impedance is equivalent here. The separate potential correction is registered as `EarthAdmittance.Formula(:DiLorenzo2023)`. |

The mixed comparisons use the same decaying/outgoing square-root branch.
The propagating lossless-air boundary value must not change because an
algebraic subtraction produces a negative signed zero. The numerical tests
include both propagation prescriptions, unequal permeabilities, and reversed
source/observation pairs.

## Insulated-wire external witnesses

| Source record | Retained formula | Equality and limits |
| --- | --- | --- |
| [Pawlik2020, external impedance](external-impedance/2020/thin-insulated-wire-two-half-space/Pawlik2020.md#identification-and-source) | `EarthImpedance.Formula(:Pawlik2018)`, same-medium self | Equations (14), (16)–(17) give the same external coefficient at the insulation outer radius. The source attenuation input maps as `k_x=jΓ`. No internal metal or insulation term is added twice. |
| [Pawlik2020, external admittance](external-admittance/2020/thin-insulated-wire-two-half-space/Pawlik2020.md#identification-and-source) | `EarthAdmittance.Formula(:MartinsBritto2024)`, same-medium self | Equations (15)–(18) give the returned potential `jω/Y_e`. The two permeability-weighted factors in (18) agree with the factored same-medium kernel. This is a scalar source witness, not a new multiconductor derivation. |

The radius remains the horizontal surface sample in both the cosine and the
image distance for these records. At exactly zero transverse propagation,
the direct-minus-image Bessel term is evaluated by its logarithmic limit.
The lossless outgoing square-root boundary is selected from frequency sign,
not from a signed zero left by arithmetic. The potential kernel is factored
to avoid cancellation at large earth-to-air material ratios. The separate
full modal solve remains outside these source-to-coefficient mappings.

## Complex frequency-dependent soil

`EarthImpedance.Formula(:DeLima2007)` supplies the overhead and buried
records using the source's full complex soil propagation constant. It reuses
the existing Carson and Pollaczek integrals where their evaluations are
identical, but does not use their conduction-only material restriction.

The buried evaluator selects the appendix's `K₀` image term, which recovers
the stated conductive limit. Both printed source witnesses remain in the
Markdown record. The joint soil law (7)–(8) maps to the existing
`Portela1999` parameterization as documented in the formula docstring;
the impedance evaluator does not apply that material transformation twice.
No mixed-pair or external-admittance expression is inferred from this source.

## Mixed image and Padé approximations

The existing `EarthImpedance.Formula(:Lucca1994)` agrees with the inspected
[Uribe2008 secondary equation (6a)](external-impedance/2008/mixed-lucca-two-step-image/Uribe2008.md#identification-and-source).
Its source table is reproduced in the registration. The same-medium branches
reuse Carson and Pollaczek; only the mixed image correction is attributed
to Lucca.

`EarthImpedance.Formula(:DeConti2024)` implements the
[Padé approximation](external-impedance/2024/sunde-pade-closed-form/DeConti2024.md#identification-and-source)
as a distinct selectable formula. Algebraic cancellation and a finite
integral of the same Padé rational expression retain low-frequency numerical
accuracy. This evaluator is checked against the printed closed form in
higher precision, not identified with the unapproximated Sunde integral.

`EarthImpedance.Formula(:Lima2012)` supplies the overhead closed form and
the distinct buried asymptote. The former Theodoulidis2015 registration
repeated the overhead special-function identity; its constructor now selects
Lima2012 with conduction-only soil, preserving that material restriction.
For unequal overhead heights, evaluation follows the height sum in the
source's parent integral (1). The buried approximation enforces the stated
horizontal-spacing restriction and shares only its Bessel/image terms with
DeConti2024, not the Padé residual.

## Modal jacket and common-pipe geometry

| Source record | Retained formula | Equality and limits |
| --- | --- | --- |
| [Wait1978](insulation-impedance/1978/thin-jacket-modal-series/Wait1978.md#identification-and-source) | `InsulationImpedance.Formula(:Ametani1980)` and `InsulationAdmittance.Formula(:Ametani2004)` | The printed thin-jacket operator is exactly `Z_magnetic + β²P/(jω)` for nonmagnetic insulation, including complex permittivity through `κ=jωε`. The propagation term belongs to the modal operator, not a second contribution to the primitive series matrix. |
| [Kane1995](insulation-impedance/1995/offset-cores-inside-cylindrical-shield-logarithmic-field-terms/Kane1995.md#identification-and-source) | Shared circular geometry used by `PipeImpedance.Formula(:DaSilva2006)` | The self and mutual logarithms in (8),(15) are identical to its geometric field terms. This equivalence excludes all metal surface and pipe-harmonic terms; it does not identify Kane's field inductance with the complete finite-wall coefficient. |

`PipeAdmittance.Formula(:Kane1995)` implements the
[lossless self and mutual potential relations](insulation-admittance/1995/offset-cores-shield-lossless-capacitance-potential-relations/Kane1995.md#identification-and-source)
for a common homogeneous dielectric. Its source table is copied exactly.
The self coefficient follows the scalar core-to-shield capacitance product;
the mutual coefficient follows (16). Their assembly and complete matrix
inverse are separate engine operations. The selected constitutive relation
supplies any conducting-dielectric extension; that extension is not
attributed to Kane's lossless equations.

## Carson normalization and self geometry

| Source record | Retained formula | Equality and limits |
| --- | --- | --- |
| [Iwamoto1958a](external-impedance/1958/homogeneous-earth-overhead-logarithmic-integral-evaluation/Iwamoto1958a.md#identification-and-source) | `EarthImpedance.Formula(:Carson1926)` | The main MKS parent integrals (2),(3) are identical after scaling the spectral variable by `sqrt(ωμ₀/ρ)` and rationalizing the square-root difference. The appendix's Duhamel operator and graphical evaluator are not implemented as additional physical kernels. |

Carson's self ground correction uses zero lateral separation. The finite
wire radius appears only in `log(2h/r)`, as the source states. The former
radius insertion in the correction integral is removed. The homogeneous
overhead branches of Carson, Sunde and DeLima2007 now share this evaluation,
with each registration retaining its own conductivity/displacement-current
assumptions. Lima2012 uses the same integral where its equivalent Struve
expression would lose numerical accuracy.

## Layered Sunde and Lee relations

The Sunde selection now uses the survey identifier
`EarthImpedance.Formula(:Sunde1949)`. Its default is conduction-only;
`displacement_current=true` retains the former engine's bulk-admittivity
extension. The old `Sunde1968` keyword selector preserves that extension.
The two-layer and recursive evaluations share one boundary relation.

| Source record | Retained formula | Equality and limits |
| --- | --- | --- |
| [Iwamoto1958b](external-impedance/1958/two-layer-earth-overhead-self-integral/Iwamoto1958b.md#identification-and-source) | `EarthImpedance.Formula(:Sunde1949)`, self | The MKS appendix integral (付1) equals the two-layer ground correction after scaling by `sqrt(ωμ₀/ρ₁)`. The complete engine coefficient adds the ideal-ground logarithm. The nonconducting lower half-space is supported; the graphical operator is not a separate physical formula. |

`approximation=:large_spacing` selects Sunde's two-term mutual correction
(4.56), not the integral. The self coefficient retains the integral.
Equal-layer and split-layer tests verify the recursive extension.
Height-scaled quadrature prevents the integral from disappearing when its
spectral contribution is confined near zero.

`EarthImpedance.Formula(:Lee2014)` retains the distinct air-referenced
roots and reuses the existing multilayer boundary solver. Its default
longitudinal input is computed from air, not the first soil layer.
Independent comparisons use the original absolute-depth F/G recurrence
for one to four earth layers. Adaptive quadrature evaluates the parent;
the paper's finite spline and tail approximations are not separate terms.

## Wait circuit reductions

| Source record | Retained formula | Equality and limits |
| --- | --- | --- |
| [Wait1972a](external-impedance/1972/full-wave-overhead-wire-interface/Wait1972a.md#identification-and-source) | `EarthImpedance.Formula(:Carson1926)`, self | The qTEM equations (31),(33) give the same conduction-only correction after rationalizing the root difference. The generalized, propagation-dependent series term and its modal root problem are not implemented by Carson. |

`EarthAdmittance.Formula(:Wait1972a)` supplies the companion qTEM
self potential, not zero external potential and not the full-wave modal
admittance. `EarthImpedance.Formula(:Wait1978)` supplies the buried
self reduction (12); `approximation=:small_argument` selects (13).
Both selections reject unsupported mutual uses.

For these records, processing applies only to the stated reductions.
The implicit full-wave modal solves remain deferred; the engine does not
replace them with a frequency-only coefficient and call it the same model.

## Snow-layer potential equivalence

| Source record | Retained formula | Equality and limits |
| --- | --- | --- |
| [Ametani2001](external-admittance/2000/snow-layer-overhead-potential-coefficient/Ametani2001.md#identification-and-source) | `EarthAdmittance.Formula(:Papadopoulos2009)` | The earlier snow-layer kernel `A₁−λA₂` equals the two-layer potential kernel `F+G` under the prescribed air propagation reference. This mapping covers the mutual potential integral, not a snow material fit or an additional series-impedance term. |

The numerical grouping uses the Japanese `b₁` in both `A₁` factors,
the denominator product explicitly supplied by the English `A₂`,
and the squared interface ratios defined by the bulk constants.
The dimensions and the independently derived potential expression agree
with that selection. The source transcriptions remain unchanged.
Absolute permittivity and conductivity are supplied separately, so a loss
term already included in complex permittivity is not added twice.

## Permeable earth layers

`EarthImpedance.Formula(:Wedepohl1966)` supplies the conduction-only
two-layer expression. Its `displacement_current=true` variant applies the
separately printed air-reference correction. `Moghram1998` supplies
the conduction-only three-layer expression and its one-/two-layer limits.
Both retain independent soil permeabilities.

These formulas and Sunde's nonmagnetic selection share the same
permeability-weighted surface-ratio calculation. Their material and
propagation assumptions remain distinct. Comparisons use the original
hyperbolic expressions, not a second copy of the numerical recursion.

## Alternative Pollaczek evaluations

The three Theodoulidis series are numerical evaluations of one physical
integral. They are implemented under `EarthImpedance.Formula(:Pollaczek1926)`,
using `evaluation=:bessel_product`, `:single_bessel`, or
`:hypergeometric`. `:finite_integral` evaluates their common finite
parent. No additional physical formula identifier is registered.

The single-Bessel evaluator includes the `n=0` term required by (9a).
The hypergeometric evaluator includes the common `exp(-kH)` factor
required by (22). These choices agree with independent quadrature and
the first Bessel series. The source's printed bounds and factors remain
visible in the records, separately from the numerical interpretation.
The `2^n` denominator and dimensionless units of the parent integral
are corrected transcription/description defects, not changes to the
author's equation.

`evaluation=:recursive` supplies Iracheta's moment-series evaluator.
The moment definitions determine its initialization and parity recurrence;
the normalized residual enters the same finite Pollaczek decomposition.
The source's zero-residual cutoff above `|D/p|=30` remains an explicit
approximation. It is tested separately from exact-integral equivalence.
All four source tables are included verbatim in the retained registration.

## De Lima single-wire image and qFW coefficients

`EarthImpedance.Formula(:DeLima2018)` and
`EarthAdmittance.Formula(:DeLima2018)` share one evaluation of
Appendix D, equations (23)–(24). The external series coefficient excludes
the conductor impedance already assembled by the internal family.
The image propagation estimate includes that separately supplied
conductor impedance.

The `approximation=:quasi_full_wave` selection uses the explicitly
supplied image wavenumber, with ``k_x=j\bar\gamma``.
It evaluates the source integrals; it does not solve the exact
full-wave dispersion equation or select a mode automatically.
Only a single bare exterior wire is supported.

The admittance inverse and material prefactor follow independently
from voltage definitions (10)–(11), and agree with the explicit image
formula (23). Both survey records now include the previously omitted
Appendix D. The printed main-text equations remain distinguishable
from their numerical interpretation.

## Xue exact four-layer coefficients

`EarthImpedance.Formula(:Xue2021)` evaluates the exact magnetic
response through the shared permeability-weighted layer calculation.
Its air-reference longitudinal prescription distinguishes it from the
default Tsiamitros2008 selection. This is not an EHEM substitution.

`EarthAdmittance.Formula(:Xue2021; evaluation=:exact)` supplies
the matching four-layer potential kernel, including the magnetic
forcing in the electric interface equations. The existing default
`evaluation=:ehem` remains distinct.

Both exact records include their source identification tables and an
evaluable elimination of the original appendix equations. Independent
comparisons use the final publication's full expressions, unequal
permeabilities, homogeneous and two-layer limits, and full matrix
assembly. The source's four-layer scope is retained; this registration
does not claim an arbitrary-layer potential formula.

## Pires numerical evaluator

The [impedance](external-impedance/2026/double-exponential-quadrature/Pires2026.md#identification-and-source)
and [potential](external-admittance/2026/double-exponential-quadrature/Pires2026.md#identification-and-source)
records map to `Xue2018b` in their respective engine families, with
`quadrature=:double_exponential`. One numerical implementation serves
both families; no duplicate physical kernel or Pires formula identifier
is introduced. The default adaptive evaluator remains available.

Lossless-medium branch points split the quadrature intervals, and the
shared outgoing root convention preserves negative-frequency conjugacy.
The source-map interpretation, precision-scaled endpoints, and numerical
stopping controls are documented separately from the printed equations.

## Core proximity and infinite pipe walls

`PipeImpedance.Formula(:Kane1995)` shares the finite-wall kernel with
DaSilva2006 and adds Kane's separately evaluated core-proximity sum.
The sum also supplies DaSilva2006 Method 2; the later publication does
not create another core-proximity implementation. Equal core radii
and materials are required by the reciprocal matrix route. The
directional unequal-core pair expression is not averaged.

Høidalen (36) is a distinct correction. Its matrix selection uses
the published symmetrical three-core factors: twice the pair term
on the diagonal and once off the diagonal. The low-frequency
subtraction is evaluated within the sum to preserve precision.

DaSilva2006 and Høidalen infinite-wall selections supply core-to-pipe
coefficients, with an infinite outer-radius argument. They do not
supply a finite exterior or transfer term and cannot silently enter
a mixed finite/infinite terminal assembly. The separate low-frequency
record explains the first-harmonic sign interpretation from its
unexpanded parent; the printed equations remain available unchanged.

## Kikuchi and Mariscotti scalar potentials

`EarthAdmittance.Formula(:Kikuchi1957)` and
`EarthAdmittance.Formula(:Mariscotti2019)` share their scalar spectral
evaluator. They remain separate registrations because the source scopes
and self-reference geometry differ. Kikuchi uses the lower surface of
one overhead wire; Mariscotti supplies same-medium self and mutual
coefficients in either half-space. Neither registration supplies a
mixed source/receiver matrix.

Both return charge-normalized potential coefficients. The source
normalizations, outgoing Hankel connection, and Mariscotti image-order
interpretation are explicit in the corresponding records. The default
air reference is distinguished from an independently prescribed
longitudinal wavenumber. The shared air-side limiting kernel does not
make the distinct earth-side expression equivalent to Xue2018b.

## Bonded conductor and semiconductor surfaces

The [Ghosh2019](internal-impedance/2019/double-layer-semiconducting-screen-effective-surfaces/Ghosh2019.md#identification-and-source)
and [Ghosh2022](internal-impedance/2022/n-semiconductor-screen-effective-surfaces/Ghosh2022.md#identification-and-source)
internal records use the interface elimination obtained from Ametani2004's
complete terminal matrix (5)–(8). Both selections are aliases of
`InternalImpedance.Formula(:Ametani2004)`. No duplicate surface kernel is
registered. The source records preserve their different assembly scopes.

Physical concentric metal layers are retained before homogenization.
Attached screens enter the appropriate inner or outer effective surface;
only the remaining insulating gap enters the separate magnetic term.
The original radial shunt layers remain intact. N-screen assembly uses the
existing loop-to-terminal transformation with the effective surface terms.
Nonradial metal reductions and unbonded screens are rejected by this route.

The Ametani transcription now retains the original printed fraction bars
and missing square in (9)–(10). Its executable interpretation follows the
preceding complete matrix and the later explicit Ghosh expressions.
Tests cover independent boundary excitations, direct full-circuit
elimination, dc and homogeneous limits, material contrast, temperature,
and retained radial shunt contributions.

## Linear ACSR reduction

`InternalImpedance.Formula(:Merkushev2015)` retains the source's
anisotropic six-strand correction and linear steel-core magnetization.
Its cylindrical Bessel functions are shared mathematical primitives, not
evidence of equivalence to a homogeneous Schelkunoff conductor. The
physical wire arrangement, two materials, and pitch remain explicit.
Morgan1965's complete current-dependent magnetic-loss model is deferred;
it is not substituted for this linear reduction.

## Rallis complex-image evaluations

The three Rallis2013 records map to `Pollaczek1926` with
`evaluation=:dcim` or `evaluation=:dcim_two_level`. The overhead,
mixed, and underground leaves share one matrix-pencil fitter but retain
their distinct sampling variables and integrated image sums. No duplicate
physical formulation is registered, and the exact spectral default is
unchanged. These are finite fitted approximations, with the controls and
large-separation limitations documented in each source record.

## Remaining registry consolidation

The former Ametani1974 impedance registration repeats Nakagawa1973's
two-layer reduction. Its selector is now an alias, while [Ametani1974](@cite)
remains an auxiliary source. The three-layer parent is not replaced by an
arbitrary-layer claim.

The Magalhaes2018 impedance selector now resolves to Xue2018b. Their
retained series kernel is identical. Xue2018b's ground-surface potential
also delegates to the Magalhaes2018 evaluator; other voltage references
remain distinct. The six recovered source records are
[Wedepohl1973](internal-impedance/1973/hollow-shell-hyperbolic-approximation/Wedepohl1973.md#identification-and-source),
[Zhao2020](internal-impedance/2020/hollow-shell-leading-asymptotic/Zhao2020.md#identification-and-source),
[Alvarado1983](external-impedance/1983/overhead-cubic-complex-image-correction/Alvarado1983.md#identification-and-source),
[Bridges1995](external-impedance/1995/buried-self-leading-logarithm/Bridges1995.md#identification-and-source),
[Theethayi2007](external-admittance/2007/impedance-derived-matrix-extension/Theethayi2007.md#identification-and-source),
and [Ametani2021](external-admittance/2021/classical-transmission-line-reference/Ametani2021.md#identification-and-source).

The two shell approximations share scaled hyperbolic factors but retain
their different curvature and transfer-radius terms. The Theethayi
potential recipe shares its impedance with the impedance family and its
normalization with Vance. Its matrix extension is attributed explicitly
to the IET book, not to the original scalar relation alone. Insulation and
semiconductor Ametani2004 selections share one constitutive evaluator.

The former Pollaczek1926 admittance registration combined a classical
overhead coefficient with an unsupported buried Bessel expression.
Its selector now resolves to the classical Ametani2021 reference:
overhead space potential is retained, and buried external and mixed
corrections are zero. Physical buried insulation remains in the separately
assembled radial potential. This corrects the earlier attribution; it
does not assert zero physical mixed-earth coupling.

`IdealGround` is a nonpublication boundary reference, not an additional
author-attributed formula. Its zero external potential must be combined
with any physical insulation potential. It differs from Ametani2021 for
overhead conductors.

Legacy constructor names resolve to the same retained identifier in both
keyword and fully specified forms. Explicit route overrides are preserved.
