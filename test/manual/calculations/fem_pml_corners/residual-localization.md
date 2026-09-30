# Saved-system residual localization — 2026-09-30

The candidate's large residual is concentrated in **axial equations in the left
earth PML strip and scalar continuity equations in the lower-left earth PML
corner**. Terminal balance equations are satisfied to much smaller absolute
residuals. The largest failing axial row has exactly the same nonzero matrix
coefficients as its counterpart in the structured reference. The candidate
solution fails that unchanged equation. This localizes an algebraic failure;
it does not identify a defective production operation or establish a faster fix.

## Evidence and mapping

This inspection used the eight saved solutions from the [fixed-system
comparison](algebraic-results.md): two meshes, original/three-correction
configurations, and two sources, all aerial at 0.1 Hz and prescribed Γ=0.
No mesh, assembly, FEM execution or solver correction was performed here.
The retained finite-Γ formulation matches the current production source; the
Γ² terms are inactive in this particular fixture.

For every solution, `r = b − A₀x` was evaluated using the saved original,
pre-scaling sparse matrix and exact binary RHS/solution values, with 128-bit
arithmetic for row accumulation. Percentages below are contributions to
`sum(abs2, r)` in the original assembled coordinates. They are **not an energy
norm**: different equation families have different physical scaling.

The [Julia inspection](localize_residual.jl) maps saved `.pre` rows through
GetDP's basis numbering and imported mesh order. GetDP duplicates elements with
multiple physical memberships and sorts the imported tags: `.pre` element
indices must not be mistaken for original `.msh` tags. Reconstruction was
checked against every saved signed edge incidence: 542,192 for the structured
mesh and 375,788 for the candidate. All 359,568/248,632 unknown rows were mapped.
The format/order was checked in the [GetDP 3.5.0 source](https://getdp.info/src/getdp-3.5.0-source.tgz),
`GeoData.cpp:409`, `GeoEntity.cpp:529,630`, `DofData.cpp:420`, and
`ProParser.y:2622,3037`.

| Equation family | `.pre` basis code | Structured rows | Candidate rows | Current definition |
|---|---:|---:|---:|---|
| Axial `a` | 1 | 89,852 | 62,118 | [perpendicular edge basis on nodes](../../../../ext/LineCableModelsGmshExt/getdp/quasi-full.pro#L195) |
| Terminal `U` | 2 | 2 | 2 | [grouped conductor basis/current balance](../../../../ext/LineCableModelsGmshExt/getdp/quasi-full.pro#L206) |
| Transverse `bt` | 5 | 180,269 | 124,801 | [edge basis and tree gauge](../../../../ext/LineCableModelsGmshExt/getdp/quasi-full.pro#L222) |
| Scalar `v` | 6 | 89,443 | 61,709 | [nodal basis](../../../../ext/LineCableModelsGmshExt/getdp/quasi-full.pro#L234) |
| Terminal `V` | 7 | 2 | 2 | [grouped scalar basis](../../../../ext/LineCableModelsGmshExt/getdp/quasi-full.pro#L234) |

The prescribed currents, outer Dirichlet values and constrained tree edges are
eliminated DOFs, not independent residual rows. Both meshes have 1,190 fixed
`a` and 1,190 fixed `v` DOFs; fixed `bt` counts are 90,827/63,093. The two known
currents enter the terminal balances. Small terminal residuals do not validate
the physical appropriateness of the eliminated boundary/gauge constraints.

Supporting regions are unions of the cells on which each test function enters
the active weak forms. Cross-region rows have explicit union categories, never
an arbitrary single-region assignment or double counting. For example,
`earth/PML/left | earth/PML/left+bottom` retains 321 candidate interface rows.
Their source-1 original squared residual is 2.3823e−9; all cross-region
categories together contribute 0.0007068% of that solution's total.

## Where the residual resides

| Mesh / configuration / source | Σ\|rᵢ\|² | Axial `a` | Scalar `v` | Transverse `bt` |
|---|---:|---:|---:|---:|
| Structured / original / 1 | 8.86773e−20 | 71.8906% | 18.4427% | 9.66665% |
| Structured / original / 2 | 9.36721e−20 | 65.1747% | 17.7799% | 17.0453% |
| Structured / three corrections / 1 | 1.03855e−20 | 0.00000794% | 80.8196% | 19.1803% |
| Structured / three corrections / 2 | 2.51855e−20 | 0.00000561% | 60.0644% | 39.9356% |
| Candidate / original / 1 | 0.0629977 | 61.1423% | 38.7512% | 0.106453% |
| Candidate / original / 2 | 0.161123 | 56.8925% | 42.9799% | 0.127631% |
| Candidate / three corrections / 1 | 0.305521 | 32.3384% | 67.1866% | 0.474978% |
| Candidate / three corrections / 2 | 0.000231542 | 13.7343% | 85.4635% | 0.802241% |

Terminal `U` plus `V` contribute less than 6.50e−7% in every structured
solution, and less than 2.68e−24% in every candidate solution. These balances
do not account for the reported large residual.

| Candidate configuration / source | Left earth strip | Lower-left earth corner | Combined |
|---|---:|---:|---:|
| Original / 1 | 60.8517% | 38.8574% | 99.7090% |
| Original / 2 | 56.6508% | 43.1074% | 99.7582% |
| Three corrections / 1 | 32.1223% | 67.6611% | 99.7835% |
| Three corrections / 2 | 11.7107% | 86.2504% | 97.9611% |

The structured residual is much smaller in absolute magnitude and has a
different distribution: its original solutions are concentrated in the upper
right air corner and right earth strip. After corrections, the right earth
strip dominates, with a physical-air contribution for source 2. Complete
family, region-union and family-by-region squared sums, counts and top twelve
rows for **each** solution are retained in
[the localization directory](../../../../.linecablemodels/fem/pml-corner-algebraic/residual-localization/).

## Dominant rows and the equations they fail

All rows below have **RHS zero**. Coordinates are mesh `(x,y)` in metres.

| Candidate row | Family / support | Node and coordinates | Original source-1 share |
|---:|---|---|---:|
| 9357 | `a`, left earth strip | 17313, (−19726.0364, −696.7934) | 45.8567% |
| 9035 | `a`, left earth strip | 16991, (−22094.4907, −74.2221) | 14.9936% |
| 227060 | `v`, lower-left earth corner | 92740, (−19140.0240, −16813.1338) | 9.71622% |
| 225945 | `v`, lower-left earth corner | 91625, (−21263.6357, −17176.5462) | 8.73506% |

Row 9357 supports eight triangles, original mesh tags 26942–26945 and
27038–27041. At Γ=0 it contains only axial-column coefficients. With
`D=sₓsᵧ`, its equation is

\[
(Ax)_i=\sum_j a_j\int_{\omega_i}
 [\nu_P\operatorname{curl}(N_j\hat z)]\cdot\operatorname{curl}(N_i\hat z)
 +(j\omega\sigma-\omega^2\epsilon)D N_jN_i\,dS=0.
\]

The terms are the current [magnetic stiffness](../../../../ext/LineCableModelsGmshExt/getdp/quasi-full.pro#L278),
[conductive mass](../../../../ext/LineCableModelsGmshExt/getdp/quasi-full.pro#L282),
and [displacement mass](../../../../ext/LineCableModelsGmshExt/getdp/quasi-full.pro#L300).
There is no conductor-drive basis on this support. The PML coefficients come
from [pml.pro](../../../../ext/LineCableModelsGmshExt/getdp/pml.pro#L13).
Native `Vol` geometry and the retained 12-point triangle rule `I1` apply.
The saved matrix combines these three `a–a` terms; their separate assembled
matrices were not retained, so the individual weak-term integrals cannot be
recovered uniquely from this sum.

For original source 1, the largest products in row 9357 are:

| Column | Aᵢⱼ xⱼ |
|---:|---:|
| 9357, diagonal | +0.280715034 − j0.675640511 |
| 9358 | −0.054984631 + j0.326681516 |
| 9356 | −0.061053649 + j0.298535020 |
| **Sum of all columns** | **+0.162228619 − j0.050700669** |

The diagonal is `(2.20032665 − j2.07984802) × 10¹²`, multiplying
`aᵢ=(2.20666421 − j0.98479875) × 10⁻¹³`. This is a finite imbalance in the
saved equation, not loss from decimal printing or binary64 residual summation.

**The matched structured row is 16121. Every nonzero coefficient is exactly
equal after matching columns by basis and physical position (maximum difference
zero).** Its original
source-1 residual is `5.80e−25 − j3.18e−24`, versus the candidate's
`−0.162229 + j0.050701`. Row 9035 likewise matches structured row 15799 with
zero coefficient difference. Thus changed local mesh/coefficient/quadrature
values in these two strip equations do not explain their changed residuals.
Changes elsewhere in the global system can still affect their solutions.

Row 225945 supports corner triangles 194656, 194658, 194724, 194745 and 194746.
Writing `κ=σ+jωε` and `T=diag(sᵧ/sₓ,sₓ/sᵧ,D)`, its scalar continuity equation is

\[
(Ax)_i=\underbrace{\int_{\omega_i}\kappa T\nabla v\cdot\nabla N_i\,dS}_{C_v}
+\underbrace{\int_{\omega_i}j\omega\kappa T b_t\cdot\nabla N_i\,dS}_{C_b}
+\underbrace{\int_{\omega_i}j\omega\kappa D a N_i\,dS}_{C_a}=0.
\]

These are exactly [lines 348–353](../../../../ext/LineCableModelsGmshExt/getdp/quasi-full.pro#L348).
Here the column blocks distinguish the three assembled terms:

| Saved solution | Cᵥ | Cᵦ | Cₐ |
|---|---:|---:|---:|
| Original source 1 | 1.217252 + j0.068928 | −1.143075 − j0.068139 | 7.75e−14 + j2.14e−12 |
| Original source 2 | 1.468343 − j0.300564 | −1.366104 + j0.258304 | 1.53e−12 + j2.55e−12 |
| Three corrections source 1 | −0.501052 − j1.944821 | 0.627755 + j1.747275 | 2.72e−12 + j3.24e−12 |
| Three corrections source 2 | 0.008790 + j0.031759 | −0.012666 − j0.035268 | 7.79e−14 + j5.55e−15 |

There is cancellation within these blocks as well: the original source-1
scalar diagonal product is `−284.290049 + j6.332094`, before adding the other
scalar columns. `Cᵥ` and `Cᵦ` are gauge-dependent contributions, not independent
physical errors. Their sum with `Cₐ` is the equation being checked.

The already completed corrections change the failing rows as follows:

| Saved solution | r₉₃₅₇, axial | r₂₂₅₉₄₅, scalar |
|---|---:|---:|
| Original source 1 | −0.162229 + j0.050701 | −0.074177 − j0.000789 |
| Original source 2 | −0.261873 + j0.131681 | −0.102239 + j0.042260 |
| Three corrections source 1 | −0.137354 + j0.188290 | −0.126703 + j0.197546 |
| Three corrections source 2 | −0.000302 + j0.005049 | +0.003876 + j0.003509 |

Full-precision RHS, every individual nonzero product, block sums, mesh support
and all four solution comparisons for the union of dominant rows are in
[`metric-1.35-row-terms.toml`](../../../../.linecablemodels/fem/pml-corner-algebraic/residual-localization/metric-1.35-row-terms.toml).
The [matched strip comparison](../../../../.linecablemodels/fem/pml-corner-algebraic/residual-localization/unchanged-strip-row-comparison.toml)
preserves both exact stencils and all reference/candidate residuals. Top twelve
candidate rows account for 96.77–98.94% of the squared residual across solutions.

## What this establishes, and what remains unidentified

The saved candidate vectors fail homogeneous axial and scalar continuity
equations, predominantly in the two identified earth PML regions. This failure
exists before voltage-path extraction or inversion of P. It cannot be explained
solely as cancellation during that output calculation. Nor is it predominantly
a terminal-current constraint failure. Its largest row lies outside the
remeshed corners, in an unchanged local stencil.

The mapping does **not** establish how much each residual row contributes to
G₁₂ error, or that the structured field/output is converged. It also does not
identify which operation in scaling, factorization, correction or solution
recovery creates the inaccurate candidate vector, or which changed global
couplings make that execution sensitive. The saved original/working matrices
and final vectors do not contain intermediate scaled solution vectors, LU
factors/pivot history or per-weak-term assembly matrices. Consequently they
cannot distinguish a faulty operation from numerical loss in executing the
coupled system, or assign a separate assembly-term defect.

**No production correction is justified by this inspection.** No further solve
or tuning is scheduled here. A reliable smaller system could eventually permit
lower runtime; neither reliability nor that runtime gain has been delivered.
Of the 105 previously recorded file hashes, 104 still match, including the
production FEM sources and original experiment files. The manual runner differs
from that older checkpoint (runner modification time 19:24:45, checkpoint
19:11:15, local time); this inspection did not edit it. New files are this note,
the manual arithmetic script and its inspection artifacts; original results
and logs are preserved.
