# Historical mathematical audit and experiment record

The [locked executive plan](unified-earth-return-implementation-plan.md) is authoritative. The user accepted the original thin-wire FEM comparison on 2026-09-09. Earlier proposed FEM gates and rollout questions below are historical and do not reopen that acceptance. The technical derivations remain implementation references.

# Unified earth-return formulation: mathematical audit and implementation strategy

Status: implementation design and exploratory numerical audit; production defaults
have not been changed. Baseline: commit
`2d694a2444f52942e9d14cd5ff835760295ff532`, which is the exact commit referenced in
[PR 47](https://github.com/Electa-Git/LineCableModels.jl/pull/47/changes/2d694a2444f52942e9d14cd5ff835760295ff532).

The supplied September 2026 manuscript is the proposed mathematical specification.
Its statements about an accompanying 65-check Python script are claims in the
document; that script was not supplied. The checks below are independent work.

The subsequent [manual earth-matrix comparison](earth-matrices-manual-comparison.md)
excludes FEM and all internal/insulation terms. It reproduces the user's 80
reference values and distinguishes their isolated primary-current factor Cf
from the complete manuscript current map L. Both conventions are calculated
explicitly; they must not be mixed when judging mutual impedance and admittance.

## What the prototype actually evaluates

The candidate earth response is transcribed from the supplied manuscript. The
historical earth formulas were evaluated separately as a comparison baseline;
their resulting Z/P/Y entries were never inputs to the candidate current map
or its final Z/Y solves. The Q and MU integrals are the manuscript's own buried
block, even though their Γ=0, equal-permeability limits also occur in earlier
literature. Their use here does not substitute an author's finalized formula
for the manuscript's full-current normalization.

| Candidate operation | Manuscript equation label |
|---|---|
| Direct/self and target circumference average | `eq:Lambda-self`, `eq:Lambda-mutual` |
| Series source kernel | `eq:Q-same`, `eq:Z-pair-same` |
| Deep-earth voltage source kernel | `eq:H-gg` |
| Scalar diagnostic, separate from the selected voltage | `eq:Mphi-same`, `eq:Pphi-pair-same` |
| Physical enclosed-current map | `eq:Ar-Fr`, `eq:D-implementation` |
| Independent enclosed-current check | `eq:D-self`, `eq:D-same` |
| Physical Z, P, Y | `eq:global-PY`, `eq:global-Z`, `eq:right-solves` |
| Copper internal impedance, in the corrected runner | `eq:internal-impedance` |

The first runner reused the existing copper internal-impedance calculation and
read material values through a default-formula workspace. The corrected runner
reads material properties directly from the case and evaluates the manuscript's
internal equation with scaled Bessel functions, including conductor displacement
current. It needs no registered impedance/admittance selection. Historical
formulas now run only with the explicit `--legacy-baseline` comparison option.

This remains an exploratory **equal-radius, equal-depth buried Γ=0 specialization**.
It does not implement or validate the entire supplied framework. In particular,
its FEM discrepancies and individual integration-check failures must not be
presented as a verdict on the complete air/earth/mixed, arbitrary-Γ formulation.

## Mathematical conclusion

The proposal requires a system-level earth-return assembly operation. Replacing
the two `default.jl` pair formulas alone would omit the physical-current
normalization and would not implement the manuscript.

Use the following equations, with `s = jω`, and solve on the right:

```math
\kappa_m^2=\gamma_m^2-\Gamma^2,\qquad
\gamma_m^2=s\mu_m(\sigma_m+s\epsilon_m),
```

```math
K=\mathcal Z-\Gamma^2\mathcal P_\phi/s,\qquad
L=A_r^{-1}-F_rK,
```

```math
PL=H,\qquad (Z-Z_c)L=K+\Gamma^2H/s,\qquad YH=sL.
```

`Pphi`, `H`, `P`, and `L` are different operators. In particular, differentiating
the path-voltage kernel to obtain leakage would implement a different current
map. Both final external impedance and potential coefficients must be referred
to total enclosed currents. Adding the internal conductor matrix happens after
this normalization; multiplying internal impedance by `L⁻¹` would count the
normalization twice.

The self part of `L` follows from the modified-Bessel Wronskian
`x*(I0(x)*K1(x)+I1(x)*K0(x))=1`; this agrees with
[DLMF 10.28.2](https://dlmf.nist.gov/10.28.E2). The regular-field part follows
from differentiating `I0(κr)`. This independently explains why `F_r` tends to
`π*σhat*r²`, while `A_r⁻¹` tends to the identity, as `κ` tends to zero.

The independent audit script solves the four original interface equations and
compares their solution with the explicit R/T coefficients, the axial Q/Mphi
split, and both vertical-field reduction identities. It uses unequal
permeabilities, nonzero air conductivity, and nonzero complex Γ. The buried
benchmark additionally constructs `L` from reflected Hertz fields and from
`A_r⁻¹-F_r*K`, through independent integrals.

The physical reciprocal quantity is `Zc + K/L`. General path-defined `P`, `Y`,
and `Z` can be directional. Preserve both ordered matrix entries. Reciprocity
tests for the equal-radius buried benchmark are appropriate, but blanket
symmetry/passivity tests on every path-voltage mixed matrix are not justified.
The present native frequency loop already preserves ordered entries; it does
not need a new symmetrization step.

The manuscript is a **two homogeneous half-space** model with air/earth/mixed
placements. It is not a solution for arbitrarily many stratified earth layers,
conductors crossing an interface, proximity redistribution, or exact thick
finite-cylinder scattering. These distinctions belong in applicability and
validation, not in numerical integration tolerances.

## Registration migration without changing PR 47's selection model

| Owner | Existing default branch | Retained author registration |
|---|---|---|
| EarthImpedance | overhead | `:Wise1934` |
| EarthImpedance | buried | `:Xue2018` |
| EarthAdmittance | overhead | `:Wise1948` |
| EarthAdmittance | buried, infinite-depth reference | `:Xue2018` |

Move the current equations, Γ prescriptions, permeability hooks, numerical
options, and applicability into the corresponding author files unchanged.
Preserve the overhead contrast approximation and the legacy self-radius
substitution in those retained formulas. Do not quietly improve their math
during the registration move. The source file continues to return its Symbol,
and each owner continues to discover files in its `formulas/` directory.

The existing Xue bibliography names a Delft thesis. The primary source is
Haoyan Xue's 2018 Polytechnique Montréal thesis,
[General Formulation and Accurate Evaluation of Earth-Return Parameters for Overhead / Underground Cables](https://publications.polymtl.ca/3190/1/2018_HaoyanXue.pdf).
Correct that citation when moving the registration; match the retained equations
to the source before stating equation numbers. Do not attribute the entire new
derived framework to Xue or Wise.

Keep `formula(:default; parameters, hooks, options, equivalent_earth)`,
`Formula{ID}`, `FormulaMethod`, and the indexed
`(kind, source_layer, target_layer)` equation declarations. The new `:default`
must declare self/mutual `(1,1)` and `(2,2)`, plus mutual `(1,2)` and `(2,1)`.
Declare exactly two half-spaces and arbitrary prescribed Γ support. An omitted Γ
can retain zero as the explicit transmission-line approximation; no modal root
solver is implied by this work.

Retain the public `air/earth/mixed` NamedTuple selection and numerical options
under the existing `integration=(method=:quad|:trapz|:cim, options=(...))` section.
Voltage reference is a physical parameter, not an integration option. Proposed
reference values are `:deep_earth`, `:interface`, and a finite depth, normalized
to a concrete value/tag before frequency evaluation. Adopt common deep earth as
the package default for all placements, and retain scalar potential as an
explicit diagnostic/convention. This is an intentional overhead-reference
change, to be documented with the default migration.

## System assembly and concrete code boundaries

| Location | Required change |
|---|---|
| `earthimpedance/formulas/` and `earthadmittance/formulas/` | Preserve authors; register the new default's indexed equations and declarations. |
| `earthkernels.jl` | Shared stateless material/branch, Bessel, and spectral primitives; no routing to another author's formula. |
| `input.jl` | Preflight the complete physical system and allocate default-kernel matrices, radius state, and factorization storage. |
| `earthreturn.jl` | Add a formulation-owned system-assembly path for the default; retain the existing entrywise path for author formulas. |
| `lineparameters.jl` | Prepare reusable default earth response once per frequency, after constitutive evaluation and before Z/P assembly. |
| `impedance.jl`, `admittance.jl` | Consume finalized external Z and P; keep cable contributions and established reduction order. |
| `integration.jl` | Share the current integrator interface, with explicit endpoint-weight/contour extensions where required. |
| provenance, formula contract tests, documentation, gauntlet | Record all new formula dependencies and the physical voltage/Γ choices; update expected registrations. |

1. **Preflight geometry and convention.** `EarthPair.layers` and `heights` are
   source first, target second, whereas manuscript superscripts are receiver
   first. `row` is the receiver and `column` is the source. New code must map this
   explicitly. The code uses horizontal `x` and vertical `y`; the document uses
   horizontal `y` and vertical `z`. Keep the code's coordinate vocabulary.
   Existing `_geometry` also reverses the document's D/d labels and evaluates a
   self interaction as a mutual one at the radius. Do not use that helper for
   the new default. For self terms, separation is zero and image distance is
   exactly `2h`; direct self is `K0(κr)` without an extra I0.

2. **Supply radii for every ordered interaction.** Mutual `EarthPair` has
   `radius=nothing`. Derive one radius vector per physical exterior assembly
   from the existing blueprint/self data, and use `r[row]` for target averaging.
   This avoids changing the author formulas' geometry contract. Require
   `abs(height)>radius` and disjoint circumferences. Mixed exponential factors
   must use separate material decay constants, never an image-distance shortcut.

3. **Prepare one concrete frequency state.** Cache σhat, γ², κ², selected
   outgoing branches, I0-related scale factors, and `I1/(κ*I0)` for each radius.
   Cache centre kernels for exchange-related pairs only when the medium state
   and geometry actually permit it. Never infer equality of directional
   endpoint terms from equal centre distances.

4. **Assemble the full auxiliary matrices.** Use indexed declarations and
   formula-owned dispatch to construct `Zcal`, `Pphi`, and `H`. The shared
   assembler forms `K` and `L`, then factorizes for the right solves. The existing
   scalar functor API cannot finalize a default entry from an isolated pair;
   require full-system prepared context for that path, with a clear diagnostic
   for context-free use. Keep author functors' scalar contract unchanged. Do not
   smuggle a two-wire inverse into a mutual-pair callback.

5. **Reuse the completed response.** `earth!(...)` copies finalized default
   entries, while author selections keep their current entrywise execution.
   A default leaf in an indexed selection requires its own complete auxiliary
   system even if only some output entries are selected. A retained author
   entry must never be inserted into `K` or `H` as though it were a new-model
   source kernel. Cross-owner caching is permitted only for identical resolved
   materials, Γ, reference, and hooks. Independent Z and P selections remain
   independent; a hybrid result is not entitled to the unified modal identity.
   A `contribution` override must have a documented stage and normalization;
   preserving the existing scalar output contract means applying it to finalized
   entries, and recording that it overrides the unified response.

6. **Keep normalized P in the native pipeline.** Native `admittance!` assembles
   potential coefficients in m/F and `_solve!` already multiplies the inverse
   by `jω`. Supply `H/L` there; supplying `H/(jω*L)` would introduce an erroneous
   extra `jω`. Use general LU, not Hermitian/Symmetric factorizations. A right
   solve `XA=B` can use `transpose(A) \ transpose(B)` with an ordinary transpose.
   Reuse the L factorization for the external Z and P solves.

The current internal-impedance law uses `sqrt(jω*μc*σc)`, without Γ or conductor
displacement current. At the benchmark's Γ≈0 this is an excellent conducting-core
approximation, but arbitrary Γ support for the complete manuscript also requires
its solid-core `κc/(2πr*σhatc) * I0(κc*r)/I1(κc*r)` relation. Extend the internal
state through the established formula interface, or explicitly retain and record
the quasi-TEM internal approximation. Do not claim that changing earth kernels
alone implements the entire arbitrary-Γ conductor boundary condition.

The paired formulas require the same material model across the full auxiliary
system. Existing pair-dependent equivalent-earth reductions do not automatically
define a single Maxwell half-space problem. For the new default, initially
admit an actual homogeneous earth (or a single globally consistent explicitly
reduced earth). Preserve the older pair-dependent reductions on retained authors;
do not silently reinterpret them as the unified field solution.

The existing coaxial pipeline sums exterior cable currents and attaches radial
insulation networks. The manuscript proves the new operator for bare circular
boundaries. Keep that distinction explicit: the first certified scope is bare
solid conductors. Before changing defaults for coated/multishell assemblies,
derive and test the mapping from physical terminal currents to the outer
circumference, including radial voltage drops. Do not promise a full coated-wire
Maxwell solution merely by inserting an insulation outer radius. This is a
rollout prerequisite, not permission to silently choose an old author for an
unsupported default case.

## Integration and numerical stability

For an ordered pair factor the common weight as

```math
e^{-a_t h_p-a_s h_q}
=e^{-\lambda(h_p+h_q)}
 \exp[-h_p\kappa_t^2/(a_t+\lambda)
      -h_q\kappa_s^2/(a_s+\lambda)].
```

The residual fits the existing `SpectralIntegral(:cosine, ...)` representation.
Using the rationalized root difference avoids cancellation at large λ. The
two distinct medium roots must stay in the mixed residual; the existing radial
Sommerfeld image identity is not a generic mixed-medium identity.

| Method | Implementation and accuracy contract |
|---|---|
| `:quad` | Primary numerical reference. Expose transition scales, branch points and extracted poles to interval/contour construction. Use dimensionless scaling and absolute tolerances for near-zero/cancelling kernels. |
| `:trapz` | Retain separate grid and tail refinement. Share the same declared outgoing sheet and admissible contour; agreement of successive grids alone cannot detect a wrong branch. |
| `:cim` | Retain spectral fitting at fixed frequency, held-out residual checks, analytic image integration, and independent quadrature validation. A fit failure remains a failure; never return quadrature labelled as CIM. |

J0 endpoint terms are the main new integration requirement for overhead/mixed
voltages. They are absent from the first buried-only benchmark. One can fit J0
as part of a cosine residual, but that can waste image rank and makes rotated
tail bounds radius-dependent. Prefer an explicit endpoint weight whose image
transform is

```math
\int_0^\infty e^{-c\lambda}J_0(r\lambda)\cos(y\lambda)d\lambda
=\tfrac12\left[((c-iy)^2+r^2)^{-1/2}
               +((c+iy)^2+r^2)^{-1/2}\right],\quad\Re c>0.
```

Continue both roots from positive real c. This keeps the circumference average
analytic for CIM, while quad and trapz evaluate the same declared J0 weight.
Extracted-pole support also needs an endpoint-weight identity or an explicitly
verified alternative fit. Do not reuse the ordinary cosine pole formula for a
J0-weighted pole. Bound complex-contour exponential growth using both separation
and receiving radius. The current angle heuristic accounts only for separation;
it is insufficient for endpoint terms.

The current trapz and CIM paths can rotate, while quadrature uses the real axis.
The new default must not blindly reuse a fixed positive rotation for arbitrary
Γ: locate branch points/poles and establish the deformation/limiting-absorption
prescription. At `κ=0`, evaluate continuous combinations before taking the limit.
The R/T expressions containing `κs²/κt²` need either nonsingular combined kernels
or a scaled solution of the original four interface equations. Setting Γ to an
air propagation approximation can put the air exactly at this degeneracy.

Use scaled Bessel products/ratios before multiplication, including combinations
such as `I0(κr)*K0(κD)`. Merely calling scaled functions and multiplying back
their separate large exponentials defeats the scaling. A shared column scaling
of the auxiliary source basis cancels from final matrices and can help when
`A_r` becomes extreme. The `1/I0` expression also needs care at zeros of I0;
the direct Hertz current-map expression supplies an independent representation.

The present `bessel_difference` uses a fixed small-argument logarithmic switch.
Do not apply it independently to divergent terms in the new kernels. Use a
tolerance-controlled combined series for direct/image differences, preserving
the self I0 distinction and the O(κ² log κ) terms when needed. A κ=0 spectral
scale cannot be passed to `SpectralIntegral`; choose a positive scale from
material contrast and geometry when the medium scale vanishes.

## Performance and scalar types

Allocate O(N²) kernel/RHS storage once per computation. Matrix elimination is
O(N³) per frequency, unavoidable for the proposed closure. Integrals should
remain the dominant cost for small systems; share Q/Mphi/MU material evaluations
and exact algebraic cancellations where the selected physical assumptions permit.
For Γ=0, `K=Zcal`, so Pphi is unnecessary outside scalar diagnostics.

Retain parametric functors and numerical buffers, with a function barrier around
each bound declaration. Avoid `Any` arrays, Symbol-based decisions, and rebuilding
material vectors in production spectral callbacks. Use caller scalar types for
π, tolerances, branch calculations, and Bessel series. The audit script is an
exploratory Float64 implementation and is not evidence of production type stability.

Measure warmed kernel calls, full frequency sweeps, allocations, and repeated
geometries at N=2 and larger N. Test `@inferred` on concrete state preparation,
kernel callbacks and solve helpers, plus the existing workspace type invariants.
Preserve Float32/Float64 for all three methods; preserve applicable BigFloat and
uncertainty routes for quad/trapz. Current CIM explicitly supports only
Float32/Float64 matrix pencils. Complex BigFloat K currently requires positive
real argument part: lossless outgoing-axis support is additional work, not an
existing capability to claim. J0/I0 support for uncertainty types needs the same
extension audit as the existing Bessel helpers.

CIM currently performs an independent quadrature calculation on every fit and
allocates SVD/pencil data. Reusing image buffers alone does not prove a speedup.
Report fit/setup/validation costs separately before optimizing or claiming
performance preservation. Preserve explicit failures when convergence is absent.

## FEM comparison: establish identical observables first

The FEM description in this section records the original coupled operator.
On 2026-09-11, the backend changed to independent magnetic and scalar
electrodynamic blocks at Γ = 0, sharing one system and factorization. See the
[current field equations and extraction](../src/fem.md#field-equations-and-matrix-extraction).
The normalization and observable checks below remain relevant to comparisons.

The supplied case has copper radius 0.0425 m, horizontal separation 1 m, depth
1 m, earth resistivity 0.1 Ω·m, relative permeability/permittivity 1, and
temperature 20 °C. The corrected requested frequencies are exactly
`[0.1, 1, 10, 100, 1000, 10000, 100000, 1000000]` Hz. Use an `ExactOverrides`
variation; retain the catalogue's original 101-point frequency declaration.

The current FEM is the coupled `A_z/u_r/phi` quasi-TEM system in
`ext/LineCableModelsGmshExt/getdp/quasi_tem.pro`; it is not simply a standalone
scalar Helmholtz electrode solve. Nevertheless it does not implement all the
degrees of freedom and voltage definitions of the manuscript:

- `model.pro` fixes Γ to `j*1e-12 m⁻¹`; the public adapter rejects problem Γ.
- `A` has only the axial component; the full transverse vector-potential
  reconstruction in the document is not independently solved.
- GetDP exports `Phi/(Gamma*I)` as raw `P`, in Ω·m. The FEM adapter directly
  inverts this to Y. To compare against analytical P in m/F, use
  `P_analytic = jω * P_FEM_raw`. Compare final Y directly without that conversion.
- GetDP exports `-U/I` as Z. It does not explicitly export the circumferential
  mean of the document's vertical electric-field path integral.
- The FEM uses a terminal equipotential constraint and spatial conductor fields;
  the manuscript retains an axial line source and mean surface projection.
  Discretization refinement does not remove that modeling difference.

Consequently, matching the existing FEM's Z and Y is a measurable acceptance
target, but is not guaranteed by these equations. First compare the new deep-earth
voltage result and its scalar diagnostic to the existing FEM, with the correct
units, material state, current normalization and reference. Do not choose whichever
voltage or normalization happens to improve agreement. If the observables differ,
derive the FEM extraction transformation, or add a separate independently
validated full-field/projection reference for the manuscript. Preserve the
existing FEM baseline and label the new reference distinctly.

For a meaningful acceptance report, retain every self and directional mutual
entry, real/imaginary parts, R, L, G and C at every requested frequency. Report
absolute errors for components near zero and norm-relative errors per matrix;
do not let a full-band RMS hide an individual failed frequency. Include raw
primitive matrices before reductions. Diagnose common/differential two-wire
responses as well as entries, especially when mutual terms decay rapidly.

Proposed numerical gates (engineering choices, not supplied acceptance limits):
quad tightened until changes are below `1e-8` relative at matrix level; trapz/CIM
agreement within `1e-5` plus scale-aware absolute terms; algebraic solve residuals
bounded by conditioning and machine precision. Require mesh/domain refinement
to quantify FEM error before fixing the analytic/FEM tolerance. A provisional
1% matrix error can be reported as a screening metric, but must not be called an
agreed acceptance criterion or silently enlarged to pass.

## Delivery sequence and release gates

1. Capture old-default numerical fixtures, then move those branches to author
   registrations. Verify their equations/hooks/results remain identical across
   supported methods; update formula inventories and provenance dependencies.
2. Implement and independently verify branch/Bessel primitives and the full
   system closure for bare buried conductors. Certify N=1, unequal-radius N=2,
   N=3/4, and the homogeneous-medium limit; equal wires alone conceal ordering
   and radius mistakes.
3. Complete the eight-frequency buried-wire comparison with the current FEM,
   resolve observable/model differences, and establish convergence/error bounds.
   Failure here is a scientific finding to resolve before switching defaults.
4. Implement overhead and both mixed directions, J0 reference endpoints, finite
   reference depths, and supported nonzero Γ contours. Validate interface
   residuals, direct field/path integrations, zero contrast, PEC/classical limits,
   the filament limit, and reference-shift modal invariance.
5. Validate the exterior-current mapping for existing coated/coaxial assemblies,
   hybrid selections, explicit equivalent-earth use, hooks, UQ, and numeric types.
   Run formula architecture/identity tests, retained formula fixtures, spectral
   integration tests, line-parameter integration tests and gauntlet checks.
6. Switch both `:default` registrations together after those gates. Document
   applicability, Γ approximation, common reference, and source-versus-terminal
   matrix roles. Keep every retained author's explicit selection available.

The exploratory runner was retired after the earth formulation entered production.
The observations below describe its historical runs and retained development
results under `.linecablemodels/qa/unified-earth-audit/`. They are not current
Gauntlet acceptance criteria. The maintained regression entrypoint is
`test/unit/engine/unified_earth_return.jl`, which reads the accepted frozen fixture
and checks quadrature, trapz and CIM:

```sh
julia --startup-file=no --project=test test/runtests.jl unified_earth_return spectral_uncertainty
```

## Observations from the initial run

Julia 1.12.7, current checkout, actual `two_bare_wires` case with the eight
frequency overrides. Case source SHA-256:
`ee15401269881036f6851eaaa505d871ca63ecaf9a84d62e6e1cbbd526f1622c`.
`--compiled-modules=existing` permits loading with the existing read-only Julia
depot without creating precompile files there. No package dependencies or
production files were changed.

The FEM executed eight GetDP invocations and completed all 16 excitation columns,
using cached geometry meshes. Its maximum potential-map condition number was
8.525 and maximum inversion residual was `6.67e-16`. These checks establish
successful linear solves, not mesh or domain convergence. This run did not
perform a mesh-refinement study.

Relative Frobenius matrix errors for the new quadrature prototype versus FEM:

| Frequency (Hz) | Z error (%) | Y error (%) |
|---:|---:|---:|
| 0.1 | 0.5454 | 1.0829 |
| 1 | 0.3990 | 1.4190 |
| 10 | 0.3789 | 2.0783 |
| 100 | 0.4258 | 3.4400 |
| 1,000 | 0.7562 | 6.3215 |
| 10,000 | 0.7673 | 10.6155 |
| 100,000 | 0.6803 | 4.1092 |
| 1,000,000 | 0.6180 | 0.5107 |

The matrix norm conceals larger errors in small mutual terms. At 100 kHz the
mutual Y error is 37.62%; at 1 MHz it is 32.87%, with mutual Z error 7.25%.
The corresponding component and entry errors are retained in
`docs/notes/unified-earth-return-initial-validation.tsv`. The largest separate R/L/G/C matrix error
is about 27.90%. There is no basis to declare the requested FEM agreement met.
The scalar-potential diagnostic does not resolve the difference: its Y matrix
error reaches about 37.54% at 10 kHz. Identifying the dominant model/reference
or discretization difference remains required work.

At 1 MHz the old defaults differ from FEM by about 12.40% in Z and 12.93% in Y.
The new normalization substantially improves those particular matrix errors;
this does not establish agreement throughout the requested band.

Verified numerical identities and integration observations:

- Forty interface/algebra checks passed; worst relative residual `2.07e-15`.
- The two independent L constructions agree within `6.44e-15` with quadrature.
  The physical `Y*P=jω*I` residual is below `4.50e-16`.
- For the required Γ=0 Z/Y kernels, trapz passes all eight frequencies and agrees
  with quadrature within `4.22e-8` in external Z or Y matrix norm.
- With `rtol=1e-6, atol=0`, CIM passes the first six frequencies (maximum Z/Y
  difference from quad `1.35e-7`). It fails its independent integral check on
  mutual Q at 100 kHz and self Q at 1 MHz. The reported absolute Q integral
  discrepancies are `3.54e-8` and `9.01e-12`, respectively. These are explicit
  failures, not substituted quadrature results. A physically scaled absolute
  error budget for small corrections is a concrete integration work item.
  A targeted `atol=1e-7` experiment passes at 1 MHz with Z/Y matrix errors below
  `5.76e-10`, but at 100 kHz the mutual MU check still fails (`1.40e-7` integral
  discrepancy). This is evidence for explicit error budgeting, not a completed
  CIM solution over the whole band.
- An earlier exploratory pass also evaluated Mphi through every method even
  though Γ=0 Z/Y do not need it. That made trapz fail at the lowest four
  frequencies and CIM fail at all frequencies. Separating this optional scalar
  diagnostic from the required kernels removed those extra failures. Nonzero Γ
  and scalar-potential operation still require robust Mphi evaluation; the
  buried Γ=0 pass must not be used as evidence that this is already supported.

Timing/allocation values in the exploratory JSON include compilation and scalar
diagnostic work. They are not a warmed performance comparison or a type-stability
certificate. The strategy's performance and general-placement release gates
remain open.

### Correction after tracing the candidate's dependencies

The corrected runner was executed without `--legacy-baseline`, with all candidate
electromagnetic equations taken directly from the manuscript, including the
internal copper relation. Compared with the first candidate, quadrature Z changes
by at most `1.26e-16` relatively and Y is unchanged. Thus the original candidate's
external results did not come from the retained author formulas. This also leaves
the reported FEM discrepancies unchanged at the displayed precision.

That historical run used the explicit dimensionless integral absolute tolerance
`cim-atol=1e-6`.

For this restricted benchmark, CIM now completes all eight frequencies, with
maximum external-Z/Y matrix differences from quad of `1.78e-7`, and maximum
relative mutual-entry difference `4.27e-6`. Trapz also completes all eight.
The previous CIM statement referred to failed relative-only checks on individual
correction integrals; it did not establish a significant final-matrix error or
incompatibility with the framework. The successful run explicitly changes the
integration tolerance, retains independent quadrature validation, and returns
CIM image integrals. It neither supplies a general error-budget derivation nor
validates arbitrary Γ, overhead, or mixed placements.

## Earth-only plots with bare PEC conductors

**Correction after tracing the user's original `thin_wire.pro`:** the custom
quasi-TEM results below do not reproduce the user's historical Electric FEM
experiment and are withdrawn as its reference. The original scalar
Helmholtz/diffusion branch has now been run directly and reproduces the archived
1 MHz Y to numerical precision. See [the provenance and replay report](thin-wire-fem-provenance.md).
The earlier numbers and plots remain diagnostic results of the different
operator described below, not evidence against the supplied framework.

The historical earth-only comparison set analytical `Zc=0` and selected the
library PEC material through an exploratory case parameter. That parameter was
removed with the runner; the two-wire catalogue case retains copper and no
insulation.

The library's `:pec` material is a finite-conductivity approximation
(`rho=eps(Float64) Ω·m`). Running it through the existing conductor-volume FEM
produced an invalid negative self reactance, approximately
`-3.188e12/f Ω/m`. A solver completion marker did not catch this failure. Those
values are preserved only as diagnostic output in `pec/fem.json` and are excluded
from the delivered plots.

The plotted reference instead uses `test/gauntlet/getdp/pec_boundary.pro`, an
explicit PEC-boundary specialization of the existing quasi-TEM equations. The
air/earth meshes and far-boundary conditions are retained. Conductor interiors
are excluded from the integration domain. On each conductor, axial vector
potential `a` and normalized scalar potential `psi=phi/Gamma` have constant
traces; their associated boundary reactions prescribe `I` and `I`, respectively.
The Γ→0 equations are evaluated after normalization:

```math
\int_\Omega \mu^{-1}\nabla a\cdot\nabla v
 + j\omega\hat\sigma a v\,d\Omega=\sum_p I_p v_p,
```

```math
\int_\Omega \hat\sigma\nabla\psi\cdot\nabla w
 + j\omega\hat\sigma a w\,d\Omega=\sum_p I_p w_p.
```

The output is `Z_pq=jω*a_p/I_q`, `Praw_pq=psi_p/I_q`, and `Y=inv(Praw)`.
No conductor conductivity or manuscript Green kernel enters this FEM operator.
This supplies exact PEC boundary conditions within the existing quasi-TEM field
model; it is not a claim that this model solves all the manuscript's full-field
degrees of freedom or uses an identical voltage functional. The reference is
audit-local and does not change the production FEM backend.

The original mesh-preparation runner is retired. The PEC postprocessor and
plotting script can still inspect the saved historical audit files; they do not
supply the production earth regression fixtures.

The Python reference runner requires NumPy. It preserves its generated GetDP
sources, commands, logs, operator hash, and reference matrices in the audit
directory. All eight solves completed; Z reciprocity was checked without
symmetrizing the output, and the worst P/Y inversion residual is `9.52e-16`.
The earth-only Z matrix error is 0.40–0.55%; the maximum Y matrix error remains
10.62%. No mesh-refinement accuracy claim is made.

The figure contains real and imaginary self/mutual Z and Y, with 281 analytical
quadrature points and eight FEM markers. Both signs are retained. PNG, PDF, SVG,
and CSV are written under
`.linecablemodels/qa/unified-earth-audit/pec/plots/earth-only-self-mutual-ZY.*`.
