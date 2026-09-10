# Unified earth return — locked executive implementation plan

**Status: LOCKED AND IMPLEMENTED, 2026-09-09.** The user accepted the completed comparisons as
sufficient to proceed. This document is the implementation contract and
supersedes the exploratory gates in the [historical mathematical audit](unified-earth-return-mathematical-audit.md).
The user authorized execution of the complete revised plan on 2026-09-09.
Implementation is complete in the working tree; validation evidence is recorded below.

The authorized 2026-09-10 performance follow-up adds reusable certified CIM fits
and local trapezoid sampling controls; see
[the implementation and measurement record](spectral-backend-optimization.md).

The user's subsequent sampling and package-reuse requirements are incorporated
below: a shared spectral contract with qualified numerical backends is an early
implementation milestone. The assurance requirements remained fixed; package
choices were resolved through the qualification recorded below.

## Accepted basis

The supplied manuscript, *Unified circumferentially averaged framework for
overhead, buried, and mixed conductor systems*, is the mathematical source.
Its SHA-256 is
`31779aee0372a5be65842869ddd7e527326e1ed41971bcc2f316c11f6bc06db0`.
The implementation follows its derived kernels and total-current elimination.
Earlier authors' finalized formulas are comparison models, not inputs to the
new default.

Accepted evidence comprises the [two-wire earth matrices](earth-matrices-manual-comparison.md),
the [three-wire full matrices](three-wire-earth-matrices.md), and the
[original thin-wire FEM replay and proposal comparison](thin-wire-fem-provenance.md).
All use frequencies `0.1, 1, 10, 100, 1000, 10000, 100000, 1000000` Hz.
The geometry has 0.0425 m radii, 1 m depth, 1 m adjacent separation and earth
resistivity 0.1 Ω·m. The analytical benchmark contains earth matrices only.

The accepted numerical snapshot is
[`test/fixtures/reference/unified_earth_return.toml`](../../test/fixtures/reference/unified_earth_return.toml):
56 cases covering two/three wires, full-current and Cf proposal variants, Xue,
and the fresh original GetDP results. It preserves complete matrices, the
proposal's K/H/L intermediates, units, source hashes and FEM provenance.

The full-current proposal's 1 MHz differences from the original GetDP are
0.22% for self Z, 0.86% for mutual Z, 0.22% for self Y and 1.16% for mutual Y.
These comparisons are **accepted engineering evidence**. No additional FEM
sweep, COMSOL rerun, mesh-convergence campaign or new FEM toy is required to
lock this plan or begin implementation. The original Magnetic branch's
active/inactive conductivity treatment and the Electric branch's scalar
Helmholtz equation remain recorded facts; acceptance does not change their
physical interpretation. The earlier custom quasi-TEM audit is not this baseline.

## Fixed mathematical and API decisions

1. **The new default is the complete manuscript closure.** For `s = jω`, assemble
   the ordered source kernels first, then solve the complete system:

   ```math
   \kappa_m^2=s\mu_m(\sigma_m+s\varepsilon_m)-\Gamma^2,
   \quad K=\mathcal Z-\Gamma^2\mathcal P_\phi/s,
   \quad L=A_r^{-1}-F_rK,
   ```

   ```math
   P_eL=H,\qquad Z_eL=K+\Gamma^2H/s,\qquad Y_eH=sL.
   ```

   At Γ=0 this gives `Ze = K/L`, `Pe = H/L`, `Ye = s*L/H`.
   Use factorizations and solves, not explicit inverses. Cf-only field averaging
   remains labelled comparison evidence; it is not substituted for L.

2. **Cover the manuscript's two half-spaces.** Implement air/air, earth/earth,
   and both ordered mixed directions; arbitrary conductor count, unequal radii,
   heights and medium properties; and caller-prescribed Γ. Default Γ is zero.
   Select the outgoing/decaying sheet before lossless limits. Use common deep
   earth as the default voltage reference, with interface and finite-depth
   references as physical parameters. Scalar potential remains a distinct
   diagnostic. No propagation-constant root solver is part of this change.

3. **Preserve PR 47's architecture.** Keep `Formula{ID}`, `FormulaMethod`, indexed
   `(kind, source_layer, target_layer)` declarations, owner-local registration
   files returning their Symbol, `formula(...; parameters, hooks, options,
   equivalent_earth)`, and the existing `air/earth/mixed` selection structure.
   Resolve physical choices and dispatch before the frequency/kernel loops.

4. **Preserve the old equations under their authors.** Move equations, Γ hooks,
   permeability treatment, self sampling and numerical options unchanged:

   | Owner and placement | Registration |
   |---|---|
   | EarthImpedance, air | `:Wise1934` |
   | EarthImpedance, earth | `:Xue2018` |
   | EarthAdmittance, air | `:Wise1948` |
   | EarthAdmittance, earth, deep reference | `:Xue2018` |

   Correct bibliography metadata without changing retained numerical behavior.
   Switch the two `:default` registrations together at the end of implementation.

5. **Earth assembly is a system operation.** Rows are receivers and columns are
   sources. Obtain target radii from the complete exterior geometry, not the
   mutual `EarthPair.radius`, which is absent. Keep ordered entries; do not
   symmetrize path-defined Pe/Ye. Preserve author formulas' entrywise contract.
   Default callbacks require a prepared full-system context.

6. **Keep existing cable composition.** Produce external Ze and Pe in the
   existing exterior cable-port basis, then retain the established internal,
   insulation and terminal-reduction stages. Native admittance assembly receives
   Pe in m/F and continues to compute `s*inv(Ptotal)`. Validate outer-radius and
   total-current mapping for coated/coaxial assemblies as a code-integration
   check. This work does not replace internal or insulation laws or claim an
   exact coated-cylinder Maxwell solution.

7. **Make applicability explicit.** Circumferences must lie wholly in one
   half-space and must not overlap. The default requires an actual homogeneous
   earth or an explicitly globally consistent equivalent earth. Keep existing
   pair-dependent reductions available through retained author selections.
   Unsupported geometry/material states fail during preflight; no silent
   fallback to a legacy formula. Arbitrary stratification and proximity surface
   harmonics are outside the supplied framework.

## Execution sequence

| Step | Work | Completion evidence |
|---|---|---|
| 1. Preserve authors | Capture current defaults in the existing retained-formula fixtures, split Wise/Xue registrations, and update inventories/provenance. | Old-default equations, hooks and results are preserved under author selections across supported methods. |
| 2. Reuse and qualify adaptive sampling | Add the feature/contour/tail and error contract; reuse QuadGK's public APIs; qualify DoubleExponentialFormulas for trapz and ExpFit for CIM fitting before building additional numerical machinery. Implement only missing coordination/assurance logic. | Candidate choices are supported by entrywise accuracy, type and cost measurements; coverage, local refinement, weighted fit error, true/fitted tails and resource-limit diagnostics pass the sampling cases. |
| 3. Build shared state and closure | Add typed branch/Bessel/spectral primitives and feature providers in `earthkernels.jl`; allocate geometry, K/H/L and solve buffers; connect integral-error budgets to `earthreturn.jl`. | One-, two-, three- and unequal-radius systems exercise complete solves and accepted matrices; small final entries drive adaptive integral accuracy. |
| 4. Complete placements and weights | Implement air/earth/mixed kernels, reference endpoints, prescribed Γ and J0/radial transform capabilities through existing declarations and the shared sampler. | Interface/path identities and limits pass; all three methods resolve each matrix entry, including the small three-wire outer coupling. |
| 5. Integrate and check performance | Prepare the default response in `lineparameters.jl` and consume it in `impedance.jl`/`admittance.jl`; preserve hooks, selective author/default routing, capture and reduction behavior. | Existing architecture and cable checks pass; scalar-type inference, allocations and warmed timings are measured. |
| 6. Cut over together | Make the complete formulation the default in both owners, update documentation and expected default fixtures, and run the affected test suites and gauntlet. | All implementation checks below pass and explicit author selections remain available. |

A default leaf in a mixed formula selection assembles its complete auxiliary
system before selecting final entries. Never place a finalized author entry
into K or H. Apply contribution overrides at the established final-entry stage.
Share Z/P preparation only when resolved materials, Γ, reference and hooks are
identical; otherwise prepare separate consistent systems.

## Shared adaptive λ sampling

The [adaptive spectral sampling design](adaptive-spectral-sampling.md) is part
of this implementation contract. Sampling policy belongs to reusable Engine
infrastructure, not to independent formula-specific grids.

- Formula kernels declare scales, branch/pole features, admissible contours and
  tail information; the controller combines these with the analytic weight.
  Finite black-box sampling alone cannot certify absence of an unseen feature.
- Reuse backend-owned partitions and workspaces through public APIs, with local
  error budgets and independent probes. Refine unresolved regions and extend
  the tail separately. Do not prescribe a new general-purpose interval tree or
  quadrature engine. A global sample count or two agreeing signed totals is
  insufficient, including when that check comes from a package.
- Qualify adaptive transformed trapezoids for trapz; retain a local composite
  rule only where a demonstrated gap requires it. CIM uses adaptive **uniform** pencil
  windows, an adaptively ranked global image fit, integrated weighted residual
  checks and control of its continuation beyond the fitting cutoff. Arbitrary
  nonuniform nodes are not fed into the uniform-grid matrix pencil.
- Return typed integral error/coverage reports internally. Propagate prefactors,
  analytic-term cancellation and matrix-solve sensitivity to final Ze/Pe/Ye,
  then tighten only the influential integrals. A matrix norm must not hide a
  small unresolved mutual entry.
- Shared metadata and workspaces serve retained authors, the new default,
  trapz/CIM and independently refined quad validation. Preserve the existing
  method/options API and truthful method labels. Distinguish estimated error
  from a bound supported by analytic feature/tail information.
- Follow the [package assessment](adaptive-spectral-sampling.md#package-reuse-decision-2026-09-09):
  QuadGK is already available; new DE/fitting packages are candidates, not yet
  benchmarked dependencies. Rational approximation is not silently substituted
  for an exponential image representation. Shared architecture means shared
  metadata and error contracts, with method-appropriate numerical backends.

## Numerical and performance contract

- Retain the public `integration=(method=:quad|:trapz|:cim, options=(...))` API.
  Quad is the independent integration reference; trapz retains grid/tail
  refinement; CIM retains fitting, held-out checks and its analytic transform.
  Never return a quadrature result labelled as CIM.
- Use rationalized root differences, combined small-argument series, scaled
  Bessel products/ratios and consistent outgoing branches. Handle κ=0 through
  continuous combinations. Mixed decay uses the two medium roots separately.
  J0-weighted endpoints require their own transform and contour-growth bounds;
  ordinary cosine image/pole formulas cannot be reused indiscriminately.
- Fix the known CIM small-coupling failure during implementation. The accepted
  three-wire quad/trapz results stand; the old CIM settings' 73% error in Ye13
  at 1 MHz is not an accepted production accuracy level.
- Check individual complex entries and small real/imaginary components, not
  only matrix norms. For the frozen Float64 analytical fixtures, target
  `rtol=1e-7` for the quadrature implementation and `rtol=1e-5` for trapz/CIM.
  Absolute floors are `1e-14 Ω/m` for Ze, `1e-14 m/F` for Pe and `1e-10 S/m`
  for Ye. These are numerical reproduction tolerances, not new FEM-agreement
  requirements. Precision-appropriate tolerances apply to other scalar types.
- Check the right-solve identities and `Ye*Pe = s*I`, with residual bounds
  accounting for precision and conditioning. Preserve the three-wire directional
  Pe/Ye entries and independently checked current-map identities.
- Allocate O(N²) kernel and RHS storage once per computation; reuse factorizations
  and bound constitutive state. The full closure requires O(N³) matrix work per
  frequency. Avoid `Any`, runtime Symbol routing, material reevaluation and
  avoidable allocations in spectral callbacks. At Γ=0 skip unused Pphi work.
- Preserve existing Float32/Float64 support and applicable BigFloat/UQ paths.
  CIM keeps its established Float32/Float64 restriction. Validate concrete
  helpers with inference checks, and measure warmed kernels, complete sweeps,
  allocations and CIM fit/validation costs at N=2 and larger N. Preserve the
  performance of retained author paths; document the new closure's measured cost.

## Definition of done

Both earth owners select the full manuscript formulation by default for the
declared air, earth and mixed cases. All retained author registrations preserve
their behavior. The accepted two/three-wire analytical matrices are reproduced;
all three integrators meet the entrywise numerical contract; interface, reference,
Γ, geometry and cable-composition checks pass. Formula architecture, provenance,
retained-equation, spectral-integration and affected line-parameter/gauntlet
checks pass, including shared adaptive sampling and its adversarial coverage,
tail, aliasing and cancellation cases. Documentation identifies units, current normalization, voltage
reference, applicability and migration from the former defaults.

Routine implementation findings are resolved within this plan. A demonstrated
contradiction in the manuscript or an unavoidable public-contract change must
be reported with the failing equation/case; changing the physics or loosening
numerical tolerances silently is not a resolution. Further FEM research is
follow-up work, not a prerequisite added to the user's accepted benchmark.

## Execution record

- 2026-09-09: Revised package-reuse strategy locked and full execution authorized.
- Wise1934/Wise1948/Xue2018 registrations retain the former default equations;
  all 204 retained-equation checks passed across supported integration methods.
- DoubleExponentialFormulas 0.1.0 was qualified with an adapter and added. ExpFit
  was rejected after scalar-type and zero-kernel failures; it is not a dependency.
  See [qualification evidence](spectral-package-qualification.md).
- Shared spectral features, admissible contours, locally verified DE subdomains,
  CIM weighted residual/tail checks, and numerical estimate propagation through
  complete K/H/L right solves are implemented. Method labels remain truthful.
- The complete manuscript implementation is active in both default registrations
  in the working tree. The accepted analytical matrices and geometry/reference
  identities passed 2,985 checks across all methods before the added intermediate
  K/H/L assertions. A subsequent combined contract run passed 3,140 checks,
  including mixed selections, prescribed Γ and contribution-wrapper choreography.
- The broad engine and line-parameter run passed 6,712 checks, including the
  original 4,096-byte retained-author allocation ceiling. An expanded run passed
  85 public uncertainty checks (including finite differences and correlation-aware
  cache equality) and 182 scalar/air-root checks. The large-argument,
  contour-feature and retained-equation run passed 283 checks.
- The production 1 MHz two/three-wire gauntlet reproduces every real and imaginary
  Ze/Pe/Ye entry under the original tolerances and records warmed backend cost.
  The direct quadrature closure allocates 784 bytes per reused workspace call;
  the public three-cable/two-frequency retained path measured 3,872 bytes versus
  3,616 bytes on HEAD. CIM is substantially more expensive; its measured costs
  are recorded in the package qualification note.
- Quality checks passed 1,185 assertions. Aqua's separate persistent-task /
  cold-precompile check also fails on the untouched HEAD in this environment;
  it first encountered the GeometryBasics/EarCut loader and still failed after
  the native-library workaround. Numerical tests and normal triangulation work.
  This pre-existing environment check is recorded separately from implementation
  failures; it has not been removed or marked as passing.
- Final validation passed 4,395 frozen-matrix and ownership/import checks, 351
  spectral/uncertainty/retained-formula checks, and 432 complete mixed-layer
  reference comparisons. The latter covers prescribed Γ, unequal permeability,
  deep/interface/finite/scalar references and every Ze/Pe/Ye component. The
  production report's source hashes match the final implementation.
- 2026-09-09: Implementation complete in the working tree. The manuscript's full
  current closure is the default for both owners, previous defaults remain
  author registrations, and all numerical acceptance gates passed. The separate
  pre-existing Aqua environment failure remains documented above.
