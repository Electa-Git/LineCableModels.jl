# Shunt model structural refactor

Status: blueprint-ownership revision implemented and focused verification complete.
Known unrelated Aqua ambiguities are recorded below.
The historical quadrature failure was not reproduced in the bounded replay
described below; its exact trigger remains unverified. This is an unpublished
API refactor, not a compatibility repair or a scientific certification.

## Locked design

- `shunt_model` selects the local field/geometry approximation independently of
  the insulation and semiconductor constitutive laws.
- `:default` means `:coaxial`; `:boundary` explicitly selects the existing
  lossless wire/tape boundary approximation. Both LineParameters and
  CableConstants consume the same blueprint coefficients.
- Formulation-aware blueprint construction owns frequency-independent boundary
  coefficients. No separate public preparation object or compute option exists.
- Completed blueprints share equivalent solved domains within one construction
  call, preserving uncertainty dependencies. Fresh compute calls rebuild their
  blueprints; no global mutable cache is introduced.
- Expensive independent boundary and derivative refinement audits are opt-in.
  Production retains convergence, rank, finite-value, and terminal passivity
  checks. A failed numerical solve is never silently accepted.
- Boundary numerical controls belong to its formula. Strict behavior is the
  default; explicitly requested coaxial fallback records its reason and effective
  model. Unexpected exceptions and invalid physical inputs are not swallowed.
- LEP/Monte Carlo comparisons select one consistent model; numerical failure is
  not a reason to reject and resample physical inputs or mix models silently.
- Material-law keywords retain their meanings. No `internal_admittance` alias,
  speculative compatibility wrappers, or equivalent-permittivity projection is
  introduced.

## Execution

1. Add the typed shunt formula family and connect both formulation constructors,
   blueprint construction, result descriptions, serialization, and material-law selection.
2. Separate the boundary audit, normalize quadrature controls, reuse numerical
   scratch, and return informative numerical failures. Reproduce the recorded
   132 kV quadrature failure or extract its smallest failing integral.
3. Construct selected coefficients in the blueprint and retain declared fallback.
   Preserve uncertainty dependencies and consistent differentiation.
4. Update numerical/API tests, documentation, and manual runners. Verify the
   default performs no boundary solve and preserves the ordinary annular model;
   verify audit invariance, explicit correction, fallback and reuse.
5. Measure warm 18 kV/132 kV calls and boundary resolution/tolerance behavior.
   Run focused ordinary, quality, serialization, UQ and Gauntlet checks, reporting
   actual outcomes and limitations.

## Numerical evidence to preserve

The previous 18 kV preparation took about 32--40 s and allocated 1.652 GiB.
About 72% of sampled time rebuilt an independent validation grid. The main
boundary matrix had 6800 rows and 3268 unknowns (169.5 MiB).

Whole-matrix differences are insufficient for choosing resolution: halving the
settings changed the matrix norm by 0.068% but the weak core-to-jacket coupling
by 6.2%; quarter settings changed those quantities by 0.25% and 22.7%.
These are differences from the existing numerical result, not certified errors.

The recorded 132 kV MC failure reports logarithmic quadrature nonconvergence,
but does not identify the realization, evaluation count, or assembly/audit
stage. Its retained previous completed attempt must not be overwritten.

## Completion evidence

### Delivered structure

- Both formulation owners expose `shunt_model`; material-law selections retain
  their existing meaning. Formula descriptions and passive scientific records
  retain the new selection. Missing historical selections are not invented.
- Default calculations use coaxial annuli and make zero boundary solves.
- `CableBlueprint` owns the selected shunt blocks and their numerical outcomes.
  `LineParameters` and `CableConstants` consume blueprint-derived coefficients.
  Boundary calculation takes place during flattening, not workspace construction
  or the frequency loop. Identical local selections share blueprints; independent
  uncertainty sources never share solved matrices.
- Blueprint construction resolves distinct domains and shares identical solved
  blocks. No workspace repeats this construction. Dense kernels reuse numerical
  scratch and release the factorization after construction.
- Boundary settings retain the original resolution. Production uses
  `(rtol=1e-8, atol=1e-10, maxevals=100_000)` for dimensionless logarithmic
  moments. Independent-grid and derivative step-refinement checks are opt-in
  with `audit=true`; absent diagnostic residuals are `nothing`, not zero.
- Recognized failures use `BoundarySolveError`, including design/terminal and
  assembly/audit context. Failed integrals additionally report the actual
  evaluation count, estimate, tolerance and face geometry. Fallback is explicit,
  warned and recorded. Invalid input and unexpected exceptions propagate.
- MC/LEP reject automatic fallback; MC also checks that model coverage stays
  fixed between realizations. Numerical failures are not resampled. Scientific
  UQ records retain model identity; shared-source LEP serialization also covers
  CableConstants without making correlated outputs independent.
- The manual runner exposes the new default. `dev/benchmark_shunt_models.jl`
  exercises ordinary compute calls using timing macros without activating an
  environment or writing campaign results. The local-shunt audit runner separates
  blueprint construction and assembly only for internal numerical inspection.

### Blueprint-ownership verification

- `test/runtests.jl shunt_model internal_shunt`: 242 assertions, 8 items,
  136 s; all passed. Includes shared local domains, global terminal remapping,
  both solvers, Float32 construction, fresh geometry-dependent coefficients,
  correlated Measurements inputs, independent uncertainty sources, strict
  failures, declared fallback, and serialization. Explicit `:lossless` material
  selections produce the same boundary result as their `:default` aliases.
- Cable constants and integration formulation grids: 177 assertions, 4 items,
  165 s; all passed. The lowering counter now observes the formulation-aware
  constructor, and annular runtime-law variants still lower each design once.
- The broader 21-item contract run passed the remaining 1,272 assertions,
  including radial dielectric physics, retained formula records, calculation
  options and UQ/grid constructors. Its four stale blueprint-field/lowering
  expectations were corrected and passed in the 177-assertion rerun above.
- `tag:quality`: 1,225 assertions, 6 items, 102 s; all passed. Display ownership
  now rejects the `Base.Multimedia` fallback while accepting the package's shared
  display implementations. Shunt formulas show their ID, and `BoundarySolveError`
  owns its rich display.
- Gauntlet `records_tests dielectric_consistency_tests`: 235 assertions,
  3 items, 76 s; all passed. No campaign or retained attempt was written.
- `docs/doctest.jl`: the package doctest passed. The complete documentation
  generator was not run because it rewrites unrelated generated pages already
  modified in this worktree.
- `tag:aqua`: 10 checks passed, including persistent-task validation in an
  offline temporary environment. The ambiguity check found two intersections
  in pre-existing literature-ID edits: `chrysochos2014` computation options
  versus generic integration options, and `ametani1980` impedance versus PSCAD
  dispatch. Neither involves blueprint/shunt methods; both were left untouched.

The separate preparation API/file, its snapshot/reuse machinery, the shunt
`selection` wrapper and the UQ metadata packaging wrapper were removed.
Julia docstring conventions informed the revised physical coefficient units,
blueprint ownership and constructor contract. No dependency, Project.toml,
Manifest.toml or consumer-environment change was made.

### Current ordinary-compute timings

Measured with `dev/benchmark_shunt_models.jl`, Julia 1.12.7, one Julia thread,
16 BLAS threads, 101 frequencies from 0.1 Hz to 10 MHz, all terminals retained
and no reductions. Three warm calls per model; no test processes from this
verification were running concurrently. These are shared-machine wall times,
not isolated-CI guarantees. Allocation is cumulative Julia allocation, not
peak memory. Values below are from the inner `@time compute(...)` calls;
the script's outer `@timed` rows additionally include timing-output overhead.

| Case | Default coaxial compute | Boundary compute, including blueprint construction |
|---|---|---|
| 18 kV trefoil | 0.2905–0.2909 s; 12.085–12.086 MiB; 137.61k allocations | 11.015–11.269 s; 295.561–295.562 MiB; 763.77k allocations |
| 132 kV horizontal | 0.5141–0.5164 s; 4.918 MiB; 92.32k allocations | 2.594–2.805 s; 107.127–107.130 MiB; 538.44k allocations |

Every ordinary boundary call above constructs fresh blueprints and includes
the selected local solves. Coefficients are reused across its frequencies and
equivalent local domains, not across separate compute calls. Default annuli
are a different approximation, not an accuracy-equivalent faster boundary solve.

Cold compilation is still substantial: the first 18 kV default compute took
49.338 s (99.39% compilation), and the first subsequent boundary compute took
22.317 s (50.71% compilation). This ownership refactor does not resolve Julia's
cold compilation latency. Restart Julia once after the changed struct definitions;
no consumer-environment update is required.

### Historical timings before the blueprint-ownership revision

Measured in Julia 1.12.7, one Julia thread, 16 BLAS threads, 101 frequencies
from 0.1 Hz to 10 MHz, all terminals retained and no matrix reductions. These
are wall times on the shared machine, not isolated-CI speed guarantees.
Preparation timings are one warm measurement; compute ranges cover three calls.
Allocation is cumulative Julia allocation, not peak memory. The inner `@time`
measurements exclude the surrounding reporting overhead.

| Case | Default coaxial compute | Boundary preparation | Compute with prepared boundary |
|---|---|---|---|
| 18 kV trefoil | 0.292–0.303 s; 12.168 MiB; 138.17k allocations | 10.674 s; 293.662 MiB; 748.77k allocations | 0.311–0.325 s; 12.176 MiB |
| 132 kV horizontal | 0.2839–0.2842 s; 4.978 MiB; 92.88k allocations | 2.670 s; 105.279 MiB; 524.32k allocations | 0.285–0.515 s; 4.955 MiB |

The former 30–40 s / 1.664 GiB warm 18 kV compute included the automatic
boundary model. The new default is deliberately the annular approximation,
not an accuracy-equivalent faster boundary solve. For the explicit boundary
model, the same-resolution preparation itself fell from about 32–40 s /
1.652 GiB to 10.674 s / 293.662 MiB. The former public reuse path in this
historical measurement was removed by the blueprint-ownership revision;
ordinary boundary compute calls now include blueprint reconstruction.

The new 18 kV boundary C differs from the previously measured tighter-integral
result by `4.39e-9` in relative matrix norm and `1.54e-8` in relative direct
core-to-jacket coupling. The latter is `1.633777251281124e-11 F/m`.
These compare numerical results at the same resolution, not against exact
physics. Resolution sensitivity of weak couplings remains as documented above.

First-call compilation remains substantial: the fresh 18 kV default compute
took 50.435 s, 99.30% reported as compilation. This refactor does not claim to
solve cold Julia/IDE latency. No Project.toml, Manifest.toml, dependency or
consumer-environment change was made. Restart Julia once to load the changed
struct definitions before rerunning an existing IDE script.

### Historical verification before the blueprint-ownership revision

Successful focused runs:

- `test/runtests.jl internal_shunt shunt_model`: 222 assertions, 8 items,
  116 s. Includes independent mathematical controls, audit invariance,
  explicit reuse, reductions, strict quadrature/budget failures, fallback,
  correlated sensitivities, UQ failure policy and serialization.
- Cable constants, formulation grids/records, dielectric laws, computation
  options, UQ options and reporting selection: 1,449 assertions, 18 items,
  432 s. No failures.
- Shared-source shunt UQ plus scientific publication serialization: 236
  assertions, 3 items, 120 s.
- `tag:quality`: 1,138 assertions, 6 items, 104 s. All passed, including
  external interface/import ownership and formula/display conventions.
- `tag:aqua`: 11 checks, 43 s. All passed. Its temporary-environment check
  required access to Julia's cache/logs and was rerun offline outside the
  read-only sandbox; no consumer environment was changed.
- Gauntlet `records_tests dielectric_consistency_tests`: 235 assertions,
  3 items, 71 s. Updated the automatic-boundary assumption to test default
  annular and explicitly selected boundary calculations separately.

The new files use the SciML formatter. Julia docstring conventions informed
the separation of physical units, material assumptions, numerical controls and
error guarantees in the API documentation. Concurrent lowercase literature-ID
edits and the user's existing inspection-script edits were preserved.

### Historical failure investigation and limits

Read the retained failure and request without overwriting any attempt. The
recorded request used 512 uniform draws and seed `0x132630`; its captured case
and benchmark source files still match the current files. The artifact does
not retain the failing draw, quadrature stage or error estimate.

A bounded replay checked the assembly and independent-audit logarithmic
integrals for the first **65** current fixed-seed realizations at the former
`rtol=1e-10`, with the original frequency normalization. All 65 completed;
the diagnostic replay was stopped before timing the new implementation.
The exact recorded failure was **not reproduced**, so neither a particular
roundoff mechanism nor exhaustion of all 100,000 evaluations is established.
New failure records contain the information needed to isolate the next such
integral, and strict/fallback behavior is tested with deliberately exhausted
numerical budgets.

The full 512-draw reference, all-references campaign, full ordinary suite,
coverage gate, documentation build and native FEM campaign were not rerun.
The existing retained result `jl_dt7ROH` and failed replacement `jl_HwDNwl`
remain unchanged. Numerical/reference acceptance is separate from this
structural refactor and the focused verification above.
