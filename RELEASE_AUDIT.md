# v0.2.0 release audit

Audit baseline: `main` at
`88f5de54202f5a97c7556b4ceac662b997f13de2`, before release-preparation
edits. The baseline test evidence was 1,750 package assertions and 19 FEM
integration assertions.

## P0 — release blockers

- **Unregistered hard dependency:** GetDP.jl was referenced through
  `[sources]`, which is not acceptable for General registration. Resolution:
  remove the dependency and vendor the compatible frontend privately from
  commit `b1d91b0d8974ea642b772462edcf6e26299fdf0a`, tree
  `2a2fe2782b1c17a235a63d3639b99ae012dab1f5`.
- **Optional systems loaded by core:** Makie, Gmsh, and GetDP were hard
  dependencies. Resolution: package extensions and explicit activation.
- **FEM coupled to routine coverage:** the external solver integration was
  discoverable by the normal test runner. Resolution: a manual release-only,
  coverage-disabled gate with a pinned GetDP archive.
- **Publication ambiguity:** two workflows could create releases. Resolution:
  TagBot is the sole tag and GitHub-release authority.
- **Registry metadata:** package version and installation guidance targeted a
  development checkout. Resolution: version 0.2.0 and clean-install checks.
- **Credential scan:** no private credentials were found. The only tracked
  token-shaped literal is the existing public Codecov badge identifier, not an
  upload credential. The repository history was not rewritten.
- **Third-party licensing:** the vendored Julia frontend retains its upstream
  BSD license and provenance. `THIRD_PARTY_NOTICE.md` also identifies the
  registered Gmsh wrapper and the licenses of the separately supplied Gmsh and
  GetDP programs. No GetDP binary or downloader is distributed.

## P1 — correctness and maintenance

- Aqua was skipped and reported three method ambiguities. The dispatch
  intersections are now explicit, and Aqua is unconditional.
- Three plotting exports had no definitions. They were corrected while moving
  implementation code to extensions.
- Binder, the TODO-to-issue scraper, the duplicate release template, and the
  accidental standalone cable-builder tree were maintained dead weight and
  were removed.
- The documentation build used `warnonly`, rewrote source files during an
  opaque build phase, and opened a browser locally. It was split into explicit
  instantiate, doctest, build, and deploy phases.
- Development instructions follows physical meaning and contains no “See also” sections.
- GetDP executable resolution could invoke a downloader. The compatibility
  frontend now checks `GETDP_EXECUTABLE`, then `PATH`, and otherwise errors.

## P2 — deliberately deferred

- Broad module or type redesign outside the dependency-extension boundary.
- Speculative performance changes and cleanup of nontrivial historical code.
- Simplification of the future FEM API beyond the compatibility facade.

## Branch evidence

The abandoned `updates-runner` tips recorded before deletion were:

- local: `47b53243b8353985665fb08e56389814d0da2410`
- origin: `56bc4ee0b3d987fc56c0b6e752180b6fa7158c0d`

Local and GitHub branch listings after cleanup contain only `main`,
`refactor/cablebuilder`, and `gh-pages`.

## Verification evidence

The following locally executable gates passed on Julia 1.12.6:

| Gate | Result |
| --- | --- |
| Routine package tests, including unconditional Aqua | 1,734/1,734 assertions; no plotting/FEM/formatter dependency resolved |
| CairoMakie extension tests | 1,744/1,744 assertions |
| Release-only FEM integration | 19/19 assertions with GetDP 3.5.0 |
| GetDP archive verification | SHA-256 `d3c28fa18f20d6147b4c7367d4dd802e9f7ddb58c608688bbb71919dbca8041d` |
| Independent Aqua run | all checks passed, including ambiguities, exports, stale dependencies, compat, piracy, and persistent tasks |
| Documentation doctests | passed |
| Strict documentation build | passed without `warnonly` |
| SciML formatting | JuliaFormatter 2.12.4 returned `true` with `overwrite=false` |
| Commit policy | scoped, lowercase, 72-character gitlint policy validated with gitlint 0.19.1 |
| Optional dependency boundary | core loaded no Makie backend, Gmsh, or GetDP; missing-backend errors were actionable |
| Registry dependency audit | all 32 dependencies and weak dependencies have registered or standard-library UUIDs |

The baseline and release trees were also compared structurally: every baseline
test expression remains represented, the 19 FEM assertions moved intact to the
manual integration harness, and two plotting lifecycle assertions were added.
The numerical totals of the routine runner changed because FEM discovery was
removed and Aqua was made unconditional, not because covered behavior was
discarded.

The prerelease-Julia job is defined in CI and must run on GitHub because only
Julia 1.12 is installed locally. The immutable release commit and its clean
commit-based installation result are recorded in the release handoff; a commit
cannot embed its own SHA.

---

## Ranking

- **R0 — release blocker:** fix, reject, or quarantine before registration.
- **R1 — should fix now:** high-value and generally inexpensive.
- **R2 — progressive:** safe to improve during later `0.2.x` releases if limitations are documented.
- **R3 — later/research:** substantial redesign or scientific validation.

Feasibility:

- **F1:** under half a day
- **F2:** roughly 1–2 days
- **F3:** roughly 3–5 days
- **F4:** 1–2 weeks
- **F5:** research or architectural work with uncertain duration

## R0 — must address before registration

| ID | Finding | Recommended v0.2 fix | Effort |
|---|---|---|---|
| SEC-1 | JSON deserialization calls `Meta.parse` and `Base.eval` on a file-controlled `__julia_type__`. A crafted JSON file can execute Julia code. See [deserialize.jl:53](/home/amartins/Documents/KUL/LineCableModels/src/importexport/deserialize.jl:53). | Replace dynamic type resolution with an explicit allowlist and schema version. Never evaluate type strings. | F2–F3 |
| SEC-2 | Public `.jls` loading invokes Julia’s native `deserialize`, which is unsafe for untrusted files and not a stable interchange format. See [cableslibrary.jl:155](/home/amartins/Documents/KUL/LineCableModels/src/importexport/cableslibrary.jl:155). | Remove `.jls` from automatic/public loading for v0.2, or quarantine it behind a deliberately named trusted-data API with a prominent warning. | F1 |
| DATA-1 | `Measurement` JSON round-tripping is catastrophically wrong: `10 ± 2` comes back as `2 ± 0`. The nominal value is discarded, and correlations cannot be preserved. See [deserialize.jl:115](/home/amartins/Documents/KUL/LineCableModels/src/importexport/deserialize.jl:115). | Safest bounded fix: explicitly reject serialization of `Measurement` values in v0.2. Full correlation-preserving serialization can follow later. | F1 now; F4 fully |
| DATA-2 | Library loading empties the destination before parsing and catches element failures, so malformed input can leave a partially loaded or empty library while appearing to succeed. | Parse and validate into a temporary collection, then swap only after complete success. Throw typed errors; never silently return partial data. | F2 |
| DATA-3 | Several exporters resolve relative paths against `@__DIR__`, potentially writing into the installed package source tree rather than the current directory. Some then swallow filesystem failures. | Resolve user paths relative to `pwd()`, create parent directories deliberately, and throw on failed export. | F1–F2 |
| HOST-1 | Merely importing the package replaces the process-wide Julia logger and changes its minimum level. Runtime probing confirmed that `using LineCableModels` changed `ConsoleLogger` into `TimestampLogger`. See [logging.jl:32](/home/amartins/Documents/KUL/LineCableModels/src/utils/logging.jl:32). | Remove logger mutation from `__init__`. Keep verbosity configuration opt-in and use `with_logger` internally. | F1 |
| SOL-1 | `ideal_transposition` has inverted and inconsistent behavior: impedance paths transpose when it is false, while one admittance path transposes when it is true. See [solver.jl:108](/home/amartins/Documents/KUL/LineCableModels/src/engine/solver.jl:108). | Correct the condition and add behavioral tests for Z/Y, Kron/non-Kron, and both option values. | F1–F2 |
| SOL-2 | The homogeneous-earth impedance documentation and integrand describe an integral over `0…∞`, but the implementation integrates only `0…1`. See [homogeneous.jl:184](/home/amartins/Documents/KUL/LineCableModels/src/engine/earthimpedance/homogeneous.jl:184). | Validate against a trusted published/reference implementation, then correct the integration domain and tolerances. Do not make this a blind one-character fix. | F2–F3 |
| SOL-3 | `EarthModel` does not retain its frequency grid. The solver checks only one vector’s length, so data sampled on different frequencies but with the same length are silently accepted. Homogeneous formulations also appear to use only the first two layers while the API advertises arbitrary multilayer earth. See [EarthProps.jl:141](/home/amartins/Documents/KUL/LineCableModels/src/earthprops/EarthProps.jl:141). | Store and validate the frequency grid and every layer. For v0.2, reject unsupported layer counts and narrow the documentation instead of attempting a full multilayer implementation. | F2 |
| SOL-4 | The primitive internal impedance matrix is stashed inside the cable loop. With multiple cables, later snapshots can include earlier earth contributions and omit later internal contributions. See [solver.jl:251](/home/amartins/Documents/KUL/LineCableModels/src/engine/solver.jl:251). | Move the snapshot to the correct completed-assembly boundary and add a multi-cable regression test. | F1–F2 |
| ENG-1 | `store_primitive_matrices=false` attempts to put `nothing` into fields typed as arrays. It fails immediately; the admittance path also writes `Pg` unconditionally. See [workspace.jl:137](/home/amartins/Documents/KUL/LineCableModels/src/engine/workspace.jl:137). | Make fields nullable or allocate only required storage; guard every optional write. Add a full compute regression test with storage disabled. | F1 |
| API-1 | `LineParameters` accepts mismatched Z/Y conductor dimensions. A modal slice is reconstructed as `PhaseDomain`, silently losing its domain. See [lineparams.jl:51](/home/amartins/Documents/KUL/LineCableModels/src/engine/lineparams.jl:51). | Enforce identical matrix dimensions and retain domain when slicing. Validate nonempty, finite, positive, ordered frequencies. | F1 |
| PHY-1 | Bare conductors reach singular potential/admittance calculations because their insulation coefficients become zero. This corresponds to open issues [#17](https://github.com/Electa-Git/LineCableModels.jl/issues/17) and [#18](https://github.com/Electa-Git/LineCableModels.jl/issues/18). | If correct bare-conductor support cannot be validated now, reject it with a clear error before matrix inversion. | F1 guard; F4 support |
| PHY-2 | `Fortescue` warns when transformed matrices are not diagonal and then discards the off-diagonal entries anyway. That silently approximates asymmetric or non-transposed systems. See [fortescue.jl:29](/home/amartins/Documents/KUL/LineCableModels/src/engine/transforms/fortescue.jl:29). | Refuse lossy diagonalization by default. Require an explicit approximation flag or return the complete transformed matrix. | F1–F2 |
| PHY-3 | Non-circular conductors are often reduced to undocumented circular/point approximations in GMD, GMR, sector capacitance, and sector conductance calculations. | For v0.2, clearly mark these paths experimental or reject unsupported shapes in stable calculations. Do not advertise validated sector/rectangular accuracy yet. | F1 guard; F5 validation |
| FEM-1 | Sector FEM geometry has a proven gap/unassigned-region defect, issue [#40](https://github.com/Electa-Git/LineCableModels.jl/issues/40). A candidate fix exists at commit `73c8ceb…`, but it is not integrated into this release branch. | Review and integrate the fix with an FEM regression, or explicitly exclude sector FEM from supported v0.2 functionality. | F1–F2 |
| SCI-1 | The test suite has extensive formula-level checks but very little independent validation against published reference cases or another trusted solver. Tests that reproduce the implementation’s own equation cannot detect a wrong equation or integration interval. | Add at least one reference benchmark for the default cable solution, earth impedance/admittance, transposition, and every shape claimed stable. | F3–F5 |
| CI-1 | PR #47’s internal checks pass, but the Codecov project gate fails: approximately 39.96% project coverage versus 48.23% previously. Patch coverage is about 76.54%, with 175 changed lines uncovered. | Cover the newly introduced critical paths. Do not solve this merely by lowering the project threshold. | F2–F3 |
| REL-1 | The current routine test count is **1,734**, not the planned baseline of 1,750. FEM’s 19 assertions have earlier passing evidence but were not rerun during this read-only audit. | Reconcile the missing 16-test discrepancy and rerun the manual FEM gate against the final release SHA. | F1 |

My pragmatic interpretation of these blockers is important: several do not need a grand implementation now. A clear rejection of unsupported data or geometry is much safer than silently returning plausible but wrong numbers.

## R1 — easy or high-value fixes for this release

### Numerical and model validation

- `calc_tubular_gmr` returns `Inf` for a tiny but nonzero inner-radius ratio instead of approaching the finite solid-conductor limit. **F1**
- `calc_equivalent_mu` can divide by zero or become badly conditioned when `r_in ≈ r_ex`. **F1**
- `calc_tubular_inductance` appears to calculate annular/coaxial field inductance while being documented as tubular-conductor internal inductance. Clarify or deprecate the name rather than silently changing the equation. **F1**
- Use factorizations/linear solves instead of explicit inverses in Kron reduction and the Levenberg transform. See [reduction.jl:74](/home/amartins/Documents/KUL/LineCableModels/src/engine/reduction.jl:74). **F1–F2**
- Validate line length, coordinates, material values, radii, thicknesses, phase labels, connectivity length, and finite values at public constructors. Negative phase numbers currently enter the “grounded/non-phase” path. **F1–F2**
- Require phase labels to be contiguous `1:N`, or normalize them. Labels such as `{1,3}` currently produce an obscure assertion. **F1**
- Replace public-boundary `@assert` checks with `ArgumentError` or `DomainError`. **F1**
- Remove `@inbounds` from public `LineParameters` indexing. Bounds mistakes should remain safe. **F1**
- Add optional validation for finite, reciprocal, passive Z/Y matrices rather than silently accepting physically impossible results. **F2**

### Cable geometry and builders

- `CircStrands(n=1, r_in>0)` ignores `r_in`; compressed layers are not area-preserving; physical fill validation is not consistently applied. See issue [#32](https://github.com/Electa-Git/LineCableModels.jl/issues/32). **F1–F2**
- `RectStrands` is effectively unfinished: almost no direct tests, inconsistent physical/effective thickness, and ParametricBuilder silently constructs circular strands for any `AbstractStrandsLayer`. See [cablebuilderspec.jl:124](/home/amartins/Documents/KUL/LineCableModels/src/parametricbuilder/cablebuilderspec.jl:124). Quarantine it as experimental unless fixed. **F1 guard; F3 repair**
- The first insulator specification’s `n_layers` is ignored by ParametricBuilder. **F1**
- Automatic cable spacing starts at roughly one radius instead of two radii, producing overlapping candidates that are later discarded. **F1**
- `get_outer_radius` assumes the last component is outermost. The model does not enforce radial ordering or non-overlap even though the solver depends upon it. **F2**
- Sector geometry hardcodes `Point2f`/Float32 despite generic numeric types. Its fixed 30-point arc approximation makes results resolution-dependent. See [sector.jl:232](/home/amartins/Documents/KUL/LineCableModels/src/datamodel/sector.jl:232). Mark experimental now; generalize progressively. **F1/F4**
- Formation helpers need positive-spacing validation and clearer trifoil radius semantics. **F1**

### Iterator and error semantics

- Cable and system builder iterators catch task exceptions and turn them into end-of-iteration or skipped designs. This can silently produce partial/empty design sets. **F1–F2**
- SystemBuilder decides whether to skip failures by string-matching error messages such as “overlap.” Use typed exceptions. **F2**
- `length(SystemBuilderSpec)` reports an upper-bound cardinality even though invalid combinations are skipped. Use `SizeUnknown` and expose a separate upper-bound function. **F1**
- Several `add!` methods return a new promoted object rather than mutating the original. A caller who ignores the return silently loses the update. For v0.2, reject incompatible element types rather than unexpectedly changing bang semantics. **F1**
- `delete!` should return the collection, consistent with Base. **F1**
- Remove logging from ordinary constructors, lookups, and `Base.get`. Library operations should not emit informational noise. **F1**

### UQ correctness

- UQ calls `Random.seed!`, changing the application-wide RNG. Accept an `rng` and use a local seeded generator. See [montecarlo.jl:46](/home/amartins/Documents/KUL/LineCableModels/src/uq/montecarlo.jl:46). **F1–F2**
- `_rand_in` asserts for `lo == hi`, mishandles empty collections, and ignores the third tuple element in some contexts. Issue [#31](https://github.com/Electa-Git/LineCableModels.jl/issues/31). **F1**
- `(lo, hi, n)` means a deterministic grid in builders but a mean/standard-deviation construction in Monte Carlo. This is an API trap. Introduce explicit range and distribution specification types progressively; document the current distinction immediately. **F2/F4**
- Normal sampling can generate negative radii, thicknesses, and physical properties. Invalid samples may then be skipped, biasing the output distribution. Add domains and deliberate rejection/resampling. **F2–F3**
- Validate `ntrials`, confidence, tolerance, bin count, and print interval. A single cable trial currently produces `NaN` uncertainty. **F1**
- PDF construction should require strictly increasing bins and finite, nonnegative densities and handle degenerate samples. **F1–F2**
- The joint `Measurement` surrogate allocates a dense trial-by-output matrix and one latent measurement per trial; this can become enormous. Make it opt-in or use a reduced-rank representation. **F3**

### API and dependencies

- Remove the unused `ForceImport` dependency/import. Aqua does not flag it because merely importing it counts as usage. **F1**
- `Colors` is loaded in core even though color handling belongs to plotting/FEM extensions. Move the import out of core. **F1**
- Remove the dead `.lis` generic/commented parser, issue [#12](https://github.com/Electa-Git/LineCableModels.jl/issues/12). **F1**
- Remove the dead `nonsensify` function and its unsuitable public-facing docstring. **F1**
- Remove the global `Base.show(::Dict{String,Material})` specialization; it changes display for ordinary dictionaries unrelated to the package. **F1**
- Make `CablesLibrary` and `MaterialsLibrary` obey the same collection interface. **F2**
- `Thickness` and `Diameter` subtype `Number` without implementing the expected numeric interface. Do not redesign it now, but at least reject nonfinite inputs and document that they are wrappers, not arithmetic scalar types. See [types.jl:1](/home/amartins/Documents/KUL/LineCableModels/src/datamodel/types.jl:1). **F1**
- Validate `ResultsView` option symbols instead of treating unknown `per` and `mode` values as defaults. **F1**
- Fix complex-value display producing forms like `real+-imagim`. **F1**

### Plotting extensions

The Makie extension split itself is sound: core import does not load Makie or Gmsh.

Remaining corrections:

- Backend extensions mutate one global backend selection during `__init__`; loading multiple backends makes the winner order-dependent. **F1–F2**
- `with_backend` cannot reliably restore a prior `:none`/unknown state. **F1**
- Calls with no active backend can report “unknown backend” instead of an actionable “load CairoMakie/GLMakie/WGLMakie first” error. **F1**
- Plot value-expression parsing uses `eval` for index expressions. It is less exposed than JSON loading but unnecessary; parse integer, colon, and range syntax structurally. **F1**

### Documentation and presentation

- There are roughly **125** `julia` code fences but only **19** executable `jldoctest` fences. Most examples can still go stale without CI noticing. Convert supported self-contained examples; remove placeholders and fixture/external-program examples. **F3**
- The mathematical-notes policy is not yet satisfied. There are approximately 46 `# Notes` sections, while central scientific files—earth kernels, insulation formulations, Bessel wrappers, Fortescue, and Levenberg—lack adequate equations, assumptions, units, validity domains, and references. **F3**
- An EarthModel example attempts to mutate an immutable `EarthLayer`, but it is hidden in a non-executed Julia fence. **F1**
- Fix unit errors: helical overlength is dimensionless, and shunt conductance per unit length is `S/m`, not `S·m`. **F1**
- Default material values need references and validity conditions. Values such as steel permeability are highly condition-dependent. **F2**
- Clean remaining jokey/profane implementation comments. This is a tiny improvement that materially affects public presentation. **F1**
- `CONTRIBUTING.md` should mention the complete test matrix, FEM release test, and security reporting. Add `SECURITY.md`. **F1**
- Validate `CITATION.cff` in CI. **F1**

## R2 — progressive hardening

These should not hold up `v0.2.0` if the affected functionality is correctly delimited:

- Full JSON schema migrations and correlation-preserving uncertain-data serialization. **F4**
- Explicit distribution/range types throughout UQ. **F4**
- Scale-aware geometry tolerances instead of fixed absolute tolerances. **F2**
- Independent uncertainty tests for the complex Bessel implementations, including analytic derivatives and correlations. **F2**
- Robust mode tracking and conditioning diagnostics in the Levenberg transform. **F4–F5**
- Shape-specific GMD/GMR calculations for RectStrands, Sector, Strip, and Tubular conductors. **F5**
- General numeric geometry without Float32 coercion. **F4**
- Resolve the `add!` promotion design with a clean nonmutating `add` API in a future breaking release. **F3**
- Clarify or eventually rename `compute!`, which returns new results rather than mutating its formulation/problem argument. **F3**
- Pin GitHub Actions to immutable commit SHAs and add automated update tooling. Current workflows mostly use version tags. **F2**
- Reduce docs workflow permissions; deployment permissions do not need to be granted to every build step. **F1**
- Reconcile local and remote `gh-pages` tips through the normal docs deployment. **F1**
- Move showcase-heavy assets out of the package if install footprint becomes important. The current source archive is only about 6 MB, so this is not urgent. **F2**

## R3 — deliberately later

- A fully validated arbitrary multilayer-earth implementation rather than the proposed v0.2 guard. **F5**
- Formal scientific validation of Levenberg normalization and mode continuity over difficult frequency sweeps. **F5**
- Production-grade rectangular and sector-strand electromagnetic models. **F5**
- Full support for bare, overhead, and mixed cable/line configurations. **F5**
- Redesigning `Thickness`, `Diameter`, and other wrapper hierarchies. **F4**
- Broad naming, signature, bang-function, export, and module-layout normalization. This would become exactly the refactor you said not to do. **F4–F5**
- Git-history cleanup. The current archive is small, but `.git` is roughly 472 MB due to historical notebooks/assets. Do not rewrite published history as part of this release. **F4 and disruptive**
- Wishlist functionality such as insulator strips ([#29](https://github.com/Electa-Git/LineCableModels.jl/issues/29)) and broader JLS support ([#3](https://github.com/Electa-Git/LineCableModels.jl/issues/3)). **F4–F5**

## Things that are already in good shape

This audit was not uniformly negative:

- The routine suite passes: **1,734/1,734**.
- Aqua passes, including ambiguities, undefined exports, and project checks.
- Julia 1.12 and prerelease CI pass.
- Cairo plotting, formatting, gitlint, docs, and clean-install jobs pass.
- Core loading does not activate Makie or Gmsh extensions.
- Gmsh and Makie are now optional dependencies through Julia extensions.
- The GetDP frontend is privately vendored with BSD license and provenance, with no binary downloader.
- GetDP executable resolution uses `GETDP_EXECUTABLE`, then `PATH`.
- Binder, accidental standalone cablebuilder, and obsolete automation are gone.
- “See also” sections are gone.
- Documentation builds strictly rather than hiding warnings through `warnonly`.
- Licensing, `CITATION.cff`, version `0.2.0`, TagBot ownership, and registry dependency metadata are broadly suitable.
- The current tree heuristic did not find an obvious credential. I did not complete a dedicated full-history secret scan, so one should still be run before final handoff.

## Recommended bounded route to v0.2.0

I would do this in four passes:

1. **Safety and state integrity:** JSON evaluation, JLS quarantine, Measurement rejection/fix, transactional loading, output paths, global logger.
2. **Core correctness:** transposition, earth integral, frequency/layer validation, Zin snapshot, primitive-storage option, LineParameters invariants, bare-conductor guard, Fortescue behavior.
3. **Scope control:** integrate sector FEM fix or disable that path; mark RectStrands, sector analytic approximations, and Levenberg experimental where they cannot be validated now.
4. **Evidence:** add independent scientific benchmarks, restore Codecov legitimately, reconcile the 1,734/1,750 test count, run the 19 FEM assertions, clean-install from the exact final SHA.

With pragmatic guards instead of implementing every unfinished model, this looks like roughly **one focused week plus however long independent earth/solver validation takes**. Fully validating sector, rectangular, multilayer-earth, and Levenberg functionality is a separate research-grade project.

The audit was read-only. I did not modify either of the already-dirty files, [RELEASE_AUDIT.md](/home/amartins/Documents/KUL/LineCableModels/RELEASE_AUDIT.md) or [docstrings.md](/home/amartins/Documents/KUL/LineCableModels/docs/src/docstrings.md). The Unslopifier quality gate could not run because this repository is not enrolled and has no `.unslopifier.toml`; that was diagnosed as missing setup, not counted as a package failure. The Julia docstring criteria were used to inspect mathematical notes and examples, and the GitHub audit supplied the issue and CI findings above.

---

You’re right. I misclassified the architectural situation.

`refactor/cablebuilder` is the future implementation, while its appearance on `main` was the accidental leak. The branch itself is not legacy. More importantly, I audited the current code as though every subsystem was intended to mature in place. That caused me to emphasize defects in code you already intend to replace while missing the much more consequential duplication in plotting and parametric/UQ.

The revised conclusion is:

> Stabilize the reusable numerical kernel, quarantine the outgoing architecture, and stop investing in abstractions that CableBuilder/Gridspace will replace.

## Revised architecture map

| Category | Components | Treatment |
|---|---|---|
| Durable kernel | Engine mathematics, impedance/admittance formulations, transforms, `LineParameters`, uncertainty-aware numerics | Fix and independently validate |
| Durable infrastructure | Makie extension boundary, `UnitHandler`, backend handling, CI/docs/release machinery | Clean and retain |
| Transitional compatibility | Current DataModel, MaterialsLibrary, CablesLibrary, ImportExport, ParametricBuilder, current UQ grammar | Deprecate, constrain, avoid redesign |
| Effectively abandoned | FEM/GetDP, sector FEM and much sector-specific machinery | Quarantine, warn, retain only recoverable knowledge |
| Future architecture | `refactor/cablebuilder`, Grid/Gridspace grammar | Develop separately and eventually replace transitional layers |

This demotes a large portion of my previous audit:

- `Thickness`, `Diameter`, collection interfaces, proxy types: **do not fix architecturally**.
- MaterialsLibrary and CablesLibrary API elegance: **do not polish**.
- Full serialization schema/correlation architecture: **do not build**.
- RectStrands/Sector integration with the old DataModel: **do not complete**.
- ParametricBuilder iterator contracts and old builder composition: **do not refactor in place**.
- FEM sector correctness: reject unsupported use rather than rehabilitating it.

One exception remains: deprecated public code must not be dangerous. The JSON `eval` issue does not become harmless because the loader will eventually disappear. The pragmatic fix is now to disable/remove that loading path, not redesign it. Likewise, silent data corruption should become a clear unsupported-operation error.

## Plotting: this is now one of the main architectural findings

Your diagnosis is exactly right. There are not merely two generations of plotting code; there are two partially interwoven engines.

### The old path

The current working LineParameters plotter consists primarily of:

- [plotmetadata.jl](/home/amartins/Documents/KUL/LineCableModels/src/engine/plotmetadata.jl): 385 lines of its own prefix/unit metadata system;
- [plot.jl](/home/amartins/Documents/KUL/LineCableModels/src/engine/plot.jl): 797 lines for LineParameters preparation, rendering, error bars, widgets, log scaling, limits, legends, and export;
- [UQ plot.jl](/home/amartins/Documents/KUL/LineCableModels/src/uq/plot.jl): 694 more lines;
- [PlotUIComponents.jl](/home/amartins/Documents/KUL/LineCableModels/src/plotbuilder/plotuicomponents/PlotUIComponents.jl): another 819-line UI framework.

It works and contains the visual behavior you like, but mixes:

- semantic quantity selection;
- unit conversion;
- data slicing;
- plotting;
- UI state;
- backend state;
- export behavior;
- uncertainty rendering.

That explains both why the output is good and why maintaining it feels awful.

### The new PlotBuilder

The declarative design is genuinely better:

```text
domain object
    → normalized plot grammar
    → RenderSpec
    → backend/UI renderer
```

The separation of `AxisSpec`, `SeriesSpec`, `ViewSpec`, `PageSpec`, and `RenderSpec` is conceptually sound, and `UnitHandler` is a much better unit authority than `plotmetadata`.

But it is not finished:

- The only concrete declarative grammar is `MCStatsPlotSpec`.
- There is no declarative LineParameters spec.
- Nothing in production calls `make_render`.
- The working `plot(::LineParameters)` still uses the old path.
- The new renderer supports only line, scatter, and heatmap primitives.
- Export code still calls an obsolete `render` function after it was renamed to `build`.
- The new `UIComponents` still imports a type from the old UI stack.
- Comments explicitly describe work to be completed “in the next iteration.”
- There are no direct PlotBuilder grammar, parity, or LineParameters plotting tests.
- Cairo CI tests cable previews, not LineParameters or UQ plotting.

So your parity concern is justified: parity has not been demonstrated, and the code strongly indicates it was never completed.

## Recommended plotting strategy

Do not ship both engines indefinitely.

The pragmatic path is a single narrow vertical migration centered on the durable output type: `LineParameters`.

1. Write a behavior inventory from the old plotter:

   - R/L/C/G;
   - real/imaginary and magnitude/angle Z/Y;
   - phase/modal indexing;
   - selected conductor pairs;
   - Measurements error bars;
   - linear/log axes;
   - unit overrides;
   - legend behavior;
   - SVG export;
   - Cairo/GL/WGL behavior.

2. Implement one declarative `LineParametersPlotSpec` using `UnitHandler`.

3. Test the generated `RenderSpec` without Makie:

   - exact numerical series;
   - unit factors;
   - labels;
   - selected indices;
   - grouping;
   - domain behavior.

4. Port only the UI behavior required by that spec into the new renderer.

5. Add Cairo smoke tests for rendering and export.

6. Wire the public `plot(::LineParameters)` through:

```julia
make_render(LineParametersPlotSpec, lp; kwargs...) |> build
```

7. Delete:

   - `plotmetadata.jl`;
   - old LineParameters spec preparation/rendering;
   - the old `PlotUIComponents` pipeline;
   - duplicated prefix/unit dictionaries;
   - generic old rendering helpers.

8. Do **not** port all current UQ plotting merely for parity. The current UQ result architecture is also transitional. Retain only plots whose result types are expected to survive the Gridspace migration.

This is probably **4–8 focused days**, not a trivial cleanup—but it is bounded because UQ plot parity is no longer part of the target.

The faster release alternative is the reverse: keep the working old plotter, consolidate its units onto `UnitHandler`, and remove/park the unfinished PlotBuilder from v0.2. That is likely **1–2 days**. It is less architecturally satisfying, but still much better than publishing both frameworks as though both were supported.

My preference, based on what bothers you, is the narrow PlotBuilder vertical slice. The important constraint is to implement it from one real use case and delete unused traits as you go. The current PlotBuilder is at risk of becoming a generic plotting framework designed before its first complete plot.

## Parametric/UQ: Gridspace changes this from “repair” to “replace”

The existing stack is approximately several thousand lines of:

- permissive tuple interpretation;
- deterministic enumeration;
- uncertainty interpretation;
- constructor reconstruction;
- nested builder handling;
- Monte Carlo sampling;
- validity filtering;
- result reduction.

Gridspace correctly observes that the first four are mostly the same operation:

```text
normalize every field into a Grid
    → Cartesian product
    → construct target
```

That is a much cleaner grammar.

The implementation in PowerBlocks and `refactor/cablebuilder` proves the architectural idea. It is small, compositional, and naturally supports nested builders. The Measurements extension is also the correct direction: deterministic combinatorics remains core, while uncertain iteration is opt-in.

However, Gridspace does not replace all of UQ. It should replace the **parameter-space grammar**, not become the solver, scheduler, constraint engine, and statistics package.

The clean separation should be:

```text
Grid / Gridspace
    parameter meaning and enumeration
             │
             ▼
CableBuilder
    validated domain-object construction
             │
             ▼
sweep / sample
    execution, RNG, rejection, parallelism
             │
             ▼
summarize
    moments, covariance, PDFs, retained samples
```

That would remove the angry tuple parser without recreating it inside another abstraction.

### Small Gridspace issues to fix before transplanting

The concept is strong, but PowerBlocks’ implementation is lightly tested:

- no explicit `rng` argument;
- `extrema(RelativeGrid)` is wrong for negative values;
- distribution selection treats Normal specially and effectively maps everything else to Uniform;
- macro expansion/default handling has no direct tests;
- nested Gridspaces, iteration order, inference, cardinality, and constructor failures are barely covered;
- the CableBuilder comment says scalar keyword input returns a concrete object, but `@gridspace` currently wraps every field and returns a `Gridspace`;
- physical-domain rejection/resampling still needs a deliberate policy;
- independent versus correlated uncertain inputs need explicit semantics.

These are F1/F2 repairs. They do not undermine the architecture.

### Migration feasibility

- Replacing tuple parsing and deterministic combinatorics: **F2, 1–2 days**.
- Connecting future CableBuilder constructors to Gridspace: substantially already conceptualized, **F2–F3**.
- Replacing the complete UQ execution/result stack: **F3–F4**, because RNG, failed samples, covariance, PDFs, parallel evaluation, and reproducibility remain.
- Rebuilding the outgoing DataModel around it: unnecessary.

One caution: I would not automatically turn frequency into another Grid dimension. Frequency is a solver analysis axis with shared matrix shapes and potential reusable workspace. Gridspace should initially enumerate designs and uncertain parameters, while the solver continues handling its frequency vector internally. Change that boundary only after benchmarks demonstrate that pointwise frequency evaluation is cleaner or faster.

## Revised pre-registration priorities

### Still critical

- Fix the global logger mutation.
- Fix `ideal_transposition`.
- Validate and fix the earth-impedance integration domain.
- Fix the multi-cable primitive-matrix snapshot.
- Fix `store_primitive_matrices=false`.
- Fix `LineParameters` dimension/domain invariants.
- Add independent reference cases for durable physical formulas.
- Disable unsafe/destructive persistence rather than redesigning it.
- Decide which plotting engine v0.2 actually supports; do not advertise both.
- Run the final routine and FEM compatibility gates.

### Transitional containment

- Mark the current DataModel, ParametricBuilder, UQ grammar, FEM, and sector paths as transitional/experimental as appropriate.
- Remove them from aspirational “stable architecture” documentation.
- Do not invest in collection purity, proxy-type elegance, or old serialization design.
- Make unsupported operations fail clearly instead of silently producing incorrect results.

### Future work

- CableBuilder/Gridspace replaces DataModel construction and ParametricBuilder.
- A separate execution/reduction layer replaces the current UQ monolith.
- Durable BaseParams equations migrate with reference tests.
- The solver consumes the new immutable physical representation.
- PlotBuilder consumes durable result containers, not the outgoing builder/DataModel.

So yes: this materially changes the audit. The largest architectural problems are now the dual plotting implementation and the obsolete parametric grammar—not the quality of the old DataModel abstractions. The old architecture needs a safe exit, not rehabilitation.