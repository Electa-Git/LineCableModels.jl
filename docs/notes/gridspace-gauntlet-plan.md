# Gridspace formulation comparisons — approved implementation plan

Inspected checkpoint: `0d464f95`, `release/v0.2.0-pscad-ext`. Approved by the user for end-to-end implementation in this worktree. The scope includes formulation comparisons, reporting, plotting, mutable drafts, accepted bundles and explicit packaging/publication commands.

## 1. The public workflow

One case, one explicit reference calculation, and one candidate formulation or formulation Gridspace form a benchmark. A formulation point means the complete calculation choice: backend, every formula family, scalar or air/earth/mixed selection, physical overrides, equivalent-earth choices and numerical settings.

The default workflow evaluates the reference once and the candidate Gridspace through the existing combinatorial `compute` overload. It compares every candidate with that reference. Reports group by complete formulation and the established frequency bands. This computes and retains comparison data; it does not open windows or generate figures. Once a user explicitly requests a plot for a selected case, that plot overlays all selected candidate formulations and the reference in each corresponding matrix cell.

Grouping does not average different formulations, merge equal curves, or reduce matrix terms without explicitly labelling a summary statistic. The reference determines comparison direction; its backend does not determine scientific authority.

## 2. Current evidence and defects

- `src/parametricbuilder/traversal.jl:226` already sends `compute(problem, formulation_space)` through `ParametricProblem` and `Combinatorial`.
- `src/parametricbuilder/results.jl:128` retains separate problem and formulation axes. `results[p, f]` already expresses the required coordinates, with problem index varying fastest in linear storage.
- `gauntlet/comparisons/saved.jl:113` still owns the collection-comparison traversal and reduces its coordinates to flat reference/candidate indices.
- `gauntlet/reporting.jl:113` forces selection of a single compared pair. `test/gauntlet/result_space_tests.jl` even requires the default plot to throw for multiple pairs. That expectation must be replaced.
- `ext/LineCableModelsMakieExt/LineCableModelsMakieExt.jl:274` and its matrix recipe already overlay multiple scalar results. Their generic labels are `Result 1`, etc.; there is no formulation-aware ParametricResult recipe.
- `BenchmarkTableDefinition` currently contains only `clip`. Band requests live in Gauntlet, while detailed benchmark table construction also lives there. The general ReportBuilder therefore cannot provide the same analysis outside Gauntlet.
- `details(result).formulations` already exists for analytical, PSCAD and FEM results, but the records differ. The analytical engine deliberately bounds those metadata types so a formulation Gridspace can yield one consistent concrete LineParameters type.
- Scalar and `(air, earth, mixed)` earth selections already exist in this worktree and in the inspected main worktree. Their applicability still belongs to indexed formula dispatch.
- `gauntlet/campaigns.jl` currently rejects a second `run_campaign` in an existing directory; `resume_campaign` permits reuse only with unchanged declarations and runtime sources. This does not yet provide the proposed replaceable draft workflow.
- The CLI exposes `lock`, but packaging and artifact binding are separate `package.jl`/`bind.jl` scripts. `package_collection` currently reads collection staging directly, without requiring a locked bundle, and permits local replacement through `force`. `.github/scripts/package_gauntlet.py` also retains the former `test/gauntlet` paths. These routes do not establish one acceptance-to-publication sequence.

## 3. Ownership

| Owner | Responsibility |
|---|---|
| Formula/backend owners | Expose declared selections and applied equation information, including triples, overrides and numerical choices. |
| ParametricBuilder | Materialize Gridspaces once; preserve problem/formulation coordinates through numerical results and comparisons. |
| Engine `compare` | Per-term RMS arithmetic, quantity-specific zero tolerances, band-to-sample selection and numerical comparability checks. |
| ReportBuilder | Default report request, five-band order, formulation grouping, complete error tables, summary tables and optional illustration. |
| PlotBuilder with the Makie extension | Formulation-aware selection, matrix-cell overlays, shared legend descriptions, stable series styling and existing UI/export behavior. |
| Gauntlet | Case/reference/candidate declarations, execution, recovery, storage, artifact locking and invoking the shared reporting operations. |

Gauntlet will not own a second formulation enumerator, band definition, legend formatter, RMS grouping algorithm or matrix renderer.

## 4. Calculation and coordinate contract

Preserve existing `compute`, `Combinatorial`, `ParametricProblem`, Gridspace product/zip semantics and `ParametricResult` indexing. A scalar candidate has one candidate formulation; it follows the same reporting rules without an extra scalar wrapper type.

Keep the full formulation selections on the retained formulation axis. Do not insert formula functors or closures of different concrete types into every `LineParameters.details`: that would break the current collection requirement that numerical results share a concrete type. Existing compact applied-equation records remain bounded; rich selection records belong to the axes/publication.

Add methods of existing `compare` for a scalar reference and a ParametricResult. For one quantity and one band, the result is an existing `ParametricResult` containing `RMSError` values, with the candidate's problem/formulation axes preserved. The reference identity is retained in the surrounding comparison publication. This does not add a nested result-space envelope or a new comparison-space type; the existing element invariant accepts a concrete RMSError.

Multiple problem points remain separate groups/pages. They are not flattened into unrelated formula curves. A deliberately supplied scalar reference can be a fixed baseline where physical output coordinates agree. Multiple reference points require explicit pairing; equal collection lengths alone must not invent a pairing. Existing saved explicit pairings remain readable.

Validate actual output coordinates: terminal identities and order, basis, domain, matrix dimensions and frequency coordinates for RMS. Equality of matrix size is insufficient. Use the retained case/problem/output mapping, including explicit reductions, rather than declaring two inputs identical by label. Different material laws or numerical models are legitimate comparisons.

A naked externally constructed LineParameters may lack terminal identities or a formulation record. It must not acquire invented identities from its type or dimensions. The standalone API accepts explicit named calculation context when that information is absent. Native result owners and Gridspace axes supply it when available. Plot overlays may use different recorded frequency grids; an RMS request still requires matched samples unless a separate explicit conversion has been performed.

## 5. Complete formulation identity and descriptions

Use one explicit record contract supplied by the owners, extending existing `NamedTuple` methods. It distinguishes:

- backend and complete requested formula selections;
- scalar defaults versus explicit air/earth/mixed selections;
- actual equation tags/cases used for the resolved problem;
- parameters and callable overrides;
- equivalent-earth rule/order and frequency/temperature constitutive selections;
- numerical choices, including integration method and controls;
- matrix reduction choices and the separate execution record.

The actual source/target layer indices and kind remain equation grammar. Display code never selects an equation or fills an unsupported case. A true multilayer formula keeps its identity and indexed cases; it is not converted into an air/earth/mixed triple.

Live axes and recovered axes must expose the same declared record schema. Normalize serialized formats once in their readers; remove live/report `hasproperty` guesses and shape-dependent default substitutions. Do not deserialize callbacks or load a solver merely to identify an old curve.

Extend the existing public `description` operation for these formulation records, with one documented compact rendering contract shared by reports and plots. Each curve retains its original formulation index and full record. Compact labels display the differences within the selected formulation collection; common settings appear once in its accompanying configuration table/caption. User labels remain an explicit display override.

For example, a combined selection can be described as:

```text
F2 · earth Z: air=Carson1926; earth=Pollaczek1926; mixed=Lucca1994
```

If earth Y, insulation behavior, constitutive laws, EHEM or integration choices also vary, the description identifies those differences too. Two `:default` selections using different integration methods remain distinct. Duplicate display labels are disambiguated with original indices; curves are never deduplicated because their samples agree. `:default` remains `:default`.

A numerical result's recorded applied choices and its declared full selection must be distinguishable. An unused branch of a triple remains part of the declared combination, and is not falsely described as an equation evaluated for the present geometry.

## 6. ReportBuilder owns the report request

Extend the existing `BenchmarkTableDefinition` with a normalized NamedTuple of comparison controls and presentation/illustration choices. No replacement policy type or comparison manager is introduced. The default line-parameter request is:

```julia
BenchmarkTableDefinition(;
    quantities=(Z, Y, R, L, G, C),
    bands=(:all, :dc, :harmonic, :narrow, :wide),
    normalizations=(:reference_rms,),
    fundamental=50.0,
    harmonics=50,
    clip=false,
    illustration=nothing,
)
```

This is proposed syntax. Public quantities use the package's existing observation grammar; serialized names are reader/writer details. Explicit tolerances and unsupported-quantity explanations retain the existing numerical meaning. Pointwise normalization is optional.

The report owns these default groups and their order. Engine remains the single owner of how each band selects stored samples. Retain the current numerical definitions for this work:

| Report group | Selector | Current requested range |
|---|---|---|
| Entire range | `:all` | All stored frequencies |
| Near DC | `:dc` | 0.1–100 Hz |
| Harmonic | `:harmonic` | `fundamental` through `fundamental * harmonics` |
| Narrowband | `:narrow` | 1 kHz–1 MHz |
| Wideband | `:wide` | Strictly above 1 MHz |

These groups overlap intentionally. Harmonic range means the stored samples within the range, not only integer harmonics. Existing ordinary endpoint snapping is retained, and requested bounds, actual bounds and sample indices are displayed. Empty bands retain an explicit no-samples result. No frequency-grid change or interpolation is part of this plan.

Gauntlet retains this report request instead of separately merging its own comparison defaults. Existing `validate` methods check the request and all known coordinate restrictions before solver execution.

The report follows its existing choreography:

```text
select → tabulate → illustrate → encode → write
```

For an explicit request over raw reference/candidate results, `select` computes the requested comparisons once through `compare`. For a completed comparison publication, `select` validates and selects the retained products without calculating RMS again. A request for an unrecorded band requires explicit reanalysis of raw operands; loading or changing display units must not silently create it.

Extend the existing `ReportArtifact` with its `published` value. That is the scientific data already produced by `select`, before display clipping, scaling or table reshaping. Its role is concrete: Gauntlet stores that exact comparison publication, and later tables/plots consume it. The current `.table`, `.illustration` and `.output` remain. This avoids reconstructing numerical truth from formatted DataFrames or calling internal report stages from Gauntlet.

The publication uses existing NamedTuples, result spaces, RMSError and observations. It retains reference/candidate calculation records, their numerical axes, the request and comparison products. Every scalar/space distinction is handled by methods of existing operations.

Tables expose the following views when requested; only the summary view is included in the standard Literate page and ordinary report display:

- formulation/calculation records keyed by original problem/formulation indices;
- complete RMS matrices by band, quantity, normalization and formulation;
- all matrix terms with physical row/column labels, units, status and reason;
- an explicitly labelled summary by case, band, formulation and quantity.

Absolute and relative RMS remain separate. Numerical-zero G retains absolute RMS and unavailable relative RMS with its reason. Summaries can show the largest per-term relative/absolute discrepancy and its coordinates; those are maxima of per-term metrics, not an RMS of the whole matrix. No default averaging across formulations or cases is introduced.

## 7. PlotBuilder owns the Gridspace-result recipe

Add methods to the existing PlotBuilder `plot` operation in the Makie extension:

```julia
plot(results, (Z, Y))
plot(results, (Z, Y); reference=reference_result)
plot(report_artifact)
```

These are proposed overloads. `results` is a ParametricResult of line parameters, not an unevaluated Gridspace. Plotting never computes a formulation. It obtains full selections from the retained axes/publication and invokes the existing matrix-cell renderer.

Plotting is explicit. For a report spanning multiple cases or problem points, the caller must select a case/problem or explicitly request a list of them. A bare `plot(report_artifact)` is convenient when that artifact contains one problem; it must not silently open dashboards for an entire campaign. Multiple selected problems get separate page groups with explicit problem context.

Default behavior within an explicitly requested plot for one problem:

- include every selected formulation and the reference once;
- produce the existing R/X dashboards for Z and G/B dashboards for Y;
- put the same terminal pair in the same cell in every series;
- show both off-diagonal terms;
- use the shared complete-formulation descriptions in legends;
- keep each formulation's style stable across quantities, bands and filtering by using its original index;
- preserve the current palette, unit options, legend controls, visibility controls, zoom and SVG export.

Direct R/L/G/C requests continue to work. Selection of formulation indices is optional; it does not renumber their identities or trigger new calculations. Full descriptions remain accessible when a large legend needs more space. The existing `pair=(i,j)` requirement is removed from the normal one-reference/many-formulation workflow. Optional terminal-pair and band selections limit the requested display using the existing observation/selection grammar.

The report's `illustrate` method calls this same recipe over its retained publication only when illustration was explicitly requested. Its default returns `nothing`, without loading a plotting backend. It does not own another plotting implementation. Saved and in-memory results must use the same selection descriptions and ordering.

## 8. Gauntlet becomes the caller of these operations

Keep the existing case and BenchmarkCalculation concepts. Add a constructor method of the existing `benchmark_definition` that accepts the loaded case, an explicit reference, and `formulations=...`, plus the shared report request. It owns materialization of the declared calculations, IDs and source-file binding; it does not create another benchmark type or computation route.

The reference may be an owned analytical formulation, PSCAD or another external backend with its own selections, or FEM. An explicit BenchmarkCalculation remains available when the reference uses an explicitly transformed problem or different execution options. Gauntlet performs no implicit homogenization or formula substitution.

Its visible run sequence is:

```text
validate declaration and report request
→ compute or recover the reference
→ compute or recover the candidate scalar/Gridspace using compute
→ report the completed operands through ReportBuilder
→ retain the report's unformatted publication and calculation evidence
```

One expensive reference calculation is shared by all candidate formulations. Recovery preserves completed calculations. Reanalysis changes only the requested comparison products; it does not run either numerical backend again.

Remove Gauntlet's formulation-comparison traversal, detailed tabulation, label construction and pair-only plot orchestration once their shared-owner methods replace them. Keep file integrity and external-format decoding at the Gauntlet storage boundary. Do not leave old/new active routes connected through forwarding aliases.

The Gauntlet documentation page consumes ReportBuilder's retained summary view according to the publication contract below. It neither computes nor defines its own grouping mathematics.

## 9. Default publication and explicitly requested detail

The standard Gauntlet Literate page is a summary of retained benchmark results. Numerical retention, page publication and interactive display have separate defaults:

| Output | Default | Explicit request |
|---|---|---|
| Retained calculation/comparison data | All completed reference/candidate matrices, formulation axes and records, and all requested per-term comparisons in the five bands | Reanalysis is needed only for comparison products not already retained. |
| Standard Literate page | Completion/availability counts, reference and formulation keys, band definitions and compact summary tables | Selected static illustrations can be added explicitly. |
| Ordinary report display | Compact summary and available detailed views | Full term tables and complete configuration records for a selected case, band or formulation. |
| Interactive matrix plots | None | Select a case/problem and quantities; the recipe overlays its candidate formulations and reference. Optional selections restrict terms, bands and formulations. |
| Static matrix plots/SVG exports | None | Explicit case/problem and quantity selections, followed by the existing export operation. |

The page has five ordered band sections. Each summary table has one row per case and complete candidate formulation, with an explicit reference key and compact quantity columns for Z, Y, R, L, G and C. Each quantity cell identifies the largest available per-term relative RMS discrepancy and its terminal pair. It also reports how many terms have unavailable relative RMS. When no relative value is available, show that fact and the largest absolute RMS with units; do not display zero percent. Absolute values and reasons for every term remain available in the detailed table. These are labelled maxima over matrix terms, never averages or rankings across formulations.

Reference/formulation keys link to a compact text description on the page. Common choices appear once, varying choices identify each complete formulation, and air/earth/mixed combinations remain explicit. Requested and actual band bounds and sample counts are shown once per distinct selection, rather than repeated inside every metric. The page includes instructions and retained-result links for reopening the detailed REPL workflow. It does not dump full configuration records or raw matrix entries.

All cases remain represented in the summary. Long text tables may be placed in collapsible sections; this does not authorize generating hidden plots. The default page contains no matrix dashboards, individual-term curves, heatmaps, plot thumbnails or embedded interactive plot payloads. It does not open a display, import Makie for benchmark illustrations, or export images.

An explicit static-publication request selects the case/problem, quantities and optionally formulations/terms/bands. The default list of requested illustrations is empty. Selected images are generated by the shared PlotBuilder recipe in that explicit operation, retained as artifacts and embedded by the page from those files. A documentation rebuild does not generate them again. Requesting one case's plot never implies plotting the other cases. An explicit batch request may select several cases; it must enumerate the intended scope.

`run_benchmark` retains the requested comparisons and summary without opening figures. `report` has `illustration=nothing` by default. `plot(...)`, an explicit report illustration request, or an explicit export request is what brings plotting into the workflow. This is the existing optional illustration stage, not a new publication controller or parallel reporting system.

Generating the Gauntlet Literate page reads retained summaries and any explicitly retained illustration files. It performs no solver calls, RMS calculations or plotting calls. Missing retained products are reported as unavailable; the page does not repair them by launching computation. This restriction concerns the benchmark page, not unrelated documentation examples that explicitly demonstrate plotting.

## 10. Drafts, acceptance, packaging and publication

Use one explicit sequence:

```text
run / resume in staging
→ inspect results and comparisons in the REPL
→ lock the accepted benchmark snapshot into the vault
→ package an explicit collection of locked snapshots
→ upload the package and bind its verified download
```

`lock` is the benchmark acceptance action. Keep that existing name instead of adding a synonymous `commit` command. A Git commit records source history independently. Acceptance means the user chooses to preserve and publish this calculation/comparison record; it does not assert agreement with PSCAD, FEM or a numerical error threshold.

The acceptance unit is one benchmark definition: case/problem, explicit reference, complete candidate formulation axis and report request, together with their completed results. Case ID alone is insufficient: the same physical case can have several benchmarks with different references, grids or formulation spaces. A campaign is a collection of these benchmarks. It must be possible to lock selected complete benchmarks without accepting or waiting for every other case in the campaign. Default locking of a whole campaign requires every declared member to be complete; filtering missing members silently is forbidden.

### Mutable drafts

`run --definition ... --directory ...` creates or replaces the drafts for the benchmark IDs supplied by that definition in the chosen staging directory. It never erases sibling benchmarks or the entire staging root. Running the same benchmark again is an explicit fresh execution request. `resume --directory ...` continues an interrupted attempt and reuses its completed calculations only after checking their inputs, implementation and output bytes. Reanalysis of retained operands remains a separate operation and does not call solvers.

Writes use temporary files/directories and a per-benchmark execution lock. Each attempt has its own identity. A completed new attempt becomes the current draft only after its result and comparison files are complete and verified. While an attempt is running or has failed, status exposes that fact; a retained previous complete result is labelled as previous and must not silently be presented or locked as the new attempt. Readers never combine old results with a new declaration. Successful replacement may discard the superseded transient attempt; staging is not an archive.

Drafts support the full shared REPL report/plot workflow. Inspection and optional exports do not accept a draft automatically. The record being inspected exposes its snapshot identity; locking can check an explicitly supplied inspected identity and must reject a mismatch. Concurrent writes and locking are mutually exclusive.

### Accepted snapshots in the vault

Extend the existing `lock_campaign` operation and CLI `lock` selection to copy the accepted benchmark set into a new immutable bundle. Keep the existing bundle reader and relative file bindings. A vault is a directory of these bundles, not a new Julia type, database or service.

Locking verifies completeness and the declared dependency inventory, copies the retained bytes, verifies the copy, records the selected benchmark IDs and their source snapshot identities, and makes the final bundle visible atomically. It does not calculate RMS, run a solver or generate figures. Source staging remains available for the next draft.

The bundle retains:

- case inputs and assets, frequencies, terminal coordinates, full requested/applied formulation selections and numerical controls;
- the reference and every candidate result, including the complete Gridspace axes;
- the report request, unformatted per-term comparison data and summary used during inspection;
- evaluated implementation/source records, dependency environment and backend version/configuration records, including captured dirty source bytes where applicable;
- raw files declared by each backend as retained evidence, plus explicitly selected exported illustrations and their display selections;
- an inventory of relative file names and checksums, format versions, snapshot identities, acceptance time and an optional user note.

No required dependency may point only into mutable staging or an original remote work directory. Transient scratch files, stale failure markers and incidental files are not swept into the bundle merely because they share a directory. The reader checks the full declared inventory. A changed, missing or unexpected payload fails verification. Moves between directories must not change scientific identity or require loading a solver to read the results.

Every Gauntlet writer rejects a vault bundle as an output target. There is no force-overwrite or unlock operation for accepted snapshots. Checksums detect later filesystem changes; the filesystem is not claimed to provide tamper-proof storage. A correction creates a new draft and a new locked snapshot, optionally recording which snapshot it replaces. New bands, revised comparisons or new publishable illustrations also produce a new accepted snapshot; viewing or exporting to a separate user directory remains read-only with respect to the vault.

### Explicit packaging command

Expose `package_collection` through `lcm gauntlet package`. Its input is an explicit list of locked bundle paths/identities, or a TOML release definition containing that list, collection name, version and description. It must not scan staging, select the newest files implicitly, accept drafts on the user's behalf or include every vault entry by accident.

Recommended proposed CLI syntax:

```sh
cli/lcm gauntlet run --definition benchmarks.jl --directory gauntlet/.artifacts/staging/study
cli/lcm gauntlet lock --directory gauntlet/.artifacts/staging/study \
    --benchmark soil_comparison --output gauntlet/.artifacts/vault/soil-r1
cli/lcm gauntlet package --definition release.toml --output gauntlet/.artifacts/releases/soil/v1.0.0
```

The new `--benchmark` selection and `package` command are proposed additions. A release definition uses an ordinary TOML document, for example:

```toml
collection = "soil"
version = "1.0.0"
description = "Accepted formulation comparisons"
bundles = ["gauntlet/.artifacts/vault/soil-r1"]
```

Relative bundle paths resolve against the release definition file. Packaging records their verified identities in the release inventory. Duplicate or conflicting benchmark snapshots must be resolved explicitly in that definition. A collection can combine accepted benchmarks from different campaigns and execution revisions; each keeps its actual execution record.

Packaging verifies locked inputs, assembles their complete retained files and release inventory, and uses the existing `Pkg.Artifacts` operations to produce the archive, archive SHA-256, artifact tree hash and package metadata. It reopens the packaged result to verify that comparisons remain readable independently of staging. Failed packaging leaves no apparently complete release. A retry with the same release definition verifies and reuses the existing package; different contents require a new version. Remove the current force-replacement path for accepted releases.

Package identity is determined by the accepted inputs and explicit release definition. The checkout used for packaging must not replace the recorded numerical implementation with its current HEAD, or require re-executing a calculation merely because that checkout has changed. If a packaging tool revision is recorded, it is separate from each benchmark's execution revision. No solver, RMS calculation, plot, Git tag, upload or artifact-binding update occurs during packaging.

### Publication and the standard page

Publication consists of uploading the exact packaged archive to the chosen host, then registering its verified download through the existing `bind_published_artifact` operation exposed as `lcm gauntlet bind`. Keep the name `bind` honest: it records an already uploaded artifact. Use the hosting service's existing upload command; a generic multi-host upload framework is outside this work. The host/URL is an explicit publication input, not an inference from the Git remote.

Binding verifies the served archive checksum and extracted artifact tree against the package, then records a version-specific download binding. A failed upload or verification leaves the accepted bundle and package intact and does not advance the public binding. An existing collection/version cannot be rebound to different bytes. A convenience current-version binding may be updated explicitly, while old version bindings remain recoverable. Retrying the same publication is harmless. Publishing never silently rebuilds the archive or creates plots.

The standard published Literate page uses explicitly pinned published collection versions/hashes. It never reads whichever staging directory was modified last. Before publication, an explicit local page preview may read a selected draft or locked bundle and labels it accordingly; it remains summary-only under section 9. The page shows the collection version and accepted snapshot identities so every displayed comparison can be reopened from its artifact.

No new lifecycle hierarchy or duplicated state controller is required. Existing run/resume, bundle manifests/readers, `lock_campaign`, `package_collection` and `bind_published_artifact` own these actions. CLI parsing supplies their explicit inputs. Move the script-only entrypoints into this CLI and remove obsolete parallel scripts/callers after updating their real consumers. Preserve only the enumerated historical file readers, not a second path that packages mutable drafts.

## 11. Implementation order and acceptance

1. Lock the output-coordinate and selection-record contracts, including scalar references and recovered axes. Keep the bounded numerical result types intact.
2. Extend `compare` over retained formulation axes. Verify scalar equivalence and explicit pairing behavior.
3. Extend BenchmarkTableDefinition and ReportArtifact; implement shared selection, complete tables and summaries. Verify one analysis execution and retention before display formatting.
4. Add the formulation-aware ParametricResult/publication recipe to PlotBuilder. Verify full overlays and shared labels.
5. Route Gauntlet declarations, execution, recovery and page summaries through those operations, and delete the replaced owners/callers.
6. Implement replaceable per-benchmark staging and selected-benchmark locking; expose package/bind in the CLI and make packaging consume only explicit accepted bundles. Update publication selection in the Literate page and remove obsolete packaging entrypoints.
7. Run the complete relevant quality, numerical, reporting, plotting/export, persistence, lifecycle and documentation checks. Record generic/type/module counts separately from added methods before committing.

Required behavioral checks:

- A real formulation Gridspace gives the same matrices and ordering as its individually selected computations; product/zip cardinality is preserved.
- One reference and N candidates produce N comparisons per quantity/band, with no repeated reference solve.
- Scalar, triple and complete generic formula selections retain exact identities. Changing only Y, an EHEM rule, a constitutive choice, an integration option or a hook remains distinguishable.
- All candidate/reference curves are present in every requested matrix cell; labels/styles match the retained selection, survive filtering and survive save/reload.
- Multiple problem points remain separate; undeclared many-reference pairing is rejected.
- Five-band groups appear in the declared order; the display records actual samples and does not interpolate.
- Lossless/numerically zero G, empty bands, asymmetric matrices, explicit unit changes and incompatible coordinates are exercised.
- Plot/table construction works without Gauntlet loaded. Loading retained comparisons performs no solve or RMS calculation. Display changes do not alter stored values.
- Generating the standard Gauntlet Literate page over the full case catalogue performs zero solver, RMS and plotting calls, and writes no figure files. Every retained case/formulation appears in its summary with the recorded availability status.
- Default benchmark/report calls produce no illustration and open no display. An explicit request for one case plots only that case, includes its selected formulations/reference, and does not recompute comparisons. Static publication embeds only explicitly retained figures.
- Missing relative RMS remains unavailable in the summary, including numerical-zero G. Summary maxima retain their quantity, normalization, units and terminal coordinates; detailed absolute values and reasons remain recoverable.
- Invalid requests fail before solver execution. Gridspace compilation still accepts the same consistent numerical result types.
- A third compatible result implementation can use the public observation/coordinate contract without changing Gauntlet's orchestration.
- Re-running one draft replaces only its declared benchmark; sibling drafts and accepted bundles remain byte-for-byte unchanged. Interrupted replacement retains a clearly identified previous result and never mixes generations. Resume verifies and reuses only matching completed calculations.
- Locking selected complete benchmarks works while unrelated campaign members are incomplete. Locking the declared whole campaign rejects an incomplete member; changing a draft during acceptance or supplying a stale inspected identity is rejected.
- A locked snapshot remains readable and plottable after moving it and deleting its original staging directory. All writers reject vault targets; missing, changed and unexpected payload files are detected.
- Packaging rejects drafts, verifies accepted dependencies and produces a standalone readable archive. Same-input retries preserve the release; changing an existing version's contents is rejected. Failures do not leave a complete-looking bundle or package.
- Publication verifies downloaded bytes/tree before binding; upload or verification failure does not change the public selection. Published version bindings remain fixed, and rebuilding the page uses the pinned versions with zero solves, RMS or plotting calls.

No new generic name, result wrapper, abstract hierarchy, module, registry, alias map or orchestration context is proposed. The changes add methods to existing operations and fields to existing owned definitions/artifacts, with the specific warrants above. UQ moment workflows remain supported under their declared current comparison capabilities; this work does not silently reinterpret moment arrays as line-parameter curves or add unimplemented moment-band behavior.

The persisted-format review found no retained local result files at the checkpoint. The existing scalar, Gridspace and UQ fixtures define the supported historical readers. The implementation writes calculation/comparison/bundle schema 2 and campaign/release schema 3, retaining the exercised earlier readers.


## 12. Completed implementation and validation

Implemented in `release/v0.2.0-pscad-ext` after checkpoint `0d464f95`.
The public examples and exact CLI arguments are in `gauntlet/README.md`.

- `compute(problem, formulation_grid)` retains its existing combinatorial route.
  Shared `compare`, `report` and `plot` methods consume the resulting formulation
  axes. Scalar references, explicit space-to-space pairing, complete selection
  records and actual reduced output coordinates are retained.
- ReportBuilder owns the five bands, complete per-term tables and labelled maxima.
  Saved analyses are selected and displayed without another RMS calculation.
  Multiple saved analyses retain distinct snapshot identities. An independent
  result implementation was exercised through public observation/coordinate methods.
- PlotBuilder overlays every selected formulation and reference in every matrix
  cell. Both off-diagonals, equal curves, filtering, multiple problems, retained
  bands and stable styles are covered. Full selection records remain available
  alongside the compact legend. SVG export preserves the selected viewport with
  and without interactive controls.
- Repeated draft runs replace only the selected benchmark. Failed replacement
  keeps the previous completed attempt explicitly accessible. Resume checks
  declarations, source bytes and the dependency environment before reusing work.
  Legacy cleanup and staging replacement reject recursive deletion of accepted
  bundles, including bundles nested beneath the requested directory.
- `lock` accepts selected complete benchmarks into self-contained checked bundles;
  `package` accepts an explicit list of those bundles; `bind` verifies downloaded
  archive bytes, the extracted tree and the release definition before updating
  an immutable version binding. Obsolete packaging scripts were removed.
- The standard Literate page reads pinned artifacts and produces compact summaries.
  Matrix figures require an explicit request. Only explicitly retained illustrations
  are embedded. `docs/gauntlet.toml` currently selects no published collection.

Validation covered:

| Checks | Result and scope |
|---|---|
| Core regression suite | 11,416 passed; numerical engine, formula dispatch, result transport, reports and displays. |
| Gauntlet | Full 856-check run passed; expanded lifecycle, record, UQ and REPL follow-ups passed after the final changes. The final affected saved-result group passed 200 checks; the final archive/vault run passed 50, including recursive cleanup protection. |
| Quality | 1,102 passed; Aqua, explicit imports, semantic economy, ownership and transport constraints. |
| Full visual suite | 880 passed; existing plots, UI and export behavior. |
| Final formulation/report checks | 234 passed; matrix overlays, retained labels, five-band reports and independent result protocol. |
| FEM numerical suite | 280 passed using GetDP, including reduced matrices and constitutive behavior. |
| Report/display suite | 1,618 passed, including workbook and text output. |
| Documentation | Doctests and the complete documentation build passed; existing HTML/search-index size warnings remain. |
| CLI lifecycle | Actual `run`, `status`, `lock`, `package` and `bind` invocations passed with temporary data and a local archive URL. The accepted results and packaged summary reopened after deleting staging. |

The CLI exercise used a local `file://` publication URL. The full remote PSCAD
catalogue was not rerun and no release was uploaded or added to the repository's
published artifact bindings. Acceptance of scientific benchmark results remains
an explicit user action through `lock`.

A runtime inventory of the package and Gauntlet compared with the checkpoint gives:

| Owned definitions | Checkpoint | Completed worktree |
|---|---:|---:|
| Generic bindings | 879 | 876 |
| Type bindings | 178 | 178 |
| Modules | 32 | 32 |
| Authored methods | 3,239 | 3,269 |

The removed generics are `Gauntlet.benchmark_comparisons`,
`Gauntlet._selection_record` and `Gauntlet._package_collision`. Optional Makie
extension methods are outside this runtime count; its new recipe extends the
existing `plot` operation. The 60 inventoried formula, earth/material and transform
implementation files are byte-identical to the checkpoint. No equation ingestion
or formula aliasing is part of this change.
