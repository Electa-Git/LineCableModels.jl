# PlotBuilder refoundation execution

Approved scope: the complete PlotBuilder refoundation specification supplied in
this conversation, including the A01–A19, B01–B18, and C01–C26 preservation
inventory. This is an execution record, not a replacement design.

## Baseline

- Branch: `release/v0.2.0`.
- HEAD: `4f72bf45aeb4558e0fd48cd5f45e7a84d86f478b`.
- Worktree clean before execution. One agent; no commits, branch changes,
  dependency changes, or external solver campaigns.
- Existing planning evidence at the same HEAD: 430 assertions for retained
  interpretation and native axis formatting, and 91 for retained boundaries.
  These are baseline evidence only, not implementation validation.
- Verified environment: Julia 1.12.7; root project with the existing visual
  environment on `JULIA_LOAD_PATH`. Dependency manifests remain unchanged.

## Decisions retained

- Effective result presentation dispatches publicly on `ObservedResult` or an
  ordinary vector. Raw conveniences acquire observations and delegate.
- `layout` alone defines nominal capacity. Quantities remain separate families;
  matrices preserve topology, explicit diagonals and preview collections flow.
- Observation owners retain scientific interpretation, grouping, units, and
  band selections. PlotBuilder owns native presentation and interaction.
- One live `UIPlot`, one shell state, native controls and guides. Window closure
  preserves a retained figure. Later content changes refit only that window.
- No historical modes, compatibility aliases, or new scientific plot records.

## Stage evidence

Stages 0–8 are complete. The entries below preserve the development sequence,
including failed checks and their corrections. Final gate results and coverage
limits are recorded at the end; chronological pending statements describe their
respective checkpoints, not unfinished implementation.

## Feature and deletion evidence

All A01–A19, B01–B18, C01–C26 remain required. The complete owner/evidence and
removed-to-surviving maps appear below the chronological stage record.

Stage 1 implementation so far:

- Typed observed `plot` now owns request/product/group/page orchestration.
  Deleted `_addon_line_pages` and its route-dependent defaults.
- `_addon_capacity`, `_semantic_line_pages`, `_addon_flow_pages` implement the
  sole capacity rule. Deleted `_semantic_line_layout_mode`, paired-quantity
  behavior, and `blocks` handling in result presentation.
- `Grammar.observation_product(...; band, reference_id)` interprets completed
  records. Removed plotting-side first-candidate band selection.
- Selected products provide extents and domains. Native scatter consumes
  assembly/index products; global short-frequency rejection is gone.
- Candidate slots precede filtering; the separate reference has its own role.
  A first-run failure caught its zero slot reaching the glyph sampler; corrected
  by giving reference glyphs their independent deterministic phase.
- Relevant preservation: A19, B10–B14, C21, C23, C24, C25, C26. Geometry,
  lifecycle, and the complete public surface are not yet claimed complete.

## Commands and outcomes

- `git status --short`: clean before edits.
- `git branch --show-current`, `git rev-parse HEAD`: baseline above.
- `julia --version`: 1.12.7.
- With `JULIA_DEPOT_PATH=/tmp/observed-plan-depot:/home/amartins/.julia`,
  `JULIA_LOAD_PATH='@:/home/amartins/Documents/KUL/LineCableModels/test/visual:@stdlib'`,
  and `LINECABLEMODELS_TEST_PLOTTING=true`:
  `julia --project=. --startup-file=no test/runtests.jl integration/retained_interpretation integration/observed_grouping`
  — first run: 59 passing assertions, one error (reference glyph index zero),
  `/tmp/plotbuilder-stage1.log`. Expanded corrected run pending in
  `/tmp/plotbuilder-stage1b.log`: 92/92 assertions passed (4 test items).
- Same environment, selectors `'retained band associations' integration/report_result_protocol unit/engine/retained_boundaries`:
  106/106 assertions passed, 8 items; `/tmp/plotbuilder-stage1-boundaries.log`.
  Includes independent source ownership, conflicting reference/band records,
  sample intersections, and all 91 retained-boundary checks.
- Same environment, selector `integration/native_makie_axis_format`: first run
  passed 377 assertions and failed 8 assertions of superseded adaptive-log
  behavior. Replaced those assertions with signed-log behavior and strict native
  `log10` page-wide preflight/rollback checks. Removed `signed_ylog` from callers.
- Same environment, selectors `'callable controls own' integration/native_makie_axis_format integration/native_makie_series_styles integration/native_makie_uncertainty_visibility`:
  `/tmp/plotbuilder-stage2-controls.log` passed all 391 axis assertions and all
  30 callable-widget assertions. Remaining failures exposed old caller defaults:
  `blocks`, implicit all-interval sampling, implicit ellipsis, and expecting no
  automatic candidate markers. These callers now state the required display
  choices explicitly. Their rerun and new native-data-update coverage are pending.

Stages 2–3 surviving owners:

- `_addon_set_axis!` resolves adaptive X/Y scales and strict native functions;
  public `axisscale!`, constructors, and controls share this operation. There is
  no `entry.signed`, `signed_ylog`, or recipe `scale_controls` switch.
- UIPlot holds count-independent native axis/colorbar vectors and the shell's
  actual status observable. It is created before controls.
- `controls.jl` owns native widget slots and removable subscriptions. Standard
  and custom widgets use `_addon_widget!`; the old toolbar builder and provisional
  handle reference were removed. Public reset/scale actions work without controls.
- Native content lifetime owns callbacks; no window-close disposal callback is
  installed. Actual GL close/redisplay coverage remains pending.
- Geometry refitting after widget changes is a remaining stage 5 dependency.
- The series rerun `/tmp/plotbuilder-stage2-series.log` passed 924 assertions;
  its four failures expected no Scatter object when all uncertainty intervals
  are drawn. Empty native glyph handles retain legend identity; the migrated
  assertion now verifies that all such handles contain no drawn points.
  All 660 uncertainty-coordinate/visibility assertions and all seven new
  native-data-update assertions passed.
- The final combined stages 2–3 run completed successfully in
  `/tmp/plotbuilder-stage2-final.log`, using selectors
  `integration/native_makie_axis_format integration/native_makie_series_styles integration/native_makie_uncertainty_visibility 'callable controls own'`.
  It also exercises the latest per-binding full-uncertainty bounds cache and
  independent native line/bar update handling.
  All 26 selected test items passed; exact totals are in the log. This closes
  the focused numeric/series/widget gate before guide refoundation. Backend
  closure and widget geometry/export integration remain later-stage gates.

## Validation status at the stages 2–3 checkpoint

At this checkpoint, the complete visual suite, affected observation/report and
extension checks, GL/WGL coverage, and documentation builds were still pending.
Their subsequent outcomes are recorded below.

Stage 4 implementation and focused evidence:

- `guides.jl` owns the current placement of complete native guide objects.
  `_addon_compose_guides!` is shared by initial construction, guide mutators,
  and viewport changes. Factories receive assigned slots and allocate no docks.
- Removed `_addon_figure_legend_target!`, the separate initial shared-dock path,
  the old live figure/panel legend implementations, and anchor parsing.
- `figurecolorbars!` removes/restores a count-independent native collection.
  Companion labels bind to the actual bar label/font/size/color attributes.
- Legend and scale removal retains native objects and independent visibility
  edits. Guide-owned subscriptions follow native block lifetime.
- Removed the now-unused `_addon_equal_matrix_cells!` and
  `_addon_refit_matrix_block!` budget/sibling machinery. Replacement frame
  calibration remains a stage 5 requirement, not a completed claim.
- Relevant preservation: B15–B17, C01–C03, C08, C18. Frame fitting and complete
  creation-order geometry evidence are still being completed.
- Same visual environment, command selectors `integration/native_makie_legends
  'shared series attributes' 'uncertainty legend shading follows'`:
  112/112 assertions passed in `/tmp/plotbuilder-stage4-legends.log`.
- Stages 2–3 final combined run: 1349/1349 assertions, 26 items, 249.96 seconds.

Stages 4–5 focused gate:

- `layout.jl` now owns nominal capacity, flow pagination, one-time native frame
  calibration, local content refitting, and panel-specific physical aspect.
- Removed `_addon_center_aspect_canvas!`, `_addon_landscape_size`,
  `_native_preview_layout`, and the GLFW aspect-ratio lock. Preview collections
  now paginate rather than reject a smaller capacity. Their original indices
  remain panel identities, and their material ranges are acquired once.
- Guide/title/widget mutators preserve current frames and views while fitting
  only their own window. Automatic flow reflow is local to existing page members;
  explicit matrix and flow layouts do not install it.
- Retained native guide objects survive removal, moves, and restoration. Explicit
  native legend bank/orientation edits release automatic wrapping. The figure
  frame region includes declared matrix selection holes.
- Same visual environment, selector `'nominal capacity calibrates'`:
  25/25 assertions passed, `/tmp/plotbuilder-stage5-capacity.log`.
- Same visual environment, selectors `integration/native_makie_legends
  integration/native_makie_axis_format 'compact preview and material scheme geometry'
  'nominal capacity calibrates' 'callable controls own'`:
  582/582 assertions, 20 items, 214.28 seconds;
  `/tmp/plotbuilder-stage5-combined.log`.
- Additional creation-order and hidden-scale restoration checks were added
  after that combined run and are not yet counted as passing evidence.
- Stage 6 implementation now separates shell presentation in
  `export_presentation.jl` from path/backend/I/O handling in `native_export.jl`.
  Removed matrix export exceptions and hard-coded chrome-row tests. Its new
  success/failure, custom-control, and exact-restoration evidence is pending.

## Integrated validation and closure work

- First complete `tag:visual` run: 2,821 passing assertions, 15 failures,
  7 errors, 61 maintained items, 862.77 seconds; log
  `/tmp/plotbuilder-visual-full.log`. The failures identified complete guide
  measurement, obsolete paired-quantity assertions, reference endpoint/glyph
  collision, and adaptive empty-axis handling. This run is not a green gate.
- The following focused run passed 2,525 assertions and failed three preview
  legend assertions; `/tmp/plotbuilder-native-convergence.log`. All axis,
  uncertainty, matrix, UQ, widget, and export sections passed. Subsequent guide
  changes require the final rerun recorded below.
- `docs/doctest.jl`: 1/1 doctest passed, 35.2 seconds;
  `/tmp/plotbuilder-docstest.log`.
- Initial `docs/make.jl` reached documentation checking and reported the five
  new UI operations missing from the canonical reference page. Added those
  entries; the full build is being repeated. No numerical campaign was run by
  the case catalogue generation.
- The first GL attempt could not resolve GLMakie from the visual-only stack.
  Adding the existing global environment resolved it, but sandboxed X11 access
  failed. The installed `dev/plotting` environment supplies the matching
  Makie 0.24.13 backend stack without changing dependency files.
- GL native check outside the sandbox: 42/42 assertions passed, including
  unloaded-Cairo refusal, late-loaded SVG export, native zoom/pan restoration,
  active-backend preservation, window closure, subsequent editing/export,
  redisplay, and exactly-once custom callbacks. Log:
  `/tmp/plotbuilder-gl-native.log`.
- WGL sandboxed import failed during Reseau's local socket precompilation.
  Outside the sandbox, the installed environment built all 16 gallery panels;
  `/tmp/plotbuilder-wgl-unrestricted.log`. This verifies construction through
  WGL; it does not claim browser-driven interaction coverage.

Backend commands use the same depot and add
`/home/amartins/Documents/KUL/LineCableModels/dev/plotting` after the visual
project in `JULIA_LOAD_PATH`:

```bash
LINECABLEMODELS_GL_ARTIFACTS=/tmp/plotbuilder-gl-artifacts \
  julia --project=. --startup-file=no dev/plotting/manual_gl.jl
LINECABLEMODELS_WGL_MANUAL=true LINECABLEMODELS_WGL_SMOKE=true \
LINECABLEMODELS_WGL_ARTIFACTS=/tmp/plotbuilder-wgl-artifacts \
  julia --project=. --startup-file=no dev/plotting/manual_wgl.jl
```

## Surviving calls and removed authorities

| Surviving call or owner | Deleted competing path | Preserved operation |
| --- | --- | --- |
| Public `plot(::ObservedResult)` and `plot(::AbstractVector{<:ObservedResult})` | Broad result-page orchestration wrapper | Raw, retained, named, parametric, UQ, and report plotting |
| Raw owned-type conveniences → observation construction → public `plot` | Route-dependent scientific preparation | Partial requests complete through Grammar; retained settings stay retained |
| `Grammar.observation_product(points, request; band, reference_id, unit, frequency_unit)` | Plot-side first-candidate band indexing | Each candidate/reference keeps its own saved sample identities |
| Selected product coordinates and `_semantic_line_pages` | First-unrelated-product extent/domain lookup; global short-sample return | Full modal/rectangular matrices, explicit diagonals, single samples |
| `_addon_capacity` and `_addon_flow_pages` | `blocks`, pairing modes, preview undersized-layout rejection | One capacity; quantity separation; compact residual pages |
| `_addon_set_axis!`, `_addon_reset!`, managed numeric formatting | `entry.signed`, `signed_ylog`, reference/recipe scale permission | Adaptive signed X/Y, strict native scales, safe page actions |
| `_addon_comparison_styles` with original candidate slots | Reference-shifted style indices and reference-dependent defaults | Stable candidates, separate reference identity |
| `_addon_finish!`, `_addon_widget!`, public widget/axis operations | Provisional handle and separate toolbar/callback paths | One live UIPlot, actual status, native subscription lifetime |
| `_addon_compose_guides!` and native factories in assigned slots | Initial shared dock, separate live dock, anchor parsers | Native guides, all placements, removal/restoration and wrapping |
| `_addon_calibrate_frames!`, `_addon_fit_frames!`, local presentation edits | Monotonic matrix budgets and sibling fitting graph | Equal initial frames, actual residual extents, local later edits |
| `_addon_fit_panel_aspects!` | Whole-canvas square fitting, landscape expansion, GLFW ratio lock | Physical preview scaling and freely resizable windows |
| `_addon_export_presentation!` plus SVG writer | Export-specific root-row and matrix blank-strip rules | Same-figure SVG, complete restoration, explicit Cairo boundary |

There is no compatibility alias for the removed capacity or anchor options.
Their explicit rejection diagnostics are not retained implementations.

## Feature-to-owner and maintained evidence map

This map identifies the maintained behavioral checks. Final executed aggregate
results are recorded separately; listing a test here does not turn an earlier
failed run into a pass. File names below are under `test/integration` unless
specified otherwise.

| Feature | Surviving owner | Behavioral evidence |
| --- | --- | --- |
| A01 | `controls.jl`, scale/reset operations | `native_makie_axis_format`, callable controls with/without toolbar |
| A02 | `_addon_set_axis!` | Signed X/Y and route-independent adaptive scale assertions |
| A03 | `_addon_scale` native reversible transform | Signed transform round trips and near-zero/large-value axis tests |
| A04 | Page preflight and rollback in `_addon_set_axis!` | Invalid strict-log/custom-scale multi-axis rollback |
| A05 | Scale action and reset view handling | Orthogonal views, full/partial configured limits |
| A06 | `_addon_axis_format!` | Pixel/font/rotation density and large-offset tests |
| A07 | Managed mantissas and `_addon_axis_label` | Live range-dependent engineering multipliers |
| A08 | `_addon_linear_tickformat` | Unique large-offset labels, negative zero, trailing zeros |
| A09 | `_addon_decade_ticks` | Short positive logarithmic spans and screen positions |
| A10 | `_addon_decade_ticks` | Two-decade transition and broad decade ticks |
| A11 | Native tick/formatter ownership | Labelled/function/numeric ticks and return to automatic |
| A12 | `_addon_constant_limits`, reset corrections | Positive/negative near-constant padding and log padding |
| A13 | Local admission and empty-axis scale resolution | Zero/unavailable products in `observable_resolution_plots` |
| A14 | Retained arrays passed unchanged to primitives | Tiny unclipped signed-signal tests |
| A15 | Native target-limit lifecycle | Native zoom/pan, limits, partial bounds, SVG zoom preservation |
| A16 | Full interval support in axis bindings | Unsampled extreme intervals and log eligibility |
| A17 | Native visible/autolimit bounds | Independently hidden intervals and native overlays |
| A18 | Native converter/scale registration | Categorical, datetime, custom native scales |
| A19 | Local trace admission | `retained_interpretation`: single frequency and local empty subsets |
| B01 | Native line/scatter owners and dependent maps | `native_makie_uncertainty_visibility` |
| B02 | Writable legend owner bindings | Native hide/show/solo/mixed-state actions |
| B03 | Visibility/color subscriptions | Independent bar edits survive guide/layout notifications |
| B04 | Native primitives in one axis transform | Projected interval endpoints under log and signed log |
| B05 | `_addon_series_styles!` | Coincident intervals, stable nested widths, explicit overrides |
| B06 | `_addon_glyph_indices` | Sparse original-sample markers across curves |
| B07 | Separate interval/marker slots | Short curves and reference endpoint tests |
| B08 | Sampling mode and native marker overrides | All-interval mode and explicit all-sample markers |
| B09 | Viewport-owned glyph sampling | Resize density changes preserve curves/support/views |
| B10 | Separate reference role in shared styles | Black/hollow defaults, endpoint markers, explicit overrides |
| B11 | Invocation defaults and shared style resolver | Candidate styling with/without reference |
| B12 | Existing perceptual palette | Stable prefix and exhaustion tests |
| B13 | Original candidate slots before filtering | Matrix/UQ/formulation overlays across quantities/pages/filters |
| B14 | Grammar descriptions and groups | `observed_grouping`, `formulation_overlays` |
| B15 | Native guide factories and common composition | Figure/panel/inside legends, rich labels, banks/orientation |
| B16 | Measured wrapping and explicit ellipsis | `native_makie_legends`, responsiveness and overflow interaction |
| B17 | Retained native guide objects | Repeated remove/move/restore preserving visibility/styles/views |
| B18 | Native data subscriptions/full-support distinction | Explicit line and interval updates without source acquisition |
| C01 | `_addon_compose_guides!` | All sides, positive grid slots, titles, native alignment/margins |
| C02 | Complete colorbar groups and endpoint extents | Endpoint text containment, companion label styles |
| C03 | Shared guide composer | Creation-order, move/remove/restore, abandoned-track checks |
| C04 | Native figure dimensions and GL adapter | Portrait sizes, manual resize, GL screen checks |
| C05 | Fixed topology versus local automatic flow | Matrix resize and automatic preview reflow identities |
| C06 | One-time native frame calibration | 3×3 residual pages and measured equal data frames |
| C07 | Panel-owned physical aspect | Circular previews, wide systems, 1×4 collections |
| C08 | Live `UIPlot` fields/native objects | Post-construction native edits and restored guide handles |
| C09 | Standard widgets through shared allocation | Icons, status, no-controls geometry, native actions |
| C10 | Callable widgets and owned slots | Once-only actions, failed builders, composite cleanup/export |
| C11 | Once-only native canvas axis discovery | Nested canvas topology, native plots, scoped panel guides |
| C12 | Optional extensions/backend dispatch | `test/extensions/makie.jl`, backend-first scripts, no-Cairo checks |
| C13 | Existing figure saved through Cairo | Native overlays, uncertainty, patterns, rich labels in SVG |
| C14 | Shell publication presentation | LaTeX font roles and restored native fonts/background |
| C15 | Presentation `finally` restoration | Success/failure exact views/size/styles; sibling isolation |
| C16 | Writer backend preflight | GL no-Cairo refusal followed by explicit late Cairo loading |
| C17 | SVG path/opening writer | Unique names, sanitization, overwrite refusal, directory/opener checks |
| C18 | Preview declaration recipes plus common shell | Grouping/pattern/title/offset/material preservation tests |
| C19 | Viewport-owned soil/sky content | Finite horizontal spans, original depths, autolimit exclusion |
| C20 | Retained statistical native drawing | `native_monte_carlo`, UQ retained histogram/CDF/Q–Q checks |
| C21 | Selected full-matrix coordinates | Modal zero/residual nine-panel fixtures and rectangular matrices |
| C22 | Native scene/block subscriptions | Composite deletion, removal, repeat actions, actual GL redisplay |
| C23 | Retained assembly products/native categorical scatter | First-seen category union, real points only, uncertainty |
| C24 | Retained scalar/vector index drawing | Negative vectors, scalar points, no frequency inference |
| C25 | Single capacity/pagination owner | Matrix, explicit diagonal and preview page counts/identities |
| C26 | Grammar retained-band interpretation | Different sample subsets, conflicting/missing bands, separate reference |

Further closure evidence:

- Observation/report/core-extension command with selectors
  `integration/report_result_protocol integration/retained_reporting
  unit/engine/retained_boundaries extensions/makie` passed **132/132** assertions,
  10 maintained items, 297.23 seconds; `/tmp/plotbuilder-observation-final.log`.
  This includes the independent test-only primary owner, source poisoning,
  explicit retained selection, empty-band diagnostics, small/zero RMS preservation,
  arbitrary precision, and core loading without any graphics extension.
- Repeated the native GL command with `--backends-first` and
  `LINECABLEMODELS_GL_ARTIFACTS=/tmp/plotbuilder-gl-backend-first`:
  **42/42** assertions passed; `/tmp/plotbuilder-gl-backend-first.log`.
- The guide minimum-size check passed all **27 endpoint-label assertions** in
  `/tmp/plotbuilder-guide-minimums.log`. It also exposed that standalone scales
  must continue filling the main content canvas rather than using a dock's
  finite natural length; that distinction now belongs to the common placement
  operation. Native length overrides remain native.
- Direct native guide/font changes enter the same local frame-preserving fit as
  public mutators. Watches are installed once on their native block scenes.
  Native canvas panel legends now use registered plot labels and the common
  outer guide composition without replacing caller-owned nested grids.
- Removed unused matrix minimum-cell constants, unused imports, and the dormant
  `export_mode` construction-theme branch. Publication typography has one
  surviving owner in temporary export presentation.

## Public calls and layout examples

The effective result methods are:

```julia
plot(observed::ObservedResult, selection=nothing; ydata=nothing, ...)
plot(observed::AbstractVector{<:ObservedResult}, selection=nothing; ydata=nothing, ...)
```

`selection` and `ydata` are alternative spellings, never simultaneous. All owned
raw conveniences acquire through Grammar and call these methods. Named tuples
retain their labels; report artifacts pass their observations and separate atomic
observed reference. See `src/plotbuilder/interfaces.jl` and the authoritative
`docs/literate/plotting.jl` for the full keyword surface and constructed fixtures.

```julia
using LineCableModels
using LineCableModels: plot
using CairoMakie
frequency = [1.0, 10.0, 100.0]
impedance = reshape(complex.([1.0, 2.0, 3.0], [2.0, 3.0, 4.0]), 1, 1, :)
raw = LineParameters(impedance, impedance .* 1e-6, frequency)
r = @observe R[:, :, :]
a = plot(raw; ydata=(r,), clip=true, length_unit=:kilo)
o = ObservedResult(raw, (r,); complete_pairs=true, clip=true, length_unit=:kilo)
b = plot(o; ydata=(r,))
pair = (@observe(R[:, :, :]), @observe(L[:, :, :]))
c = plot(raw; ydata=pair, layout=(1,1))
d = plot(ObservedResult(raw, pair); ydata=pair, layout=(1,1))
```

The last two calls each return separate R and L handles. The sole capacity rule:

| Selected population | `layout` | Result per quantity |
| --- | --- | --- |
| Full 3×3 matrix | `(1,1)` | Nine one-coefficient figures |
| Full 3×3 matrix | `(2,2)` | Four figures: 2×2, 2×1, 1×2, 1×1 |
| Three explicitly diagonal coefficients | `(1,2)` | Two panels, then one residual panel |
| Four cable previews | `(1,2)` | Two figures with two original preview identities each |

Gridpoints add traces. Quantities/statistical meanings remain separate figure
families. No `blocks` or quantity-pairing mode survives.

The added UI operations are declared/exported by PlotBuilder and the root:

```julia
figurecolorbars!(p::UIPlot; position, group_attributes, native_bar_attributes...)
axisscale!(p::UIPlot, dimension::Symbol, scale; panel=nothing)
resetview!(p::UIPlot; panel=nothing, x::Bool=true, y::Bool=true)
addwidget!(builder, p::UIPlot, key::Symbol; event=nothing, callback=nothing, success=nothing)
removewidget!(p::UIPlot, key::Symbol)
```

Guide arguments omitted at mutation preserve state; explicit `nothing` hides the
guide. Existing `figurelegend!`, `panellegend!`, `figuretitle!`, `paneltitle!`,
`export_svg`, previews, material primitives and retained statistical verbs remain.
`addwidget!` returns its native block/grid, `figurecolorbars!` the current native
collection, and scale/reset/removal return `p`.

Final geometry gate:

- Selectors `'local native edits and flow' 'guide creation order'
  integration/native_makie_legends 'nominal capacity calibrates'
  'colorbar endpoint labels'`: **197/197** assertions passed, six maintained
  items, 262.69 seconds; `/tmp/plotbuilder-geometry-final.log`.
- Covers full endpoint text, complete native scale groups, fractional placement,
  guide creation order, native font edits preserving frames, long custom controls,
  automatic wide/tall reflow with stable identities, nested native canvas panel
  guides, equal nominal frames, compact residuals, and sibling isolation.
- The guide composer now allocates and sizes both leading and trailing spacer
  tracks. Leaving the trailing span implicit had allowed native grid alignment
  to distribute space a second time; the measured fractional-placement test
  detects that defect.
- Existing fixed canvas dimensions are refreshed when decorations change.
  Temporary SVG presentation also snapshots/restores these native dimensions,
  in addition to the figure size and views.
- Automatic legend wrapping uses the existing presentation guard while updating
  native banks. A resize-driven automatic bank change is not a user content edit
  and must not resize the outer window.
- Empty-axis scale changes use `Makie.defaultlimits` and validate those limits
  before mutation. A temporary linear transition lets the native empty-data
  lifecycle adopt a valid new default interval, including strict positive log.
  No scientific values or eligibility rules change.
- The earlier geometry check failed on a missing native trailing track; its
  allocation was fixed immediately. The numeric run that loaded that intermediate
  source is a failed development run, not numerical regression evidence or a
  passing gate. The fresh complete suite validates the corrected source.
- The completed documentation build `/tmp/plotbuilder-docs-final.log` passed all
  example, cross-reference, and missing-doc checks. Documenter reported its normal
  large-HTML PNG fallback and search-index size warnings. A fresh build follows
  the final geometry/export corrections.
- Ignored historical planning notes under `docs/notes` were inspected as history
  and left unchanged. Supported public documentation and maintained callers
  use the current API exclusively.

Final dispatch and backend checks:

- Split the observed scalar and vector `Makie.plot` forwarding methods. Their
  previous union signature was ambiguous with the supported primary collection
  method for a vector of observations. The surviving methods forward directly
  to public `plot`; ordinary native array plotting retains its native dispatch.
- Selectors `'observed public forwarding' 'uncertainty legend shading follows'`
  passed **17/17** assertions, two maintained items, 171.67 seconds;
  `/tmp/plotbuilder-dispatch-shading-final.log`. This includes actual scalar,
  vector and tuple observed calls, native array dispatch, and pixel equality
  after a guide notification preserves a hidden legend entry's appearance.
- Repeated final GL checks passed **42/42** assertions in each import order:
  `/tmp/plotbuilder-closure-gl.log` and
  `/tmp/plotbuilder-closure-gl-first.log`. Both execute actual display closure,
  retained-handle mutation/export and redisplay with once-only callbacks.
- Final WGL smoke check built all **16** gallery panels;
  `/tmp/plotbuilder-closure-wgl.log`. The separate localhost embedding check is
  recorded below; smoke construction alone is not browser interaction evidence.
- `git diff --check` passed. The residue search found removed option names only
  in explicit rejection diagnostics, negative tests, this execution history,
  and unrelated scientific uses of the ordinary word `blocks`.
- Final documentation command
  `JULIA_DEPOT_PATH=/tmp/observed-plan-depot:/home/amartins/.julia julia --project=docs --startup-file=no docs/make.jl`
  exited successfully; `/tmp/plotbuilder-closure-docs.log`. The authoritative
  Literate source and generated plotting page were reviewed. Documenter retained
  its large-HTML image fallback and search-index size warnings; no failed
  examples, cross-references, or missing public docstrings remain.
- An additional direct GL statistical script check exposed its synthetic
  fixture's missing explicit `Measurements` import, before plotting began.
  Added that import to `dev/plotting/manual_gl_monte_carlo.jl`; the package's
  optional dependency boundary remains unchanged. Initial failure:
  `/tmp/plotbuilder-gl-consumers.log`; corrected run recorded below.
- The corrected GL consumer command used the backend environment above and
  `LINECABLEMODELS_GL_GALLERY_SMOKE=true`, with
  `julia --project=. --startup-file=no -e 'include("dev/plotting/manual_gl_monte_carlo.jl"); include("dev/plotting/manual_gl_cable_collection.jl")'`.
  It exited successfully: both statistical handle/axis assertions and all six
  cable-collection assertions passed. Five retained statistical views and the
  automatic 2×3 / explicit 1×4 physical previews were constructed;
  `/tmp/plotbuilder-gl-consumers-final.log`. This additional command was a
  construction check; the separate 42-assertion gates exercised real windows.

WGL embedding and browser evidence:

- Started the existing gallery with the installed backend environment,
  `LINECABLEMODELS_WGL_MANUAL=true`,
  `LINECABLEMODELS_WGL_ARTIFACTS=/tmp/plotbuilder-wgl-embedding`, and
  `LINECABLEMODELS_WGL_PORT=18082`:
  `julia --project=. --startup-file=no dev/plotting/manual_wgl.jl`.
- `curl --fail --silent --show-error --max-time 60 http://127.0.0.1:18082/ -o /tmp/plotbuilder-wgl-gallery.html`
  succeeded. The existing Bonito app served all 16 figure cards.
- An isolated Chrome profile used native headless WebGL with SwiftShader:
  `google-chrome --headless=new --no-first-run --no-default-browser-check --user-data-dir=/tmp/plotbuilder-chrome-profile --remote-debugging-address=127.0.0.1 --remote-debugging-port=19223 --enable-unsafe-swiftshader --use-gl=angle --use-angle=swiftshader --disable-background-networking about:blank`.
- Disposable Node scripts used the installed Node WebSocket client and the
  localhost Chrome debugging protocol, without installing dependencies:
  `node /tmp/plotbuilder-wgl-browser.mjs`,
  `node /tmp/plotbuilder-wgl-browser-ready.mjs`,
  `node /tmp/plotbuilder-wgl-actions.mjs`, and
  `node /tmp/plotbuilder-wgl-log.mjs`.
- All **16 canvases initialized**, with zero remaining native loading spinners
  and **zero JavaScript exceptions**. Inspected the rendered screenshot at
  `/tmp/plotbuilder-wgl-browser.png`. The browser's only resource error was the
  gallery's absent `/favicon.ico`, not a plotting asset.
- Actual browser mouse actions changed the first four-panel page from log X to
  linear X, reset its view, and returned it to log X. Inspected the resulting
  coordinates, ticks, unchanged Y presentation, and status text in
  `/tmp/plotbuilder-wgl-reset.png` and
  `/tmp/plotbuilder-wgl-log-settled.png`. The latter action did not use reset.
  The other displayed quantity retained its original logarithmic view.
- Browser evidence is limited to the stated initialization and actions under
  headless Chrome/SwiftShader; it does not establish every interaction in every
  browser or hardware WebGL implementation. Closed the owned browser and the
  temporary localhost server after inspection.

## Final integrated closure

The complete visual command exited successfully:

```bash
JULIA_DEPOT_PATH=/tmp/observed-plan-depot:/home/amartins/.julia \
JULIA_LOAD_PATH='@:/home/amartins/Documents/KUL/LineCableModels/test/visual:@stdlib' \
LINECABLEMODELS_TEST_PLOTTING=true \
julia --project=. --startup-file=no test/runtests.jl tag:visual
```

**2,884/2,884 assertions**, 62 maintained items in 18 files, 1,848.26 seconds;
`/tmp/plotbuilder-closure-visual.log`. This final run replaces the earlier failed
development runs as the complete visual gate. The forwarding regression was
added after that run started and passed separately: **17/17** assertions including
its 12 new dispatch checks and five repeated legend-shading checks. These totals
overlap and must not be summed as distinct coverage.

Selected final suite sections:

| Evidence | Passing assertions |
| --- | ---: |
| Numeric axes, scales, formatting, precision and rollback | 391 |
| Native guide creation/restoration/alignment | 128 |
| Matrix capacity, residual dimensions and frame calibration | 194 |
| Preview controls and callable/composite widget lifetime | 66 |
| Preview materials, patterns and viewport-following earth | 173 |
| Selection preservation and SVG success/failure restoration | 192 |
| Responsive geometry, native edits and local flow identities | 56 |
| Native series styling | 56 |
| Uncertainty projection, bounds, sampling and legend interaction | 872 |
| Retained UQ overlays/statistics/bands | 231 |
| Explicit retained statistical verbs | 15 |
| Ordinary/benchmark reporting resolution and raw-value preservation | 79 |
| Full modal/rectangular coordinates, ragged samples, assemblies and index views | 72 |

The remaining sections cover backend dispatch, native primitives/material scales,
source-independent formulation overlays and scientific grouping. Their exact
test-item results are in the complete log and the feature map above.

Additional completed gates:

- Observation/report/core-extension boundary: **132/132** assertions in 10 items.
- Dispatch/legend focused regression: **17/17** assertions in two items.
- GL normal and backend-first actual-window checks: **42/42 each**.
- GL retained statistical and physical-preview construction: **8/8** assertions.
- WGL: 16 native gallery panels, all 16 browser canvases initialized, actual
  scale/reset actions, no JavaScript exceptions; scope described above.
- Documentation doctest: **1/1**; full `docs/make.jl`: successful.
- Final `git diff --check`: successful. No Project or Manifest changes.

| Stage | Surviving callable/owner | Competing path removed | Executed gate |
| --- | --- | --- | --- |
| 0 | Existing observations/native APIs | Complete predecessor inventory identified | Clean baseline and installed environment verified |
| 1 | Typed public observed `plot`, Grammar retained products | Catch-all page orchestration, first-product assumptions, plot-side bands, pairing/capacity duplication | Retained interpretation, report/source boundaries, modal/assembly/index and band checks |
| 2 | Shared numeric and native series operations | Signed-log permissions, reference-dependent defaults and shifted slots | 391 axis and 872 uncertainty assertions, styles/UQ overlays |
| 3 | One live UIPlot and shared widget allocation | Provisional handle and separate callback/toolbar construction | Callable/composite tests plus actual GL close/redisplay |
| 4 | One guide composition operation | Initial/live docks and managed anchor parsers | 128 guide assertions, endpoint geometry, native edits |
| 5 | Single capacity and local frame fitting | Whole-canvas aspect, residual dummy tracks, sibling budgets, landscape/GLFW lock | 194 matrix and 56 responsive assertions, physical previews |
| 6 | Shell presentation plus SVG writer | Hard-coded chrome rows and matrix export exceptions | 192 selection/export assertions and GL late-Cairo checks |
| 7 | Complete supported consumers and public docs | Obsolete caller keywords and superseded documentation | Full visual/observation suites, retained GL scripts, docs build |
| 8 | Reviewed unified call graph | Dormant methods/imports/state/aliases and duplicate authority | All final gates above and residue/diff audit |

Intentionally replaced assertions were limited to superseded requirements:
quantity pairing, full-size residual padding, forced landscape dimensions,
recipe-dependent adaptive-log rejection, implicit interval/ellipsis choices,
and absence of empty marker handles that preserve native legend identity.
Projected coordinates, full uncertainty bounds, native visibility, physical
aspect, source independence, and exact export restoration remain tested.

No implementation failures remain in the executed gates. The documented limits
are browser/hardware coverage beyond the stated WGL check and ordinary
Documenter size warnings. No solver campaign, dependency upgrade, branch change,
commit, push, reset, or historical compatibility layer was introduced.
