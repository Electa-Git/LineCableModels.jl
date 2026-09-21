# Manual scripts after the engine / ObservedResult / PlotBuilder refactors

Audit date: 2026-09-21. Checkout: `release/v0.2.0`, base HEAD
`4f72bf45aeb4558e0fd48cd5f45e7a84d86f478b`, with the PlotBuilder changes in the worktree.
This inventory includes Git-ignored scripts. No campaign archive was changed.

## Running the scripts

Use the existing Julia environment in the IDE. The Gauntlet runners add its local
project to `LOAD_PATH`; they do not activate a different project or install packages.
For the installed checkout, the checks used this stack from the repository root:

```bash
export JULIA_DEPOT_PATH=/tmp/observed-plan-depot:/home/amartins/.julia
export JULIA_LOAD_PATH="@:$PWD/gauntlet:$PWD/test/visual:$PWD/dev/plotting:@v1.12:@stdlib"
julia --project=. --startup-file=no dev/inspect_saved_gauntlet.jl
```

The inspectors expose `plot_backend=:gl` / `:cairo`, `make_plots`, and
`display_plot` near the top. GL needs a display. `enable_svg_export=true` explicitly
loads Cairo. `show_full_report=true` prints the full-width report; all tables are
available in the IDE regardless of this display setting. Default legends show
owner-described differences; `series_labels` remains an optional explicit override.

Before/now comments accompany changed calls. In particular:

- Result plots use `fig_size`; previews and `plotwindow` keep `size`.
- `layout` is the sole panel capacity. Quantities have separate figure families.
- Retained reports/observations support compatible display-unit changes directly:
  `plot(report; ydata=(R,X,G,B), length_unit=:base)`. The observation owner
  re-expresses retained curves; clipping decisions and recorded RMS units remain fixed.
- Raw numerical conveniences construct observations and delegate to observed plotting.
- Results carry `ComputationDetails`, not a plain details NamedTuple.
- XLSX exports return one filename per quantity; append those paths to a flat list.
- `formula(:unified; options=(Γ=...,))` owns prescribed Γ under `options`.

Several personal runners remain Git-ignored as before and were repaired in place.
The shared inspector reader and this inventory are explicitly unignored so the
tracked inspectors do not acquire an omitted dependency.

The direct compute-to-plot example needs only the package and GLMakie. In the
installed checkout, this invocation keeps a REPL open with its plot windows:

```bash
JULIA_DEPOT_PATH=/tmp/observed-plan-depot:/home/amartins/.julia \
JULIA_LOAD_PATH='@:/home/amartins/Documents/KUL/LineCableModels/dev/plotting:@stdlib' \
julia --project=. --startup-file=no -i dev/run_two_bare_wires.jl
```

The same file can be included in an IDE session that already provides those
dependencies. Its declarations, computation, four quantity figures, controls,
and close/edit/redisplay were executed on GL; it leaves the active project and
`LOAD_PATH` unchanged. It does not read saved runs or execute a FEM/MC campaign.

## File disposition

| File under `dev/` | Current operation and check |
| --- | --- |
| `inspect_saved_gauntlet.jl` | Reads saved numerical operands; explicitly builds current comparisons/observations/tables. Default PSCAD archive produced 1,620 term rows without solver execution. The full 81-panel R view passes with automatic scientific labels; earlier two-quantity/four-coefficient scale/reset checks also pass. |
| `inspect_saved_gauntlet_mc.jl` | Same path for full UQ results. Default uses the readable 320 kV armoured DC case; the previous 30 kV default requires regeneration. The full 36-panel R view passes with owner-supplied LEP/Monte Carlo labels. Retained statistics (58,176 values), comparisons (288 rows), and two-quantity/four-coefficient scale/reset checks also pass. |
| `inspect_saved_gauntlet_inputs.jl` | Shared checksummed numerical reader and saved timing handoff. Rejects missing attempts and removed marginal-only records; previous attempts require explicit selection. |
| `run_two_bare_wires.jl` | Direct declaration → `compute(problem, formulations)` → `plot(results; length_unit=:base)`. Two buried bare copper wires with constant versus Longmire soil properties; automatic labels/colors, no Gauntlet/FEM or environment mutation. |
| `run_two_bare_wires_fem.jl` | Current parametric constructors; flattened per-quantity XLSX filenames. Constructors, full CSV writes and four actual XLSX workbooks checked with retained fixtures. New FEM sweeps were not run. |
| `run_18kv_trefoil_fem.jl` | Preserves active project; current case/formulation/options constructed. New FEM solve was not run. |
| `run_quasi_full.jl` | Calls the current FEM voltage-path owner directly. Problem/model/mesh-plan construction checked. Meshing and GetDP execution were not run. |
| `run_three_bare_wires_baseline.jl` | Removed mandatory dependency on a historical local plan; keeps actual runner/fixture capture. Synthetic 3×3×9 CSV round trip and checksums checked. No baseline campaign run. |
| `run_line_parameters.jl` | Revise is optional. Default analytical case executed: finite 9×9×101 Z/Y arrays. FEM branch was not run. |
| `benchmark_shunt_models.jl` | Uses current catalogue `.problem`, so edited frequency overrides are applied; preserves project. Both public coaxial cases executed, and an explicit three-frequency override was verified; no performance improvement is claimed. |
| `benchmark_local_shunt.jl` | Current Γ options and seven-argument private completion call. Boundary preparation, diagnostics and a three-frequency completion executed. Full timing/comparison campaign was not repeated. |
| `prototype_wire_screen.jl` | Ordinary result collection replaces fabricated ParametricResult; candidate/reference labels follow input order. Current geometry and independent kernel self-checks executed. Full six-refinement experiment was not run. |
| `prototype_tape_element.jl` | Current independent mathematical helper, retained. Its self-checks execute through the parent prototype; it is not a standalone runner. |
| `gridspace_test.jl` | Current declarative Grid syntax; two-point problem construction and materialization checked. |
| `surprised_pikachu_cable.jl` | Current polygon preview; actual GL window checked. This remains preview geometry, not a coaxial solver claim. |
| `verify_getdp_artifact.jl` | Current archive checker; parsed and inspected. Network downloads were not repeated for this API migration. |
| `plotting/_gallery_fixtures.jl` | Current synthetic/retained fixtures shared by GL and WGL galleries. |
| `plotting/manual_gl.jl` | Native controls/window/export gate; normal and backend-first runs passed in the refoundation validation. |
| `plotting/manual_gl_comparison.jl` | Corrected `ComputationDetails` and `fig_size`; actual four-window GL execution and assertions passed. |
| `plotting/manual_gl_gallery.jl` | Current gallery construction; retained from the validated refoundation. |
| `plotting/manual_gl_cable_collection.jl` | Current automatic 2×3 and explicit 1×4 previews; validated in the refoundation. |
| `plotting/manual_gl_monte_carlo.jl` | Uses a completed synthetic UQ fixture; no fresh MC run. Native statistical views validated in the refoundation. |
| `plotting/manual_wgl.jl` | Current browser embedding; browser initialization/actions validated in the refoundation. |

Deleted:

- `dev/test_manual_gauntlet.jl`: orphaned pre-refactor test-case loader and removed
  `test/gauntlet` paths. Its underground formula sweep also invokes current
  `not yet implemented` Pollaczek/Saad branches; replacing those would change the
  experiment, not merely its syntax. Current Gauntlet entry points are the runners/inspectors above.
- `dev/quasi_full_paths.jl`: redundant forwarding helper; its caller now invokes
  the existing extension operation directly.

`gauntlet_progress_plan.md` is a historical design note, not a runnable example.
Existing binary exports, captured solver data, dependency manifests and prior
execution notes were not removed as obsolete code.

## Do saved Gauntlet runs need regeneration?

**A complete FEM/PSCAD campaign rerun is not required by the plotting change.**
The local `gauntlet/.work/all-references` inventory contains 62 benchmark states:
61 selected saved attempts and one without a completed attempt. All 61 selected
analysis snapshots lack `observed_data`; they need new comparisons/observations
before current retained-report loading. The inspectors do this explicitly in memory
from readable saved operands without writing an analysis or launching a solver.

| Saved numerical representation | Cases | Required work |
| --- | ---: | --- |
| `gauntlet_calculation` / `gauntlet_result_space` | 41 | Current supported representation. Confirmed reuse for the default 220 kV PSCAD and two-bare-wire FEM cases; rebuild their comparisons/reports. Other cases are classified from headers, not a completed bulk read. |
| `gauntlet_uncertainty` | 8 | Same; full UQ readers pass with `Measurements` loaded. No resampling required for those retained results. |
| `gauntlet_moments` | 12 | Regenerate UQ calculations for the current workflow. These files retain only R/L/G/C marginal mean/std arrays, not a current full scientific result or its uncertainty dependencies. |
| No selected completed attempt | 1 | `benchmark_220kv_milliken_1x2500_252_trefoil_lep_montecarlo` needs an initial completed calculation. |

The twelve marginal-only cases are:

- `benchmark_132kv_630mm2_flathor_lep_montecarlo`
- `benchmark_132kv_cigre_tb880_case0_630cu_trefoil_lep_montecarlo`
- `benchmark_18kv_1000mm2_trefoil_homogenized_lep_montecarlo`
- `benchmark_18kv_1000mm2_trefoil_lep_montecarlo`
- `benchmark_220kv_eaxecew_1x2500_252_trefoil_lep_montecarlo`
- `benchmark_30kv_na2xs2y_630mm2_trefoil_lep_montecarlo`
- `benchmark_380kv_2000mm2_flatver_lep_montecarlo`
- `benchmark_525kv_1600mm2_bipole_lep_montecarlo`
- `benchmark_640kv_2000mm2_bipole_lep_montecarlo`
- `benchmark_solid_1000mm2_single_lep_montecarlo`
- `benchmark_two_bare_wires_lep_montecarlo`
- `benchmark_two_insulated_wires_lep_montecarlo`

These statements concern archive readability and retained information. Reusing
old numerical operands does not claim they were computed by today's formulas.
Several top-level states are pending/failed/running while retaining an older
`current` attempt; neither inspector silently treats that older attempt as the
latest successful campaign. Set `use_previous_complete=true` deliberately when
that is the intended source. No solver liveness is inferred from those state words.

For a persistent new analysis, use the current `Gauntlet.compare_saved` operation
with explicit numerical-file bindings and checksums; its TOML format and CLI are
documented in `gauntlet/README.md` under saved comparisons. This differs from
`run --force`, which actually recomputes the selected benchmark.

For example, when ready to regenerate the old 30 kV UQ calculation into a new
campaign directory (this command was **not executed** during migration):

```bash
./gauntlet/lcm gauntlet run \
  --definition gauntlet/benchmarks/uq/benchmark_30kv_na2xs2y_630mm2_trefoil_lep_montecarlo.jl \
  --directory gauntlet/.work/current-uq
```

Then set the inspector's `campaign_directory` and `benchmark_id` to that run.
The current default remains an existing readable case so the inspector can be
used immediately without an automatic campaign launch.

## Examples, Literate sources and evidence

- All 23 current dev Julia files, three tutorials, two Literate sources and two
  Gauntlet example definitions parsed: 30 checks.
- Constructor/export checks: 16 passing assertions; prototype checks: 40 passing.
- Normal and UQ inspector plots both passed; the rendering checks selected the first 2×2 coefficient population and two quantities to keep window construction bounded. Full comparison/statistical tables were built.
- All three tutorials executed through their final plots and exports using byte-identical
  copies under `/tmp/dev-tutorials-current`. Tracked example output artifacts were
  preserved. Tutorial 2/3 comments and names now use explicit `ObservedResult`.
- `docs/literate/plotting.jl` explains old/new layout syntax.
- `docs/literate/gauntlet.jl` documents retained-only selection, explicit Cairo
  loading, optional physical-point filters and saved-analysis regeneration.
- Authoritative Literate generation ran through the existing `docs/make.jl`
  generation functions. `docs/doctest.jl`: 1/1 passed.
- Earlier complete graphics/docs validation is recorded in
  `plotting/refoundation_execution.md`; it is not substituted for the actual
  script checks listed here.

The disposable checks used the environment above, with these exact commands:

```bash
julia --project=. --startup-file=no /tmp/check_dev_scripts.jl
julia --project=. --startup-file=no /tmp/check_manual_inspector.jl
julia --project=. --startup-file=no /tmp/check_inspector_mc.jl
julia --project=. --startup-file=no /tmp/check_inspector_plots.jl
julia --project=. --startup-file=no /tmp/check_inspector_uq_plots.jl
julia --project=. --startup-file=no /tmp/check_dev_gl.jl
julia --project=. --startup-file=no /tmp/check_dev_analytical.jl
julia --project=. --startup-file=no /tmp/check_shunt_override.jl
julia --project=. --startup-file=no /tmp/check_baseline_helpers.jl
julia --project=. --startup-file=no /tmp/check_tutorial_files.jl
julia --project=. --startup-file=no /tmp/check_saved_uq.jl
julia --project=. --startup-file=no /tmp/check_saved_numeric.jl
julia -O0 --project=. --startup-file=no /tmp/check_saved_numeric.jl
julia --project=docs --startup-file=no /tmp/generate_dev_docs.jl
julia --project=docs --startup-file=no docs/doctest.jl
```

Logs are `/tmp/dev-script-checks.log`, `/tmp/dev-inspector-check.log`,
`/tmp/dev-mc-inspector-check.log`, `/tmp/dev-inspector-plots-final.log`, `/tmp/dev-inspector-uq-plots-final.log`,
`/tmp/dev-gl-script-check.log`, `/tmp/dev-analytical-check.log`, `/tmp/dev-shunt-override.log`,
`/tmp/dev-baseline-helpers.log`, `/tmp/dev-tutorials-final.log`,
`/tmp/dev-uq-readers.log`, `/tmp/dev-saved-numeric.log`, `/tmp/dev-saved-numeric-o0.log`,
`/tmp/dev-literate-generation.log`, and `/tmp/dev-docstrings-check.log`.
The complete archive header inventory is `/tmp/dev-gauntlet-inventory.tsv`.

Validation limitations: no new FEM/GetDP/PSCAD campaign, no new MC sampling,
no full prototype refinement/timing campaign, and no repeat network artifact
verification. Full GL/WGL/manual interaction coverage remains the earlier
refoundation evidence except the two actual GL scripts rerun here.

The first disposable harness attempts exposed Julia world-age errors in the
harness itself (after completed script work), and an incorrect assumption about
example output filenames. Subsequent checks use top-level inclusion and actual
script output paths. The analytical log contains a world-age error only in its
final, separate Gauntlet-definition check, after the analytical and local-shunt
script checks completed. Those two definitions subsequently passed with
`julia --compile=min --project=. --startup-file=no /tmp/check_other_callers.jl`.
That interpreter setting is not suitable for the JLD2 archive sweep: it produced
`invalid struct allocation` inside JLD2; it is not evidence of damaged solver data.

The broader deterministic archive sweep was stopped during prolonged JLD2/compiler
work (including an `-O0` retry). It is **not** counted as passing evidence for all
41 deterministic cases. The two representative deterministic cases and all 16
operands of the eight full-UQ cases were read successfully in ordinary compilation
mode. A remaining deterministic record that cannot be read needs individual
inspection; its header alone does not establish a need to rerun its solver.

Final syntax scan: 30/30 (`/tmp/dev-final-parse.log`). `git diff --check` passed.
No commits, pushes, dependency updates or campaign regeneration were performed.
