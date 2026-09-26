# Manual tests

This tree contains human-operated investigations only. `JuliaTestItems.toml`
excludes it from automated test discovery, so TestItemRunner neither parses nor
runs these files. Nothing here is a CI gate or a maintained feature API.

Use the active Julia environment when including scripts from an IDE or REPL.
The scripts do not select or install a project on the caller's behalf, except
that `plotting/` has an explicit environment for its optional backends.

| Area | Purpose |
| --- | --- |
| `performance/` | Local timing and numerical inspection of shunt models. |
| `calculations/` | Editable analytical, grid-space, export, and FEM comparison runs. |
| `fem/` | Direct Gmsh/GetDP investigations and baseline captures. |
| `plotting/` | Human-inspected GL/WGL galleries and plotting experiments. |
| `prototypes/` | Disposable scientific prototypes that do not extend production code. |
| `pscad/` | A manual PSCAD station run and its private local configuration. |
| `verification/` | Manually invoked packaging checks. |

The filename states the runnable experiment. In particular:

Run `calculations/modal_analysis.jl` with the existing
`test/manual/plotting` project when the active REPL project does not already
provide GLMakie.

| Script | Manual purpose |
| --- | --- |
| `performance/benchmark_local_shunt.jl` | Inspect local-shunt preparation, allocation, and selected saved references. |
| `performance/benchmark_shunt_models.jl` | Compare warmed public-API shunt-model timings. |
| `calculations/gridspace.jl` | Inspect a small constructed `Gridspace`. |
| `calculations/run_line_parameters.jl` | Run one editable catalogue case through the analytical or FEM backend and export PSCAD input. |
| `calculations/run_two_bare_wires.jl` | Compare two analytical earth-property formulations interactively. |
| `calculations/run_two_bare_wires_fem.jl` | Sweep the two-wire FEM case and inspect CSV/XLSX/plot output. |
| `calculations/modal_analysis.jl` | Direct REPL study of nine-terminal armored cables through phase, modal, segment, observation, report, and interactive GLMakie APIs. |
| `fem/run_18kv_trefoil_fem.jl` | Run the current 18 kV catalogue case through quasi-TEM FEM. |
| `fem/run_quasi_full.jl` | Exercise the disposable quasi-full GetDP formulation directly. |
| `fem/run_three_bare_wires_baseline.jl` | Capture analytical/FEM baseline data for human comparison. |
| `prototypes/prototype_wire_screen.jl` | Investigate the disposable wire-screen boundary model. |
| `plotting/surprised_pikachu_cable.jl` | Exercise arbitrary-polygon preview geometry. |
| `pscad/run_pscad.jl` | Run a two-wire problem on the configured PSCAD station. |
| `verification/verify_getdp_artifact.jl` | Download and verify the packaged GetDP archives. |

The remaining plotting scripts are described in `plotting/README.md`.

## Output policy

Generated output does not belong in the versioned tree. Scripts that explicitly
write reports or images default outside the checkout, under
`joinpath(tempdir(), "linecablemodels-manual", ...)`. Set
`LINECABLEMODELS_MANUAL_OUTPUT` to an absolute directory when the output should
survive temporary-directory cleanup. Backend-specific plotting variables remain
available as overrides and must likewise point outside the checkout.

FEM and PSCAD engines may retain their own ignored runtime directories according
to their explicit execution configuration; those directories are runtime state,
not manual-test fixtures.

For PSCAD, edit the ignored `pscad/local-pscad.toml`, then run
`include("test/manual/pscad/run_pscad.jl")` from the repository root. A
sanitized configuration template is available at
`examples/pscad/remote.example.toml`. Each execution creates a fresh run in the
configured exchange directory and leaves `result`, `Zresult`, and `Yresult`
available for inspection.

See `plotting/README.md` for the interactive plotting procedures.
