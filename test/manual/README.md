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
| `calculations/` | Editable analytical, grid-space and export runs. |
| `plotting/` | Human-inspected GL/WGL galleries and plotting experiments. |
| `prototypes/` | Disposable scientific prototypes that do not extend production code. |
| `pscad/` | A manual PSCAD station run and its private local configuration. |

The filename states the runnable experiment.

Run `calculations/modal_analysis.jl` with the existing
`test/manual/plotting` project when the active REPL project does not already
provide GLMakie.

| Script | Manual purpose |
| --- | --- |
| `performance/benchmark_local_shunt.jl` | Time and inspect the boundary shunt stage of one known case. |
| `performance/benchmark_shunt_models.jl` | Compare warmed public-API shunt-model timings. |
| `calculations/gridspace.jl` | Inspect a small constructed `Gridspace`. |
| `calculations/run_line_parameters.jl` | Run one editable datasheet case through the analytical backend and export PSCAD input. |
| `calculations/run_two_bare_wires.jl` | Compare two analytical earth-property formulations interactively. |
| `calculations/modal_analysis.jl` | Direct REPL study of nine-terminal armored cables through phase, modal, segment, observation, report, and interactive GLMakie APIs. |
| `prototypes/prototype_wire_screen.jl` | Investigate the disposable wire-screen geometric boundary model. |
| `plotting/surprised_pikachu_cable.jl` | Preview geometry with arbitrary polygons. |
| `pscad/run_pscad.jl` | Run a two-wire problem on the configured PSCAD station. |

The remaining plotting scripts are described in `plotting/README.md`.

## Output locations

Generated output does not belong in the versioned tree. Scripts that explicitly
write reports or images default outside the checkout, under
`joinpath(tempdir(), "linecablemodels-manual", ...)`. Set
`LINECABLEMODELS_MANUAL_OUTPUT` to an absolute directory when the output should
survive temporary-directory cleanup. Backend-specific plotting variables remain
available as overrides and must point outside the checkout.

A PSCAD run can leave ignored runtime directories, as its execution configuration
sets. Those directories are runtime state, not manual-test fixtures.

For PSCAD, edit the ignored `pscad/local-pscad.toml`, then run
`include("test/manual/pscad/run_pscad.jl")` from the repository root. A
sanitized configuration template is available at
`examples/pscad/remote.example.toml`. Each execution creates a fresh run in the
configured exchange directory and leaves `result`, `Zresult`, and `Yresult`
available for inspection.

See `plotting/README.md` for the interactive plotting procedures.
