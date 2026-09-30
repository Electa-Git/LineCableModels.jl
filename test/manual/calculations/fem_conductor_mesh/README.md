# Conductor mesh fixtures and comparisons

[fixtures.jl](fixtures.jl) defines passive screen, tubular-shell and sector
constructors. [fixtures.md](fixtures.md) records their dimensions, simplifications
and historical qualification. Including the constructors performs no computation.

[run_cable_fixtures.jl](run_cable_fixtures.jl) uses the public Julia `compute`
and `export_data(:onelab, ...)` APIs. It runs serially, using packaged native
executables without an external wrapper. Its default is quasi-fw at 0.1, 50,
10 000 and 1 000 000 Hz, plus one explicitly refined 1 MHz comparison per design.
These are prescribed comparisons, with no automatic acceptance or refinement.
The existing 24-skin-depth / 192-layer mesh preset is retained; this cleanup
does not qualify a cheaper preset.

From the repository root, choose a new output directory:

```sh
julia --project=. test/manual/calculations/fem_conductor_mesh/run_cable_fixtures.jl /tmp/cable-fixtures > /tmp/cable-fixtures.log 2>&1
tail -F /tmp/cable-fixtures.log
```

Restart the same command after interruption. The backend checks compatibility
before resuming retained solves. Existing detached exports are left untouched
so caller edits survive. Use a new output directory when changing inputs.

For a smaller selection, include the file and prescribe controls explicitly:

```julia
include("test/manual/calculations/fem_conductor_mesh/run_cable_fixtures.jl")
run_cable_fixtures("/tmp/cable-fixtures-50hz";
    frequencies=[50.0], refined_frequencies=Float64[])
```

The `physics`, `mesh_options` and `refinements` keywords permit explicit studies.
The detached bundle for each design is `<output>/<design>/detached/study.pro`;
it contains the normal scan and runs in ONELAB without Julia. Each managed scan
saves raw complex Z/Y in `matrices.csv`, primitive Z/P tables, prescribed options
and wall/compilation timings. No values are clipped or signs corrected.

[plot_cable_fixtures.jl](plot_cable_fixtures.jl) reads saved normal/refined CSVs
through the public observation/plotting API and exports R/G/B without solving:

```sh
julia --project=. test/manual/calculations/fem_conductor_mesh/plot_cable_fixtures.jl /tmp/cable-fixtures
```

The plotter supports both historical physics selections and exports missing
SVG/PNG files only. [plot_bare_wire_mvp.jl](plot_bare_wire_mvp.jl) and
[plot_spectrum.jl](plot_spectrum.jl) likewise remain useful readers of the saved
bare-wire and spectrum evidence; they launch no FEM computations.

The completed campaign's [delivery assessment](delivery-assessment.md) and
[execution journal](../fem_conductor_mesh.md) retain its scientific findings and
limitations. Obsolete campaign launchers and Python tooling were retired during
the [development cleanup](../fem_development_cleanup.md); recorded results remain
under `.linecablemodels/fem/conductor-mesh-qualification/`.
The original procedure and experiment descriptions remain in
[qualification-history.md](qualification-history.md) as a historical record.

Automated native geometry, integration, export and resume tests remain under
`test/extensions/fem_*.jl`. This manual directory is excluded from automatic
test discovery and contains no production acceptance policy.
