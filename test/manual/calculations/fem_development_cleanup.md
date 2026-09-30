# FEM development cleanup — 2026-09-29

The maintained FEM implementation uses Julia orchestration and shared native
Gmsh/GetDP sources. Managed computation and detached `:onelab` export remain
supported. This cleanup changes development tooling and documentation only;
production equations, mesh defaults, integration, source-column reuse, runtime
resume and automated regression tests are preserved.

## Retained entry points

| Purpose | Entry point |
|---|---|
| Interactive two-wire study and one representative detached export | [run_two_bare_wires_fem.jl](run_two_bare_wires_fem.jl) |
| Three-conductor source/receiver comparisons across layers | [run_three_bare_wires_baseline.jl](../fem/run_three_bare_wires_baseline.jl) |
| Screen, tubular shell and sector comparisons | [Cable fixture guide](fem_conductor_mesh/README.md) |
| Explicit PML count comparison with raw component and cost reporting | [run_fem_fixed_pml.jl](run_fem_fixed_pml.jl) |
| Small detached export | [export_onelab_toy.jl](../fem/export_onelab_toy.jl) |
| Native perfect-conductor diagnostic | [run_quasi_full.jl](../fem/run_quasi_full.jl) |
| Detached execution and caller ownership checks | [validate_onelab_export.jl](../fem/validate_onelab_export.jl) |

The cable runner now directly calls public `compute` and `export_data(:onelab)`.
It no longer includes the old MVP campaign, requires a solver wrapper, patches
native sources, or repeats completed native parity/source-edit campaigns.
It defaults to quasi-fw, with explicit frequency and mesh-control keywords.
Normal and refined controls remain independent prescribed computations.
Raw CSV output remains compatible with the retained saved-data plotters.
The 192-layer preset is unchanged; cleanup supplies no new accuracy claim.

The standalone Julia plotters for the cable, bare-wire, spectrum and historical
voltage-reference CSVs remain available. They do not launch FEM jobs. The
automated tests under `test/extensions/` retain native geometry, conductor
grading, integration, algebra, export, import, worker and resume coverage.
Qualification and scientific acceptance remain outside production.

## Retired tools and preserved evidence

Removed 53 obsolete executable files: 30 Python scripts and 23 associated
Julia campaign/export launchers, Gmsh prototypes and GetDP experimental problems.
They covered retired voltage-preparation paths, source-patched boundary trials,
Cartesian compactification/static-boundary trials, isolated wave controls and
completed conductor qualification grids. Current tests keep the native
integration fixtures under `test/fixtures/data/fem/`.

All scientific notes are retained. Historical records now identify their source
era and retired commands. The old conductor qualification guide is preserved in
[qualification-history.md](fem_conductor_mesh/qualification-history.md);
its replacement README documents current Julia usage. Current PML controls and
timing evidence remain in [fem_fixed_pml_controls.md](fem_fixed_pml_controls.md),
[fem_quasi_fw_performance.md](fem_quasi_fw_performance.md) and
[fem_pml_resolution_diagnostics.md](fem_pml_resolution_diagnostics.md).

No existing computed runs, meshes, exports, matrices, logs or plots under
`.linecablemodels` were removed or rewritten. Because many retired files were
uncommitted, their exact contents were archived locally before removal:

```text
.linecablemodels/fem/development-cleanup-20260929/
    before.tar.gz          original manual files and plans
    before.sha256          original file checksums
    archive.sha256         archive checksum
    retired-files.txt      exact removal manifest
    preserved.sha256       production and automated-test checksums
    status-before.txt      pre-cleanup Git status
    check.log              bounded cleanup verification
```

The archive was compared with the originals before deletion. It is ignored
local evidence, not a dependency of the maintained code. To inspect the old
tools without overwriting current files:

```sh
mkdir -p /tmp/lcm-fem-history
tar -xzf .linecablemodels/fem/development-cleanup-20260929/before.tar.gz -C /tmp/lcm-fem-history
```

This cleanup does not implement arbitrary finite complex Gamma. The existing
quasi-fw formulation and its approximation remain unchanged. Historical failed
approaches are retained as evidence, not alternative production implementations.

## Cleanup verification

All 16 retained Julia files under the FEM/calculations directories passed native
JuliaSyntax parsing; manual scripts remain excluded from automated discovery.
The current guide's local file links resolve. Production and automated-test
file checksums match their pre-cleanup values.

The simplified cable runner completed three serial 50 Hz quasi-fw solves
(screen, tube, sector; nine source columns), exported three detached bundles,
and wrote the expected raw Z/Y CSVs, primitive P tables and compilation timings.
This execution check explicitly used eight PML layers and two skin depths to
bound cost; it provides no new scientific accuracy evidence. No historical
qualification grid was rerun. Logs and the exact verification callers are in
the cleanup evidence directory above.
