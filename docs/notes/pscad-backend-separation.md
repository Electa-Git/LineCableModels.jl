# PSCAD backend and Gauntlet separation

Implementation on `release/v0.2.0-pscad-ext`, based on checkpoint
`2d694a2444f52942e9d14cd5ff835760295ff532`. The other release worktrees are
independent. This change does not add scientific formulations.

## Ownership and entry points

| Owner | Entry point | Responsibility |
|---|---|---|
| PSCAD | `ext/LineCableModelsPSCADExt/PSCAD.jl` | Native formula methods, constitutive conversion, export/import, explicit station configuration, execution, readback and native recovery. |
| Gauntlet | `gauntlet/Gauntlet.jl` | Case ingestion, declared calculations, campaign scheduling, retained comparisons, reports and artifact locking. |
| LCM | `compute`, `validate`, `compare`, `details` | Existing computation, applicability and result grammar used by both entry points. |

`using LineCableModels` exposes `LineCableModels.PSCAD`. This module is statically
included; PSCAD has no Julia dependency whose loading would trigger an automatic
extension. Import and formulation construction perform no station or license
action. Gauntlet consumes the backend; production code does not load Gauntlet,
its catalogue, a test file or a local configuration file.

Use the [Gauntlet guide](../../gauntlet/README.md) for standalone REPL examples,
explicit station mapping, benchmark declarations, CLI commands and migration.
`test/gauntlet/` now contains tests rather than application implementations.

One `BenchmarkCalculation` retains its declared problem, formulation and options.
One benchmark/campaign path invokes the same public `compute` available in a
REPL. It does not reconstruct a backend formulation, replace frequency grids,
inject reductions or alter constitutive laws. Scalar, product/zip Gridspace and
existing UQ declarations use this path. The existing UQ moment settings continues
to accept one parameter point; this migration does not expand that settings.

## Physical and frequency boundaries

PSCAD keeps the shared formula tags and indexed self/mutual method grammar.
Unsupported cases, incompatible native switches, EHEM and unsupported material
laws fail at their existing applicability boundary. Each native default retains
its identity. Temperature correction uses the selected constitutive law.

`base_frequency` is a physical PSCAD formulation option, defaulting to 50 Hz.
It is supplied to conversion and to the native line's `Freq` field, with readback
verification. Requested/exported dielectric losses and the cap of ten are
retained. The native `enablf` control enables cable dielectric loss tangent;
passivity sampling is not a frequency-grid selector.

The verified unfitted phase-matrix scan accepts 101, 201, 501 and 1001 logarithmic
frequencies. The corresponding native menu has 100, 200, 500 and 1000 increments.
Roundoff-equivalent requests are accepted, and requested and native frequencies
are retained separately. The 0.1 Hz floor remains explicit. See the native
[frequency scan settings](https://www.pscad.com/webhelp-pscad-v5.1.0-ol/Master_Library_Models/Transmission_Lines_Cables/Distributed_Line_Models/fd_phase_options.htm)
and [detailed output description](https://www.pscad.com/webhelp-pscad-v5.1.0-ol/EMTDC/Transmission_Lines/Line_Constants_Program_Output.htm).

Arbitrary-grid execution remains an **adapter limitation**. Inspection of the
installed master library and its separate single-frequency and Bergeron settings
did not establish a route delivering the same full unfitted phase Z/Y matrices
with the required dielectric controls. This is a bounded feasibility finding,
not a claim that PSCAD cannot compute at other frequencies. Unsupported requests
fail without clamping, interpolation or substitution.

## Recovery and comparisons

Native input/completion records use schema 3. Reuse checks exported input,
requested frequencies, physical settings, installed solver identity, bundled
worker sources, native readback and raw-output hashes. An incomplete or altered
run is not accepted as a completed calculation.

Campaign numerical identity includes actual source and environment bytes and
explicit declaration/configuration sources. Fresh-process resume restores the
recorded dependencies and verified case builders before reading typed execution
declarations. Numerical source changes reject stale execution reuse.

Portable records retain plain matrix data, axes, port ordering, comparison
bindings, source bytes and copied native evidence. Locked bundles can be moved
and read without the original checkout or native work directory. Analysis-only
changes consume retained operands without launching solvers.

Reference direction defines the RMS denominator. LCM, PSCAD and FEM remain
independent calculations. Large cross-model differences do not fail a campaign.
Incompatible coordinates fail comparison. Numerically zero reference traces
retain absolute RMS and an unavailable relative RMS with its reason.

Two representative schema-1 archives were read without changing their hashes.
A historical typed PSCAD result was converted in the pinned checkpoint
environment using `gauntlet/migrate_typed.jl`, retaining its original file and
source evidence. The current reader recovered its 4×4×101 matrices without an
old-module alias or a solver rerun.

## Validation record

Validation logs and native evidence are retained locally under
`.linecablemodels/pscad-ext-verification/`; station configuration is not committed.

- Ordinary core, integration and non-graphical extension tests pass.
- Deterministic FEM numerical checks, package quality checks and the core-only
  dependency boundary pass.
- Cairo visual tests and GLMakie/WGLMakie activation checks pass.
- The documentation build and doctests pass.
- The full Gauntlet suite passes. Campaign checks cover definition authority, an independent backend, product/zip
  result spaces, copied native evidence, corruption rejection, exclusive locking,
  source-change rejection and archive packaging.
- Fresh-process UQ resume passes. An interrupted Monte Carlo candidate was
  recomputed with its declared seed and reproduced the retained numerical digest.
- Local worker tests use PythonCall with an automation double to exercise native
  field readback, solver-identity changes, terminal retention, raw-file collection,
  diagnostic retention and cleanup. These matrices are protocol fixtures, not
  electromagnetic reference values.
- Native scans completed for all four supported sample counts. A later standalone
  201-frequency run verified schema 3, 60 °C conductors, 60 Hz native base-frequency
  readback and retained raw-file evidence without loading Gauntlet.

The additional live standalone-versus-Gauntlet comparison is pending explicit
station/payload approval after automatic approval review rejected that action.
The local entry-point invariance checks pass; they do not establish that this
additional live comparison has run.

Production coverage passed at **15010/15796 lines (95.02%)**. No source-adjacent
`.cov` files remain; `lcov.info` and validation logs are retained. The source-amended 95%
threshold still includes all Julia production files under `src/` and `ext/`,
including the remote worker. Quality/Aqua, documentation and golden regeneration
are not used to satisfy that threshold. Collection cleans source-adjacent traces
after instrumented workers finish, including a failed gate, and retains LCOV.
