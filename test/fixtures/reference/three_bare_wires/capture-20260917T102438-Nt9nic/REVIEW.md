# Step 0 baseline review

Captured on 2026-09-17, before any production changes in this execution round.
All five physical placements completed with both current unified and existing
quasi-TEM FEM. Each scan has nine frequencies, three distinct terminals, and
full ordered 3×3 Z/Y matrices. No fixture, frequency or matrix entry was dropped.

## Evidence

- 10 successful public `compute` calls; 90 frequency points across the two backends.
- 1,620 retained complex entries, saved with exact Float64 round-trip checks.
- 75 capture-file and 1,480 persistent FEM-artifact checksums verified after completion.
- 91 focused physical-fixture tests passed.
- Julia 1.12.7; Gmsh 4.15.0; GetDP 3.5.0 (PETSc 3.14.4).
- QuadGK controls remain the current resolved defaults: `rtol=1e-8`,
  `atol=0`, `maxevals=10000000`; the existing production acceptance logic was
  not changed for the capture.
- No warnings or failures in the captured public-compute logs. Native solver
  logs and inputs are retained in the referenced FEM directories.

`manifest.toml` identifies each run, source revision, dirty patch, environment,
terminal order, timing and file checksums. Each backend directory includes the
complete serialized physical problem, formulation record, computation options,
trace/details, timing/log, and Z/Y CSVs. `verification.toml` records the checksum
verification. Original execution sources are under `sources/`; the subsequently
extended comparison script is retained separately as `comparison-runner.jl`.

## Observed differences, not acceptance thresholds

For each quantity, this table selects the entry with the largest **absolute**
complex difference over all nine frequencies. The percentage is the relative
difference at that same entry, not a normwise error or the largest relative error.

| Placement | max abs ΔZ [Ω/m] | Relative at that entry | max abs ΔY [S/m] | Relative at that entry |
| --- | ---: | ---: | ---: | ---: |
| all air | 1.94850723 | 14.8994% | 1.26094674e-5 | 21.0683% |
| all earth | 0.01110244 | 0.137638% | 0.12638102 | 0.278505% |
| wire 1 in air | 1.87231519 | 3.88664% | 0.12586657 | 0.277378% |
| wire 2 in air | 1.70714875 | 3.32557% | 0.12511125 | 0.275711% |
| wire 3 in air | 2.13082152 | 4.44982% | 0.12560263 | 0.276792% |

All these Z maxima occur at 10 MHz. The all-air Y maximum occurs at 10 MHz;
the other Y maxima occur at 1 MHz on buried-wire diagonal entries.

The largest all-air entry-relative Z difference is **27.6661%**, at 10 MHz,
response wire 3 / source wire 1: absolute difference 1.93230389 Ω/m against
a FEM magnitude of 6.98437655 Ω/m. This is not a tiny-denominator effect.
Below 10 MHz the all-air complex entry-relative differences stay below 0.781%
for Z and 1.402% for Y on the sampled grid.

Some mixed-layout relative discrepancies are much larger but involve small
cross-interface mutual entries. For example, wire-1-in-air Y at 1 MHz,
response wire 1 / source wire 3, differs by 8.21753790e-6 S/m against a FEM
magnitude of 6.01129668e-9 S/m. This is retained explicitly, not hidden in a
matrix norm or portrayed as the error of the large buried-wire diagonals.

`comparison.toml` records both kinds of extrema, coordinates and denominator
magnitudes. Each placement's `Z-differences.csv` and `Y-differences.csv` retain
absolute and relative complex/component differences for every entry. Relative
errors are marked NaN at arithmetic-scale zero denominators; their absolute
errors remain present. This reporting convention is not an acceptance floor.

## Disposition

The computations completed, but this is **not a demonstrated matching
baseline**. No final comparison tolerances have been approved. The cause of
the discrepancies is not established by these measurements alone.

Per §0.3 of the locked plan, production cleanup awaits the user's baseline
review/direction. No production refactor, tolerance relaxation, FEM alteration
or commit was performed in this execution round. Preserve these payloads;
do not regenerate them to make a subsequent refactor pass.
