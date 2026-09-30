# Conductor mesh delivery assessment — 2026-09-29

The conductor mesh feature is implemented and its bounded delivery assessment
is complete, superseding pending-batch statements in the historical execution
notes. This closes the conductor mesh scope, not the separate bare-wire shunt
sign investigation. Numerical qualification stays in this manual harness; no
scientific acceptance, retry or refinement policy belongs in the production
engine.

Quasi-fw is the user's primary target. Its coupled solves and comparisons are
retained separately; quasi-TEM and the local electrostatic controls isolate
discretization effects but do not substitute for quasi-fw evidence. In the final
quasi-fw results, screen Y differs from the qualified prototype by 0.003099%,
and sector Y and loop R change by 0.067236% and 0.069526% under refinement.
Small complex-Y differences do not establish the sign of a much smaller real
part; the raw quasi-fw G discrepancies remain open.

## Completed evidence

The public conductor controls, round/annular boundary layers, sector strips,
frequency reuse, cache identity and editable native export are implemented.
The retained native measurement contract is unchanged: grouped terminal
quantities and conforming GetDP field-edge integrals, with no external
measurement preparation or contour averaging.

The bare-wire selection completed 138 cases / 374 source columns. Its maximum
self-R difference from the matched analytical reference is 0.9292%. Its raw G
sign disagreements remain open; conductor resistance convergence is not a
certificate of shunt accuracy.

All 36 cable cases / 108 source columns now completed: 12 public scans and six
same-mesh detached comparisons. Primitive Z/P are identical in all six parity
comparisons. The following maxima include both formulations:

| Fixture | Self-R change under the prescribed 1 MHz refinement | Largest complex Y-entry change | Largest 0.1 Hz DC loop deviation |
|---|---:|---:|---:|
| 49-wire screen and thin foil | 0.062015% | 11.05621% | 0.070258% |
| Tubular Pb sheath | 0.00078242% | 0.24354% | 0.071447% |
| Three sectors and neutral | 0.086865% | 3.32555% | 0.042288% |

Refinement differences are not independent reference errors. In particular,
the prescribed refinement changed conductor geometry and normal grading;
it did not independently qualify the surrounding dielectric discretization.
The 12 public scans took 2696.6 seconds in total (four frequencies in each
normal scan, one in each refined scan). This includes first-use Julia work
and is not a warmed performance comparison or the total campaign cost.

The cable plot failure was a manual caller defect: `ObservedResult` requires a
complete primary representation, so `(R,G,B)` omitted `X`. It now retains
`(R,X,G,B)` and selects `(R,G,B)` for display, matching the other manual
qualification plotters. The production API was not weakened. Postprocessing
alone produced all 18 cable figures as SVG and PNG; the screen B and sector R
pages were visually inspected. No completed numerical case was rerun.

Evidence root: `.linecablemodels/fem/conductor-mesh-qualification/`.
The cable tables are `cable-fixtures/{refinement,dc-loops,costs}.csv`, and plots
are `cable-fixtures/<design>/plots/<physics>/`. The captured numerical-source
manifest still matches the live tree after the plot recovery.

## Screen shunt root cause established locally

The suspect B entry has virtually the same value and refinement change in
quasi-fw and quasi-TEM. To isolate it, extract only the PE inside the foil from
each retained mesh. Set the foil to zero potential, drive core/screen with
`(1,0)`, `(0,1)` and `(1,1)` V, and let GetDP integrate `grad(v)^2`. These
three energies recover the two-terminal capacitance matrix referenced to the
foil, hence `C23 = -C22-C21`. No earth, PML, skin effect or line integral is
present in this control. See [the caller](diagnose_screen_capacitance.py) and
[the native formulation](screen_capacitance.pro).

| Internal PE mesh | B23 at 1 MHz [S/m] |
|---|---:|
| Saved normal mesh | -0.0371340199576 |
| Saved normal, each triangle subdivided once | -0.0309838120293 |
| Saved normal, subdivided twice | -0.0286993406560 |
| Saved refined mesh | -0.0334371406833 |
| Saved refined, subdivided twice | -0.0279813683365 |

The unsubdivided local controls reproduce the full quasi-TEM B23 within
approximately 3.1e-10 relatively. Subdivision keeps the loaded polygon
boundaries fixed. Thus field discretization alone accounts for a substantial
part of the discrepancy; changing the physical exterior cannot cure it.
This control does not diagnose the separate bare-wire G sign disagreements.

Two concrete sizing interactions explain the deficient insulation mesh:

1. `_configure_conductor_mesh!` sets circle edge counts from angular geometry
   tolerance alone. The normal foil consequently has roughly 1.96 mm edges,
   despite its existing region size target of 0.225 mm. That target is applied
   to the metal bulk, but transfinite edge counts override ordinary edge sizing.
2. `Mesh.MeshSizeExtendFromBoundary=0` avoids spreading fine sizes throughout
   the exterior. Conductor normal layers correctly exclude the dielectric,
   but no replacement field extends the local conductor boundary resolution
   into the internal insulation. Its bulk target instead comes from the outer
   jacket (3.675 mm at these controls), across a 0.4 mm minimum wire/foil gap.

## Qualified local candidate and completed coupled confirmation

A native prototype caps the foil edge length by its existing local target and
adds an `Extend` field from dielectric boundary curves, restricted to the
passive cable surfaces. The exterior and metal normal-layer prescriptions are
unchanged. This uses ordinary native
[Gmsh Extend/Restrict fields](https://gmsh.info/doc/texinfo/#Gmsh-mesh-size-fields),
not adaptive acceptance, field-error estimation or an external mesher.

| Local control | B23 [S/m] | Full mesh nodes |
|---|---:|---:|
| Foil edge cap alone | -0.0312799029417 | 303216 |
| Dielectric extension alone | -0.0319917583135 | 320957 |
| Both, 0.225 mm foil cap | -0.0274549055952 | 339484 |
| Both, 0.1125 mm foil cap | -0.0273915393832 | 353231 |

The two combined controls differ by 0.23133% in B23. Independent subdivision
of the finer candidate's fixed polygon gives -0.0273299456360 S/m, a further
0.22537% change. This is convergence evidence, not an exact reference. The
normal combined candidate adds 15.6576% nodes over the saved normal full mesh
(293525). Its first-order triangles have positive Jacobians and every retained
voltage-path segment is still an electric-field mesh edge. Local three-drive
solves took about 4 seconds per combined mesh. Full coupled costs are below.

Reproducible native inputs, mesh snapshots, scripts and CSVs live in
`screen-capacitance-diagnosis/`. All four prescribed coupled controls completed
before implementation: two meshes, both formulations, twelve source columns.
Their matrices and `coupled-summary.csv` are retained under that directory.

| Coupled control | Quasi-fw B23 [S/m] | Quasi-TEM B23 [S/m] | GetDP wall seconds, fw / TEM | Peak reported memory, fw / TEM [MB] |
|---|---:|---:|---:|---:|
| 0.225 mm cap + Extend | -0.0274549082397 | -0.0274549056118 | 252.024 / 97.8171 | 5270.34 / 2440.95 |
| 0.1125 mm cap + Extend | -0.0273915420199 | -0.0273915393997 | 257.515 / 102.447 | 5392.89 / 2519.64 |

Full B23 agrees with the local electrostatic values within 9.64e-8 relatively
(quasi-fw) and 6.05e-10 (quasi-TEM). Between the two meshes, the largest complex
Y-entry change is 0.295955%, primitive P change 0.031573%, complex Z change
0.000756%, and self-R change 0.000423%.

Self-R contains the common earth-return contribution, so it can conceal larger
local resistance errors. Assessing every net-zero two-terminal loop gives a
maximum R refinement change of **0.112875%** for the corrected pair of meshes.
For comparison, the original normal/refined screen controls changed loop R by
up to **9.42621%**, despite small self-R changes. The corrected dielectric mesh
therefore matters for both shunt coefficients and local magnetic coupling.
These differences are convergence evidence, not independent exact errors.
The old native runs also exported field maps, so their wall times are not a
matched performance baseline for the map-free controls above.

## Implementation and completed delivery check

The qualified correction is implemented in existing owners:

- `geometry.jl` supplies the round arc's length alongside its angular fraction.
- `mesh.jl` respects both the angular bound and the existing local edge-size
  target, and restricts a native boundary-size Extend field to each cable's
  passive surfaces. Mesh fingerprint version is 12.
- `export.jl` writes matching edge constraints and scales the Extend distance
  consistently with native `MeshScale` edits.

The existing user controls are sufficient; no option, solver equation, voltage
extraction method, averaging step, scientific acceptance gate or fallback was
added. Bare-wire cases have no passive cable surface for this Extend field.
At the tested angular resolution their existing round edge counts remain 96.

All **125 focused conductor-mesh checks pass** in 100.36 seconds, including
high/low/high frequency reuse, native material/frequency edits, screen path
conformity, local foil edge limits, native MeshScale edits, sector partitions and
both faces of an annular wall. This is integration evidence, not a scientific
accuracy guarantee.

Retained current/loss fields were rendered and visually inspected for all three
fixtures under source 1, quasi-fw, 1 MHz. They show surface concentration in
solid cores, thin foil and thick tubes, and current crowding around sector
boundaries and neutral wires. Those images use the original graded-conductor
meshes, before this interface correction. They do not certify integrated power
balance or the remaining shunt quantities. See
`cable-fixtures/retained-current-loss-maps.png` and its rendering script; colours
span seven decades, with lower values omitted from the display only.

The bounded delivery selection checked the implemented public route at
1 MHz: screen and tube normal meshes, sector normal and refined meshes, both
formulations. This is **8 cases / 26 source columns**, not another spectrum or
shape campaign. Fresh native exports of all four inputs are structurally checked
before any solve. Completed original runs stay untouched; new source identity,
meshes, raw Z/P/Y and costs go to `interface-delivery/`.

All eight cases completed. Saved raw matrices were assessed without another
solve; `interface-delivery/assessment.csv` contains both formulations separately.
The following maxima include both formulations. Each percentage compares an
entry against the corresponding entry of the second named result; loop R uses
`real(Zii + Zjj - Zij - Zji)` over every terminal pair.

| Delivered comparison at 1 MHz | Maximum complex Z difference | Maximum complex Y difference | Maximum self-R difference | Maximum loop-R difference |
|---|---:|---:|---:|---:|
| Screen normal vs qualified cap + Extend prototype | 0.00000434% | 0.003100% | 0.00000487% | 0.001283% |
| Tube normal vs original normal | 0.105428% | 0.686620% | 0.000208% | 0.053297% |
| Sector normal vs delivered refined | 0.003742% | 0.067236% | 0.014691% | 0.069526% |

The screen implementation reproduces the qualified native recipe. Tube numbers
measure the change from the previous implementation, not an independent error
or a new refinement estimate. Sector Y refinement decreases from 3.32555% to
0.067236%; its loop-R comparison also stays below 0.070%. The interface sizing
correction therefore addresses the observed cable shunt discretization problem
without introducing a new extraction method or production validation policy.

Public normal-case wall times were approximately 158/63 seconds for the screen,
104/40 seconds for the tube and 131/58 seconds for the sector (quasi-fw/TEM).
The refined sector took 198/102 seconds. These include Julia work and are
recorded costs on this host, not matched warmed speedups. Source identities,
fresh native mesh checks and each run's cost table are retained with the results.

The public controls, Julia compute route, detached export, manual runner wiring,
documentation and bounded fixture evidence are delivered. No further historical
grid is required for this feature. Bare-wire G sign disagreements remain a
separately documented limitation; this dielectric correction does not resolve
them. The earlier cable plots and current/loss maps retain their original-mesh
meaning; they have not been relabelled as plots of these final eight cases.
