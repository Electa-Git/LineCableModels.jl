# Retired anisotropic PML corners: historical results

Removed from production on 2026-09-30 at the user's request. The modest measured
speed gain and degraded conductances did not justify retaining the feature.
The controls and usage below document the retired implementation and are no
longer accepted. Julia and detached ONELAB now use structured corners only;
the two-stage mesher, tensor fields and corner-specific UI controls were removed.
Computed results remain preserved under `.linecablemodels/fem/pml-corner-mesh/`.

Removal verification: 3,016 focused assertions passed; one existing detached
versus managed PML-count assertion failed identically before and after removal
(4,896 versus 4,864 triangles). Its tolerance was not changed. The current manual
runner's representative mixed case passed all ten ordered frequency meshes,
and native detached `BuildMesh` passed at 129.1549665014884 Hz. These were mesh
checks, not a new numerical solve campaign. Logs are `structured-removal-tests.log`,
`structured-removal-baseline.log` and `structured-removal-runner.log` in that
results directory. Mesh cache identity advanced to partition version 18.

Implemented 2026-09-30 following the user's explicit decision to accept the
near-DC magnitude tradeoff for a manual engineering comparison. The fastest
of the two earlier sign-retaining candidates was BAMG with metric scale 1.35.
The earlier failed 2% qualification remains recorded in [results.md](results.md);
this implementation does not retroactively pass that criterion.

## Historical use

```julia
pml_resolution = (interpolation_cells=72, coefficient_change=0.12),
pml_corner_mesh = :anisotropic,
pml_corner_scale = 1.35,
```

These controls were selected in `../run_two_bare_wires_fem.jl`, including
its single-case detached export, during the feature comparison. They have
since been removed from that runner and the public options.

Both public `compute` and `export_data(:onelab, ...; mesh_options=...)` support
the option. In the exported ONELAB project, use **Run**, optionally with
**Mesh only**. The ordinary Mesh/2D shortcut does not invoke the two-stage
construction. **PML corner metric scale** is editable in ONELAB; changing the
structured/anisotropic geometry choice requires another export. Detached
execution requires native Gmsh/GetDP, with no Julia or Python interpreter.

## Construction

The four Cartesian corner interiors become polygons with the same graded
perimeter segments as the adjacent strips. Physical-region meshing, conductor
grading, strip prescriptions, PML stretch, equations and voltage paths retain
their existing prescriptions. Native Gmsh first meshes the ordinary regions,
then fills the four empty corners using BAMG and a tensor background field.
This is a fixed construction, with no scientific acceptance, retries, clipping
or adaptive solution refinement.

The metric uses the existing frequency/material/prescribed-Gamma modal density
and each orthogonal direction's finite-layer amplitude envelope. The static
mode is included. The 1/8 amplitude floor and 2/3 exponent bound the additional
cross-direction length relaxation by four. Scale 1.35 multiplies requested
interior lengths. This is an interpolation heuristic, not an error bound on
terminal quantities. See [Gmsh t17](https://gmsh.info/doc/texinfo/gmsh.html#t17).

Compact one-dimensional samples are retained in mesh plans and cache identity;
the tensor view is expanded only when constructing or exporting a mesh. A
native `MinAniso` field wraps the tensor `PostView` so Gmsh's scalar queries
also receive a size, without the earlier prototype's missing-scalar warnings.
The feature consequently need not reproduce the earlier isolated-corner
prototype's exact triangulation or values. Independent fresh unstructured
meshes also need not be identical between the Julia and GEO routes.

## Direct feature measurements

All three pilots were solved through public Julia-managed execution and through
independently meshed native detached execution: six solves, twelve source
columns. No new R/X/G/B signs appeared relative to the retained reference in
either route. Soil resistivity is 0.1 ohm m; wire radius is 4.25 cm. Finite Gamma
means 0.99 times the complex earth propagation constant.

| Pilot | Retained corner triangles | Detached feature corner triangles | Feature total triangles | Maximum detached change from retained FEM |
| --- | ---: | ---: | ---: | ---: |
| Air, 0.1 Hz, Gamma=0 | 103824 | 51326 | 128398 | Very weak G: about 10006 times the reference magnitude |
| Mixed, 0.1 Hz, finite Gamma | 47112 | 14740 | 147453 | 0.00371% |
| Mixed, 1 MHz, finite Gamma | 44104 | 14144 | 160872 | 0.00170% |

The aerial reference has 180892 total triangles: the detached feature removes
29.0%. Its G12 is -2.3283e-21 S/m; the initial managed feature run gives
-2.9295e-21 S/m. The retained reference is -2.3423e-25 S/m and the analytical
value is -5.9537e-25 S/m. The negative sign is retained, but magnitude accuracy
of this very small conductance is sacrificed. These discrepancies have not
been established to be pure floating-point noise. Nothing is clipped.

Managed aerial R/X/B changes are at most 0.0279%. Managed mixed low-frequency
changes are at most 0.00371%. For the mixed high-frequency case, managed X
changes by up to 2.065% against the saved reference; its independent native
and managed solutions differ by up to 2.064%. Thus the two routes are not
claimed to have componentwise 2% parity. Their full signed values, including
the analytical comparisons, are in each pilot's `components.csv`.

Native detached timings for the feature, measured separately from export:

| Pilot | Meshing [s] | Assembly [s] | Factorization/solve [s] | Native solver total [s] | Native peak RSS [GiB] |
| --- | ---: | ---: | ---: | ---: | ---: |
| Air, 0.1 Hz | 3.13 | 13.67 | 5.73 | 22.75 | 1.28 |
| Mixed, 0.1 Hz | 2.01 | 18.68 | 5.57 | 27.28 | 1.45 |
| Mixed, 1 MHz | 2.28 | 18.06 | 6.62 | 28.53 | 1.58 |

These are individual native measurements, not repeated workload medians.
Assembly and factorization/solve come from native GetDP event timings.

## Matched public runtime

Two warmed, serial public `compute` repeats per method on aerial 0.1 Hz include
model preparation, meshing, solving and extraction. Each run forces a fresh
mesh, with one frequency worker and one solver thread. Methods alternate.

| Corner method | Warm repeat 1 [s] | Warm repeat 2 [s] | Median [s] |
| --- | ---: | ---: | ---: |
| Structured | 28.5521 | 28.6613 | 28.6067 |
| Anisotropic, scale 1.35 | 25.8879 | 26.4646 | 26.1762 |

The measured end-to-end reduction is **8.50%**. It is useful but does not meet
the original 20% target. Fewer triangles do not translate directly into the
same elapsed reduction: native BAMG construction costs more than constructing
the rectangular corner grid. This timing does not establish a speedup for an
entire multi-frequency/manual parameter sweep.

The first pair was retained as warm-up evidence and excluded from these
medians. Structured first use took 32.5058 s, including 4.0562 s of Julia
recompilation. The first managed feature verification took 44.3777 s, including
18.7048 s of Julia compilation. All four reported warmed runs recorded zero
compilation/recompilation time. Raw records are in `feature-execution/benchmark`.

## Verification and evidence

The focused native construction test checks conforming shared edges, correct
material groups, compact metric storage, distinct cache identities and
unchanged physical element counts between the two detached meshing stages.
It also exercises rebuilding the managed exterior for a different frequency.
The initial corner test passed 769 assertions, but did not cover the later
manual-runner failure described below. An initial public two-frequency aerial
computation at 0.1 and
1 Hz completed in 51.09 s, including 2.15 s of recompilation. All eight raw
conductance entries stayed negative. Its run directory and values are retained
in `feature-execution/final-frequency-reuse.toml`; this verifies the actual
managed frequency sweep as well as mesh-only reconstruction.
An initial native macro synchronized GEO after the first pass and remeshed
physical regions; it was corrected by completing GEO definitions before that
pass and retaining surface visibility restrictions. The failed initial runs
remain under `feature-execution-initial`; they are not counted as final results.

The broader focused run recorded 1076 passing assertions and one failure.
The failure is the existing detached frequency-scan test's byte-for-byte
comparison of independently remeshed numerical tables at 10 kHz. Differences
are about 1.5e-8 relative in the displayed resistance values. Repeating that
test with the original structured `Mesh 2` command also gives 33 passes and
the same strict text-equality failure. Its tolerances were not loosened and
its tables were not altered. Native scan ordering, publication and mesh-only
checks pass; this is not reported as an entirely green test suite.

Evidence root: `.linecablemodels/fem/pml-corner-mesh/feature-execution/`.
`verify_feature.jl` and `benchmark_feature.jl` are resumable Julia manual
harnesses. The existing full interactive manual sweep was not executed.
No new all-frequency/all-geometry sign guarantee is inferred from three pilots.

```bash
tail -n 80 -F .linecablemodels/fem/pml-corner-mesh/live.log
```

The numerical changes are opt-in. Existing equation sources and native GetDP
measurement integration were not modified for this feature.

## Correction after the actual manual runner failed

The user's first resistivity case failed at its fifth frequency, 129.1549665 Hz:
the saved mesh had no elements on air surface 303. The run record confirms
zero GetDP invocations. A fresh Julia process reproduced the exact failure;
the accompanying Revise warning was not its cause.

The corner pass restricted visible surfaces but left all curves visible.
An empty discrete scratch curve left by native conductor boundary-layer
meshing made Gmsh regard its 1D mesh as unfinished. Before building the
corners, Gmsh restarted 1D generation and erased the physical 2D mesh.
The earlier attempt to remove scratch curves from GEO did not prevent the
failure through the complete frequency sequence.

The correction limits the second pass to the four corner surfaces **and their
already-meshed perimeter curves**, then restores visibility. Julia and native
GEO implement the same sequence. The ineffective scratch-curve deletion was
removed, and the mesh cache version increased to 17. No physical settings,
mesh densities, equations, voltage paths or numerical tolerances changed.

The regression now deliberately supplies an empty discrete curve and verifies
that physical groups and shared edges survive. All 770 corner assertions and
35 conductor-grading preservation assertions pass. The dedicated manual check
also passes every mesh of the actual runner: four resistivities and three
radii, each at ten frequencies. Counts and timings are retained in
`.linecablemodels/fem/pml-corner-mesh/runner-mesh-check.csv`.

The formerly failing 129.1549665 Hz case also completed through both public
Julia compute and the detached native bundle. Managed compute took 34.07 s,
including 15.72 s of compilation; native solver wall time was about 16.2 s.
The completion records are in `runner-repair-solve.toml`, and the detached
bundle is `runner-repair-native/study.pro`, under the same evidence root.

`verify_runner_meshes.jl` reads the runner's setup before its export/compute
calls, omits interactive imports, and checks all 70 meshes in their actual
frequency order. It runs no solver grid and opens no plot windows. It appends
to the existing live log. This mesh-construction check is not a new scientific
qualification of all 70 solutions.
