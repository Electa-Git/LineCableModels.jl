# Corner coarsening: bounded qualification stopped

This records the original strict-preservation decision. A later user request
explicitly accepted its near-DC magnitude tradeoff for an optional manual-test
feature; see [feature-results.md](feature-results.md). The numerical evidence
and original decision below are retained unchanged.

Completed 2026-09-30. Both prescribed candidates failed the agreed 2% component
preservation limit. The three pilots were completed with both candidates:
nine native solves, eighteen source columns. No production implementation or
default change follows this result. The broader fixture sweep and repeated
timing campaign were not started because neither candidate passed accuracy.

Evidence: `.linecablemodels/fem/pml-corner-mesh/`. The live log, native logs,
frozen sources/meshes, raw matrices, signed analytical comparisons,
`pilot-assessment.csv` and `pilot-selection.toml` are retained. All numerical
work was serial and used Julia/Gmsh/GetDP. No Python or subagents were used.

## Measured results

The metric prescription and two fixed settings are described in [README.md](README.md).
Finite Gamma means 0.99 times the complex earth propagation constant. The
reference and candidates use the same corrected quasi-fw equations.

| Pilot | Metric scale | Corner triangles, before → after | Maximum component change | New signs | Native seconds, before → after |
| --- | ---: | ---: | ---: | ---: | ---: |
| Air, 0.1 Hz, Gamma=0 | 1.00 | 103824 → 87692 | 538595% | 0 | 27.80 → 30.18 |
| Air, 0.1 Hz, Gamma=0 | 1.35 | 103824 → 48356 | 314359% | 0 | 27.80 → 20.90 |
| Mixed, 0.1 Hz, finite Gamma | 1.00 | 47112 → 24694 | 0.000001538% | 0 | 33.03 → 29.42 |
| Mixed, 0.1 Hz, finite Gamma | 1.35 | 47112 → 13942 | 0.000002500% | 0 | 33.03 → 26.69 |
| Mixed, 1 MHz, finite Gamma | 1.00 | 44104 → 23506 | 0.0004768% | 0 | 33.68 → 30.80 |
| Mixed, 1 MHz, finite Gamma | 1.35 | 44104 → 13404 | 0.002912% | 0 | 33.68 → 28.05 |

Every entry of R/X/G/B is compared individually without clipping or relative
error floors. The failure is the very small aerial conductance. For G[1,2]:

| Aerial 0.1 Hz result | Conductance, S/m |
| --- | ---: |
| Retained FEM | -2.3422506890550553e-25 |
| Metric scale 1.00 | -1.2534291686657373e-21 |
| Metric scale 1.35 | -7.316754025275699e-22 |
| Analytical | -5.95372593928867e-25 |

Both candidates retain the negative sign but substantially worsen magnitude
agreement with both the retained FEM and the analytical result. This is an
error in raw solver output, not plotting or postprocessing clipping. The
other components of the aerial result change much less; agreement in dominant
components does not establish accuracy of its very weak conductance.

The coarser aerial mesh has 125424 total triangles versus 180892 (30.7% fewer)
and 248632 DOFs versus 359568. Its one native solve is 24.8% faster, but building
the replacement corners takes another 3.99 s. The conservative candidate's
first corner construction takes 6.54 s, including 0.80 s of Julia compilation;
its native solve is slower despite fewer elements. Other corner-construction
times range from 1.03 to 1.54 s without recorded compilation.

These are individual native timings plus prototype construction timings, not
repeated warmed end-to-end performance measurements. The coarser candidate's
median native-time ratio over the three pilots is 0.808, before adding meshing.
The 20% warmed total-time objective is therefore not established. No further
timing solves were spent on scientifically rejected candidates.

## Isolation and mesh validity

For all six candidate meshes, the physical-region coordinates/connectivity,
noncorner PML coordinates/connectivity and native voltage-path coordinates
retain exactly matching hashes. Corner perimeter edges match the frozen
reference node-for-node, with no hanging nodes. Post-run checks confirm positive
triangle areas, at most two triangles per interior edge, complete rectangular
area coverage and exact perimeter connectivity in all four corners. Each
candidate has an `integrity.toml` record.

The first preparation attempt merged the reference mesh into independently
exported CAD. It changed entity classification and was interrupted without
solving. The successful method meshes the four corner polygons in separate
native Gmsh models and inserts only their interior triangles into the frozen
mesh. The unsuccessful preparation remains in the append-only log.

Gmsh emits `No vector or scalar value found in PostView field` during scalar
size queries on the tensor-only background. Its `PostViewField` has a separate
tensor evaluation operator consumed by BAMG; there were no missing-tensor
warnings in successful construction. See the
[Gmsh source](https://raw.githubusercontent.com/live-clones/gmsh/master/src/mesh/Field.cpp).
The tensor metric was used; this warning is distinct from the rejected initial
CAD/mesh import. The conforming meshes and numerical comparisons above are
retained, rather than treating the warning as a solver failure.

The copied `.geo` files describe the structured reference geometry. Candidate
experiments use their explicitly supplied `study.msh`; they are not released
ONELAB remeshing options. Detached feature support was conditional on a
surviving candidate, and was not implemented.

## Decision and scope

Stop this prescribed corner-metric approach at qualification. Keep the working
structured PML in production and the manual runner. Production source hashes
in `production-before.sha256` remain unchanged. All new code is confined to
this manual qualification directory.

The finite-Gamma pilots show that substantial corner coarsening can preserve
those particular results. They do not justify a Gamma-dependent automatic
switch, a general corner option, or a claim that other corner prescriptions
cannot succeed. This experiment identifies the low-frequency aerial Gamma=0
conductance as the blocking observable for these two prescriptions.

The subsequent [fixed-node diagonal and soil diagnostics](diagonals-and-soil-results.md)
passed: reversing corner diagonals changes at most 0.0349%, and localizing
soil in a fresh export changes at most 0.832%, with no new signs. The displayed
full-width soil band came from an older detached export. These diagnostics
do not qualify the rejected corner-coarsening candidates.
