# Both-media localization: qualification outcome

Execution on 2026-09-30 used the frozen FEM reference and the authorized 2%
limit for each individual R/X/G/B entry, with no new component sign reversals.
The initial unchanged-equation comparison below exposed a finite-Gamma
metal-drive cancellation defect. The subsequent investigation and exact
substitution are recorded in [finite-gamma-diagnosis.md](finite-gamma-diagnosis.md).
All six finite-Gamma pairs then passed with the correction on both meshes,
and the Gamma=0 identity control was bitwise unchanged. Production now applies
both-media localization and the qualified exact substitution; public execution
and regression verification are complete. Original FEM references remain untouched.

## Results

Compact footprints passed all nine Gamma=0 cases and the mixed 0.1 Hz finite-
Gamma case. The mixed 1 MHz case at Gamma=0.99 gamma_earth failed. The wider
footprint was tested on that failing case first and also failed. The four
remaining finite-Gamma aerial/buried points were not run, following the fail-fast
plan. Ten new qualification solves ran; the initial two compact results and all
baseline results were reused. All twelve comparisons retained component signs,
including the two failures.

| Placement, Gamma=0 | Frequency | Maximum component change | Nodes before → after |
| --- | ---: | ---: | ---: |
| Air | 0.1 Hz | 1.62393% | 100844 → 91044 |
| Air | 1 kHz | 0.30559% | 76117 → 66425 |
| Air | 1 MHz | 0.59236% | 39454 → 36215 |
| Soil | 0.1 Hz | 0.08456% | 101869 → 92098 |
| Soil | 1 kHz | 0.14217% | 76674 → 66954 |
| Soil | 1 MHz | 0.78667% | 40488 → 37356 |
| Mixed | 0.1 Hz | 0.73969% | 101290 → 91518 |
| Mixed | 1 kHz | 0.36202% | 76402 → 66694 |
| Mixed | 1 MHz | 0.40172% | 40234 → 37097 |

At aerial 0.1 Hz, air/soil triangles decrease from 16235/14535 to 6451/4719,
removing 63.7% of physical medium triangles. The PML retains 169100 triangles.
All completed candidate meshes preserve PML, conductor-contour and voltage-path
coordinate hashes exactly. The mesh images show the broad full-width bands
replaced by refinement concentrated around the cables.

Screen, tube and sector mesh checks passed the same geometry/prescription
contracts. Nodes decreased 181632→171823, 67746→58017 and 147518→137767.
Conductor controls are unchanged; unstructured metal interiors need not be
bitwise identical. Their triangle counts were 163506→163502, 8756→8756 and
137188→137192. These are mesh contract checks, not new cable accuracy claims.

## Blocking mixed finite-Gamma case

At 1 MHz, Gamma=0.99 gamma_earth, baseline and both candidates have 66248
physical air/soil triangles. Total nodes are 95963, 95964 and 95965 respectively.
The largest changes are in the source-2 mutual column:

| Component | Reference FEM | Compact | Wider | Compact / wider relative change |
| --- | ---: | ---: | ---: | ---: |
| R12 [ohm/m] | 0.0026022901 | 0.0033504476 | 0.0034726188 | 28.75% / 33.44% |
| X12 [ohm/m] | -0.0024531667 | -0.0038637183 | -0.0040903292 | 57.50% / 66.74% |
| G12 [S/m] | -1.4507844e-7 | -1.7941002e-7 | -1.8505757e-7 | 23.66% / 27.56% |
| B12 [S/m] | 1.9231773e-7 | 2.9485328e-7 | 3.1134797e-7 | 53.32% / 61.89% |

This is not a marginal failure of a 2% limit. No tolerance was relaxed and no
unchanged-mesh identity control was substituted for these failed comparisons.
Earlier independent remeshing of the untouched reference prescription changed
X12 by 56.1%.

`check_constant_bounds.jl` confirmed that this case's wave-size prescriptions
are unchanged throughout the physical domain. Air has equal threshold minimum
and maximum (0.0426286 m). The physical rectangle is x in [-4.5,5.5] m and z in
[-5,5] m. Over the entire soil rectangle, distance to a compact interface
footprint is bounded by 6.69846 m, below the 6.76934 m threshold onset. Thus the
soil wave-size target remains 0.299164 m everywhere for both prescriptions;
the wider footprint can only decrease that distance. This is an analytic bound
using the saved native geometry and field values, not sampled error acceptance.

The failed comparison therefore exposes existing remeshing sensitivity without
changing those size targets. It does not justify either relaxing the 2% gate or
claiming that localization requires the original full-width band.

## Cost and artifacts

The fresh matched mixed 0.1 Hz repeat used the frozen reference and compact
meshes, with hashes checked. It took 32.43→32.60 s; GetDP CPU was
31.5566→31.7507 s. Peak RSS decreased 2004652→1816200 KiB (9.4%), and DOFs
400546→361458 (9.8%). Assembly was 22.10→22.83 s; factorization/solve was
7.45→7.22 s. The original pair was 30.21→28.29 s. Mesh and memory savings are
established; a dependable runtime improvement is not. Julia startup/compilation
is excluded from these native process timings. Other runs showed substantial
timing variation and should not be compared to old baselines as speedup proof.

Evidence root: `.linecablemodels/fem/local-interface-mesh/`.

- `both-media-selection.toml`: failed general qualification; no chosen candidate.
- `both-media-assessment.csv`: FEM preservation and analytical errors separately.
- Case `two-percent-comparison*.toml` and `components*.csv`: all raw entries.
- `mixed-f1.0e6-gamma0.99/constant-target-bounds.toml`: unchanged-field proof.
- `shapes/*/preservation-both.toml`: cable constraint checks.
- `warmed-repeat-localized/*/solve.toml`: fresh matched timing and memory.
- `both-media-figures/air-f0.1-gamma0.0/localized/`: before/after PNG/SVG meshes.

The single log is `live.log` under that root. `qualify_both.jl` resumes completed
candidate meshes/solves without rebuilding or solving frozen baselines. The
initial failed comparisons remain evidence of the defect, not successful
preservation checks. `both-media-corrected-selection.toml` records the subsequent
15-case qualification with the metal substitution held fixed on both meshes.
