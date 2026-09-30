# Native FEM performance qualification

Execution plan: [fem-native-performance.md](../../../../docs/plans/fem-native-performance.md).
Historical screening: [results.md](results.md). Its thresholds do not authorize
withholding numerical choices. Harmonic consolidation is now in feature code;
physical quadrature and quadrilateral PML are explicit computation options.
The existing triangle and integration defaults remain unchanged.

## Feature comparison through the customary runner

```bash
# Serial, resumable; one mixed pair, ten frequencies; headless plots at the end.
julia --project=test test/manual/calculations/run_two_bare_wires_fem.jl --compare-native --save \
  > .linecablemodels/fem/native-feature-review/live.log 2>&1
tail -n 80 -F .linecablemodels/fem/native-feature-review/live.log
# Open the completed comparison with GLMakie; no computations.
julia -i --project=test test/manual/calculations/run_two_bare_wires_fem.jl --review-native
```

The comparison explicitly uses the retained reference's mesh factors 3/8,
ρ=0.1 Ω·m, radius 4.25 cm, heights −1/+1 m and Γ=0. It includes consolidated
triangles with 12/12 and 3/12 physical/PML points, plus consolidated PML
quadrangles with 3/9 points. The ordinary runner retains its 5/10 factors and
offers the same discretization options directly. `pml_quadrature=16` is also
supported and exercised by the managed/detached implementation test.

Results, raw component differences, exact settings and native/managed timings
are in `.linecablemodels/fem/native-feature-review`. The `onelab-quad9/study.pro`
bundle runs independently of Julia. First-use compilation is recorded separately.
No error threshold stops this comparison; the researcher evaluates the plots.

## Review the saved results

Open interactive GLMakie figures, with the same public plotting controls as
the two-bare-wires runner. This **does not execute the runner or invoke FEM**:

```bash
julia -i --project=test test/manual/calculations/fem_native_performance/plot_results.jl screen
```

Replace `screen` with `air`, `soil`, `mixed`, `air-finite`, `soil-finite`,
`mixed-finite`, `three`, `tube`, or `sector`. `air-low`, `screen-low` and
`tube-low` isolate 0.1 Hz, keeping the tiny conductance differences visible.
The default is `air`. `review_plots` retains the returned handles without
printing Makie internals. Raw observations use `clip=false` and zero cutoffs.

To export PNG/SVG figures without opening windows:

```bash
julia --project=test test/manual/calculations/fem_native_performance/plot_results.jl --save --all
```

Exports are in `.linecablemodels/fem/native-performance-20260930/review-plots`.
Existing exports are retained when resuming; `--replace` regenerates only these
named plot files. Panels containing both signs use symmetric native y limits
so tiny opposite-sign samples remain visible across large dynamic ranges.
Titles and legends state the actual
saved frequencies; cable fixtures have only two endpoints, and several pilot
variants have only one frequency. Markers denote computed samples; joining
lines do not supply missing intermediate results. Analytical series are shown
where analytical results were retained. None were computed for the cable-shape
fixtures in this campaign.

## Qualification record
This folder contains manual Julia qualification only. No acceptance policy,
retry, automatic refinement or solver search belongs to the production engine.

All evidence goes to `.linecablemodels/fem/native-performance-20260930`.
The active execution appends native output to one `live.log` there. Runs are
serial and completion files are reused. Original studies and runner settings
are not modified. Sources before feature edits are retained in `snapshot/`.

```bash
tail -n 60 -F .linecablemodels/fem/native-performance-20260930/live.log
```

`qualify.jl` compares separate changes on three retained meshes: aerial 0.1 Hz
at Gamma=0, mixed 0.1 Hz at Gamma=0.99 gamma_earth, and mixed 1 MHz at the same
fraction. `warm.jl` checks the combined assembly change and repeats matched
native executions twice, reversing execution order in the second repetition.
`mesh_io.jl` measures ASCII/binary MSH read/write costs separately.
`spectrum.jl` and `fixtures.jl` are the subsequent qualification stages, selected
only after inspecting the pilot results. They do not promote production code.
`separate_assembly.jl` isolates the two assembly changes on the failing screen,
reusing the saved baseline and stopping each candidate at its first failure.

Every nonzero raw R/X/G/B entry must change by at most 2%; exact zero and signs
are preserved. Pilot CSV files also contain the analytical values. No absolute
floor or clipping is used. Julia compilation is excluded from native GetDP
process timing; mesh preparation records its compilation separately. Initial
pilot times are screening measurements, not matched warmed performance claims.

## Assembly changes under qualification

For the exp(+j omega t) convention, define
`kappa_z = (sigma + j omega epsilon) s_x s_y`. The old finite-metal axial block
assembles conductivity and permittivity separately. The four combined entries
are `j omega kappa_z (a,a)`, `kappa_z (ur,a)`,
`j omega kappa_z (a,ur)` and `kappa_z (ur,ur)`. This changes neither the
finite-Gamma coupling nor the terminal constraints or voltage extraction.

For affine linear physical triangles with constant material coefficients,
the basis-function products have degree at most two. Native three-point Gauss
integration is exact for those polynomials in exact arithmetic. The candidate
uses GetDP's `Integration Criterion` to keep the original 12-point rule in the
nonpolynomial PML. It leaves the line integration rule unchanged.

## Initial quadrangle finding

The qualification copies pair triangles along their existing rectangular grid
diagonals. All node coordinates, physical elements, paths and perimeters stay
fixed. Gmsh's general Blossom recombiner was unsuitable for this preparation:
it can pair across grid lines on stretched cells. Native transfinite/recombine
construction would require a separate parity check before feature delivery.
No quad option has been added to production.

At aerial 0.1 Hz the 9-point quadrangle solve took 17.17 s versus 27.74 s for
triangles, with 275,018 versus 359,568 unknowns. However G21 moved from
−2.3266e−25 to −2.8749e−25 S/m (23.56%). The analytical value is −5.9537e−25 S/m.
All signs stayed unchanged; the 16-point rule still changed this component by
23.47%. Both fail the agreed preservation gate. Moving closer to the analytical
value does not establish convergence. Neither candidate proceeds to the
expanded fixture qualification.

References: [GetDP integration objects](https://getdp.info/doc/texinfo/getdp.html#Integration),
[Gmsh transfinite and recombined meshes](https://gmsh.info/doc/texinfo/gmsh.html#t6).
