# Local interface mesh grading

**Current implementation:** physical-medium wave fields use finite cable
footprints in both air and soil. `interface_refinement_factor` widens those
footprints in Julia and detached ONELAB. Conductor grading, Cartesian geometry,
PML nodes and native voltage paths retain their existing prescriptions.

The finite-Gamma failure was traced to cancellation between the metal axial
drive and Gamma-squared scalar potential. An exact change of unknown removes
that cancellation before assembly, with no additional unknown or solve. See
[the diagnosis and controls](finite-gamma-diagnosis.md) and
[both-media qualification](both-media-results.md).

All nine Gamma=0 and six finite-Gamma mesh comparisons passed the authorized
2% component limit without new signs. The Gamma=0 substitution control is
bitwise unchanged. Original references are retained; finite-Gamma localization
comparisons apply the correction to both meshes. Public managed/detached/resume
verification and six maintained regression items passed. The high-frequency
managed check uses a captured pre-edit managed reference; its pre-existing
2.0655% difference from the detached mesh is reported separately in the diagnosis.
The engine has no scientific acceptance, retry or refinement loop.

Serial live log:

```sh
tail -n 80 -F .linecablemodels/fem/local-interface-mesh/live.log
```

`qualify_terminal_shift.jl` reproduces the finite-Gamma controls on frozen
meshes. `verify_execution.jl` checks the public managed/detached/resume routes.
`mesh_images.jl` writes static before/after figures. Results remain under
`.linecablemodels/fem/local-interface-mesh/`; these scripts are manual tools.

## Earlier conservative implementation (superseded)

The production feature uses the existing Gmsh size fields. A medium's full
interface refinement source is replaced with finite projected cable footprints
only when its remote element target already satisfies the existing wave-size
bound `mesh_size_factor/(8abs(q))`, with
`q = sqrt(im*omega*mu*(sigma+im*omega*epsilon)-Gamma^2)`.
Otherwise the full interface band is retained. Constant fields stay unchanged.
This prescribed rule runs once during construction; no solution acceptance,
retry, fallback, clipping or adaptive refinement is added to production.

For a cable at `(xc,yc)` with outer radius `r`, the finite interface segment has
halfwidth `interface_refinement_factor*(abs(yc)+r)`. Its distance field combines
with distance to the actual cable exterior, retaining the existing local sizes,
growth, medium restrictions and remote ceilings. The dimensionless factor
defaults to 1 and must be at least 1; increasing it retains more interface
refinement. Detached ONELAB exposes `InterfaceRefinementFactor` under Mesh.
Domain size, PML prescription, conductor targets and voltage paths are separate.
The same existing `q` includes the prescribed complex Gamma in both paths.

## Qualification performed before implementation

The frozen evidence is under `.linecablemodels/fem/local-interface-mesh/`.
`qualify.jl` prepares native bundles and changes only qualification meshes.
It never replaces production methods or performs voltage integration in Julia.
GetDP owns the complete field solve and voltage/matrix extraction. The existing
interactive `onelab-two-bare-wires` bundle was copied, never overwritten.

The final selection has 15 two-wire cases: air, buried, mixed; 0.1 Hz, 1 kHz,
1 MHz for Gamma=0, plus 0.1 Hz and 1 MHz for Gamma=0.99*gamma_earth. All use
0.1 ohm m earth, 4.25 cm wires and the existing domain/PML/conductor preset.
All raw R/X/G/B entries and analytical values are in each `components-*.csv`.
Where localization changed the fields, maximum component changes relative to
the preceding prescription were below 0.55%, with no new sign reversals.
This establishes preservation on these fixtures, not absolute physical accuracy.

Representative observations (native process time, one serial comparison each):

| Mixed case | Nodes before → after | DOFs before → after | Native seconds before → after | Largest component change |
| --- | ---: | ---: | ---: | ---: |
| 0.1 Hz, Gamma=0 | 101290 → 96410 | 400546 → 381026 | 30.21 → 29.05 | 0.0916% |
| 1 kHz, Gamma=0 | 76402 → 71553 | see saved solve records | 23.43 → 21.77 | 0.0359% |
| 0.1 Hz, Gamma=0.99 gamma_earth | 99556 → 90500 | 393673 → 357449 | 40.09 → 34.95 | 0.3565% |

Gains are modest because the PML remains unchanged. In these fixtures, 1 MHz
offers essentially no coarsening. Single-run timing differences on unchanged
meshes are measurement variation, not optimization gains. `assessment.csv`
records all 15 comparisons, including the cases with no change.

Screen, tubular and sector fixtures were checked at 1 MHz in 100 ohm m earth.
Their PML, contour and path coordinates matched exactly; all prescribed
conductor sizing/grading values matched. Total nodes changed 181632→176753,
67746→62887 and 147518→142655 respectively. Unstructured metal interiors were
not bitwise identical: triangle counts changed 163506→163502, 8756→8756 and
137188→137186. This is not a new numerical accuracy qualification of those
cables; it checks the retained geometry/resolution contracts.

## Important finite-Gamma limitation

In the mixed 1 MHz, Gamma=0.99 gamma_earth case, independent remeshing of the
**unchanged** size prescription changed a weak mutual X entry by 56.1% and G12
by 23.1%, with two extra nodes. PML/path/contour coordinates, medium triangle
counts and all size targets were unchanged. A pointless rewrite of a constant
field gave a 73.9% maximum component change. These failed comparisons remain
in `localized-wave-active` and `localized-wave-constant-rewrite`; they were not
counted as successful accuracy results or used to relax tolerances.

When the feature makes no field change, the qualification verifies matching
sources and uses the exact baseline mesh. All three finite-Gamma 1 MHz cases
then produced identical coefficients. These are **unchanged-mesh preservation
checks**, not independent remeshing/convergence checks. This feature does not
fix or certify the baseline finite-Gamma sensitivity. Tiny Gamma=0 conductances
also retain pre-existing relative disagreement with the analytical reference.

## Rejected constructions

Unconditional localization passed the initial mixed case but changed an aerial
mutual G by 1.62%, exceeding the declared 1% preservation target. Air-only
localization passed at 0.55%; soil-only localization changed G by 2.23%.
Widening both footprints by the existing decay radius still changed G by 1.45%.
An anisotropic BAMG attempt grew one physical face to 900000 vertices and was
stopped before GetDP. A fixed-mesh AMF ordering control changed components by
only 6.23e-7 relative, so it did not explain the low-frequency mesh sensitivity.
These experiments remain exclusively in the manual qualification harness.

## Run and inspect

One append-only live log contains Julia and line-buffered native GetDP output:

```sh
tail -n 80 -F .linecablemodels/fem/local-interface-mesh/live.log
```

`qualify.jl wave-spectrum` resumes the frozen qualification bundles and complete
meshes/solves. It uses the saved pre-feature native sources; do not delete them
and silently substitute current production exports as the old reference.
`assessment.jl` reads results without loading the FEM backend.
`preserve_shapes.jl` checks the saved cable bundles without solving.
`verify_execution.jl` exercises current public managed computation, detached
GetDP on the same mesh and explicit completed-run resume for zero and finite
Gamma. It records first-use compilation separately from elapsed time.

Maintained regressions are in `test/extensions/fem_interface_mesh.jl`, plus
the existing export, physical-PML, Gamma-transport and option-ownership tests.
No qualification gates or experimental meshing algorithms enter feature code.

## Implemented execution checks

The maintained selection passed 2408 assertions: options 160, existing editable
mesh/export controls 45, PML/path nodes 2135, Gamma transport 36, and the new
interface grading/native controls 32. The new test covers both zero and finite
Gamma and native edits from footprint factor 1 to 1000 while preserving PML
counts and conductor/path nodes.

Two public managed computations, their same-mesh detached solves and explicit
completed-run resume passed another 20 assertions. Production Γ=0 matched the
qualified coefficients within 1.71e-13 relative; production Γ=0.99γearth on a
fresh mesh differed by at most 0.1267%, retaining all component signs.

| Public mixed 0.1 Hz run | Elapsed | Julia compilation | Recompilation |
| --- | ---: | ---: | ---: |
| First call, Γ=0 | 61.6505 s | 31.7356 s | 0.0735 s |
| Warmed call, Γ=0.99γearth | 35.0255 s | 0.0204 s | 0 s |

These two calls solve different Γ cases, so their elapsed times are not a
before/after speedup comparison. `benchmark.jl` repeats the matched frozen
Γ=0 native pair, excluding Julia compilation. Results are saved under
`warmed-repeat/`. This repeat took 30.12 → 29.18 s (3.1% less native process
time), with peak RSS 2005408 → 1904052 KiB (5.1% lower). Assembly took
20.60 → 20.23 s and factorization/solve took 6.88 → 6.53 s. The original pair
took 30.21 → 29.05 s. These two observations support a small saving on this
case; they do not establish a general speedup across geometries/frequencies.
Public-run records and inspectable native `study.pro` bundles
are under `feature-execution/gamma-0.0/` and `feature-execution/gamma-0.99/`.

The usual `run_two_bare_wires_fem.jl` now lists `interface_refinement_factor=1.0`
and forwards it to its detached export. Restart Julia after this struct change
and re-export to update an older detached bundle; existing files are not edited
by the new package code until explicitly exported again.

Native field reference: [Gmsh mesh size fields](https://gmsh.info/doc/texinfo/#Gmsh-mesh-size-fields).
