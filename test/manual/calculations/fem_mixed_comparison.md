> Historical investigation record. Equations, observations and limitations are
> retained at their recorded source versions. Referenced prototypes and campaign
> launchers may have been retired; their commands and pending-work statements
> are not current execution instructions. See the [cleanup record](fem_development_cleanup.md)
> and [current PML controls](fem_fixed_pml_controls.md).

# Mixed-case analytical/FEM comparison diagnosis

## Confirmed cause of the large disagreement

The manual runner compares different matrix treatments. Its FEM formulation
sets `reduce_bundle=false`, `kron_reduction=false` and
`ideal_transposition=false`, but its analytical formulation is `Formulation()`.
The analytical defaults enable all three options. For this two-wire fixture,
the consequential difference is ideal transposition.

Source locations:

- `run_two_bare_wires_fem.jl:103`: explicit untransposed FEM settings.
- `run_two_bare_wires_fem.jl:247`: default analytical formulation.
- `src/engine/options.jl:15`: defaults, including ideal transposition.
- `src/engine/matrixops.jl:1`: cyclic averaging of matrix entries.
- `src/engine/lineparameters.jl:95`: transposition of Z and, after inversion, Y.

For a two-conductor matrix A, this operation returns

```text
[(A11+A22)/2   (A12+A21)/2
 (A12+A21)/2   (A11+A22)/2].
```

It therefore assigns half the buried self admittance to the aerial diagonal
and mixes the two different mutual responses. It also averages the self
resistances. Equal wires in the same half-space largely conceal this mismatch;
the mixed placement exposes it immediately. No reciprocity requirement was
used in this diagnosis.

Fresh analytical calculations from the exact saved FEM problem snapshots
confirm that the default result equals this explicit averaging of the
untransposed result **exactly**, for both Z and Y at every checked frequency
and resistivity. This isolates the option mismatch without a new FEM solve.

The narrow runner correction is:

```julia
analytical_formulation = Formulation(; options = (
    reduce_bundle = false,
    kron_reduction = false,
    ideal_transposition = false
))
```

Use the same physical mixed inputs when recomputing the analytical results.
At investigation time, the saved runner's two placements both used `vert`,
while the retained mixed runs below have heights (-1,+1) m. The diagnostics
decode each retained snapshot rather than rebuilding from the current runner.
The runner and production equations were not edited in this investigation.

## Exact retained cases and numerical comparison

All four runs use quasi-fw, copper radius 0.0425 m, positions (0,-1) and
(1,+1) m, relative earth permittivity/permeability 1, and ten logarithmic
frequencies from 0.1 Hz to 1 MHz. Their settings include 24 soil skin depths,
192 PML layers, conductor mesh factor 3 and exterior mesh factor 8.

| Soil resistivity, ohm m | Retained native run | Maximum complex-entry Y error |
|---:|---|---:|
| 0.1 | run-a3RDNf | 4.861% |
| 1 | run-RV2bFu | 2.723% |
| 100 | run-8a7AVQ | 1.540% |
| 1000 | run-OPUcin | 1.370% |

The error is `abs(Y_FEM[i,j]-Y_ref[i,j])/abs(Y_ref[i,j])`, maximized over
the four entries and ten frequencies in each run. It is not a componentwise
conductance tolerance or a convergence bound. For example, at rho=0.1 and
0.1 Hz, G22 is -6.33942e-25 S/m in FEM versus -9.69607e-25 S/m analytically:
the relative G error is 34.62% despite a small complex-Y error. At rho=1,
the sampled G21 near its zero has a 19.05% relative error. At rho=100 the
maximum componentwise G error is 2.049%, and B error is 1.523%.

An illustrative diagonal comparison at rho=100 ohm m and 1 MHz:

| Quantity, S/m | FEM | Analytical, matched settings | Analytical, original default |
|---|---:|---:|---:|
| G11, buried | +0.0111242252 | +0.0110465567 | +0.00552313463 |
| G22, aerial | -2.91234031e-7 | -2.87450224e-7 | +0.00552313463 |
| B22, aerial | +9.15072764e-5 | +9.06174591e-5 | +0.00160243109 |

## Returning sign changes

With matched settings, all **160/160 G signs and 160/160 B signs** agree
at the saved FEM samples. This covers four resistivities, ten frequencies
and four entries. No values are clipped, symmetrized or sign-corrected.

A fresh 359-frequency analytical sweep per resistivity resolves smooth
mutual zero crossings. Examples for G21 are:

| Soil resistivity, ohm m | Analytical zero bracket, Hz |
|---:|---:|
| 0.1 | 19054.6 to 19952.6 |
| 0.1 | 398107 to 416869 |
| 1 | 158489 to 165959 |

There are also G12 and mutual B crossings in the low-resistivity cases.
There are no mutual G/B crossings in this dense reference sweep for rho=100
or 1000. Ten samples over seven decades miss entire lobes of some curves;
connecting those samples makes smooth behavior look like abrupt switching.
The default transposition further mixes the two directions and moves their
apparent zeros. The supplied screenshot therefore does not establish a
return of the previous FEM-only sign defect or inconsistent root selection.

This is evidence of agreement with the analytical model, not independent
physical validation of every zero or convergence of the FEM zero locations.
The newly dense curves are analytical; the FEM points are the original ten.
These results concern quasi-fw, not the separate quasi-TEM model discrepancy.

## Remaining series self-resistance error

The transposition mismatch also affects the displayed analytical R diagonals.
After correcting it, a real low-frequency error of about 4.56% remains.
Reading the actual rho=100, 0.1 Hz mesh gives, for each conductor:

```text
meshed area       = 0.00541875000000000 m²
intended pi*r²    = 0.00567450173054657 m²
area deficit     = 4.507034%
copper rho       = 1.7241e-8 ohm m
rho/A correction = 1.43401662e-7 ohm/m
```

The observed aerial self-R difference is
`3.28028501184e-6 - 3.13712223823e-6 = 1.43162774e-7 ohm/m`.
The simple polygon-area estimate agrees with this difference to 0.17%.
This identifies conductor geometry resolution as the next
target for low-frequency R accuracy. It does not justify changing the PML,
voltage reference or GetDP global quantity extraction. The remaining general
Y errors still require a controlled local-mesh convergence study if tighter
percentage accuracy is required.

## Reproduction and artifacts

Evidence is saved at `.linecablemodels/fem/mixed-reference-evidence/`:

- `check.jl`: decode the four native snapshots and evaluate default and
  untransposed analytical Z/Y; no FEM solves.
- `compare.py`: read native P/Z TSVs, invert each full complex P, compare
  against the fresh analytical matrices, and write `comparison.csv`.
- `dense.jl`: evaluate the dense untransposed analytical spectrum.
- `plot.py`: render saved signed values with FEM markers and dense reference
  curves. The signed logarithmic axes retain nonzero values.
- `analytical.csv`, `dense.csv`, `comparison.csv`, and the corresponding logs.
- `mixed-G-matched.png`, `mixed-B-matched.png`, `mixed-crossings.png`, and
  `mixed-self-R.png`.

From the repository root on this host:

```sh
JULIA_DEPOT_PATH=/tmp/lcm-fem-depot:/home/amartins/.julia \
JULIA_LOAD_PATH=@:@stdlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
julia --startup-file=no --compiled-modules=existing --project=test \
  .linecablemodels/fem/mixed-reference-evidence/check.jl

python3 .linecablemodels/fem/mixed-reference-evidence/compare.py

JULIA_DEPOT_PATH=/tmp/lcm-fem-depot:/home/amartins/.julia \
JULIA_LOAD_PATH=@:@stdlib OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
julia --startup-file=no --compiled-modules=existing --project=test \
  .linecablemodels/fem/mixed-reference-evidence/dense.jl

MPLCONFIGDIR=/tmp/fem-pml-matplotlib /usr/bin/python3 \
  .linecablemodels/fem/mixed-reference-evidence/plot.py
```

The four retained native run directories are required for recollection.
The saved CSVs suffice for plotting. The Julia diagnostics intentionally
use the native snapshot decoder: `input/problem.json` is a native run
payload rather than the wrapped document produced by public JSON export.
