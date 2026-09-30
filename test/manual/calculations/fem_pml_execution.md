> Historical investigation record. Equations, observations and limitations are
> retained at their recorded source versions. Referenced prototypes and campaign
> launchers may have been retired; their commands and pending-work statements
> are not current execution instructions. See the [cleanup record](fem_development_cleanup.md)
> and [current PML controls](fem_fixed_pml_controls.md).

# FEM exterior correction: execution record

References below to prepared voltage paths and contour averaging describe the historical extraction implementation. Retained scientific results have not been rewritten. See the [native voltage contract](../../../docs/plans/fem-native-voltage-extraction-plan.md).

Historical status at 2026-09-26: implementation and acceptance testing in progress.
The full broadband correction is **not yet validated**. The dedicated checkout is
`/home/amartins/Documents/KUL/LinecableModels-fem-investigation`.
The acceptance specification is [fem_pml_design.md](fem_pml_design.md).
Existing unrelated plotting, documentation and test changes are preserved.

## Established issue and reproduction

The retained old quasi-fw run is `.linecablemodels/fem/runs/run-scf77b`.
At 100 ohm m and 1 MHz, its aerial G11 is −0.129370 microS/m versus
−0.288927 microS/m from the matching analytical formulation. Changing only
integration in the old air infinite shell moves G11 to −0.163416 microS/m;
changing only the earth shell gives −0.128895 microS/m. The real spherical-shell
Jacobian compresses an unbounded oscillatory air solution into a finite annulus.
These tests localize the sensitivity to the exterior discretization. They do
not establish a material-law defect or a square-root branch switch.

Both outer scalar boundaries currently impose v=0. Air already participates in
the electromagnetic domain through displacement current. Adding it to DomainC
with zero conductivity did not alter the diagnostic solution. Neither a ground
constraint on y=0 nor a physical air conductivity is part of this correction.

## Implemented replacement

The physical rectangle extends to ±L vertically and L on either side of the
layout centre, with L=max(layout clearance, domain_skin_depths*soil skin depth).
A finite Cartesian PML uses s=1+(1−j)b(d/t)^3 for the positive-time phasor.
The x stretch is shared across the air/earth interface. Conductivity,
permittivity, permeability and every affected weak-form term are transformed;
ordinary volume Jacobians replace the infinite-shell Jacobian.

Buried voltage paths remain vertical through the bottom PML. They integrate
the pulled-back edge field using mesh-coordinate line weights. Surface
references and circumferential quadrature for aerial receiving rows are
preserved. The PML endpoint is an approximation to the deep-earth reference,
whose convergence still has to be measured.

Execution controls are owned by the existing compute options: PML thickness,
thickness multiplier, layer count, reflection target, physical mesh-size factor
and volume quadrature. Immutable GetDP sources, including pml.pro, remain captured
per run; all effective controls participate in resume checks. Mesh-affecting
controls participate in mesh-cache identity.

Native field output additionally records material tags, a physical/PML mask,
transverse current components and their transformed mesh flux. PML values are
labelled as analytic continuations, not physical loss densities.

## Failures found while testing the candidate

1. At 0.1 Hz the initial coupled matrix produced NaNs. Its 2,664,855 assembled
   coefficients were all finite, with magnitudes from 6.31e−30 to 3.28e11.
   Disabling the stretch made the solve finite; changing only MUMPS pivoting did
   not. PETSc diagonal equilibration solved the original equations without NaNs.
   `-ksp_diagonal_scale_fix` restores the original matrix and RHS before residual
   evaluation and the next source solve. Native PETSc traces confirmed one LU
   factorization and two solves. Unscaled residual norms were about 1e−5 for a
   RHS norm of sqrt(2); these are reported, not interpreted as error bounds.
   See [PETSc diagonal scaling](https://petsc.org/release/manualpages/KSP/KSPSetDiagonalScale/).

2. A purely imaginary stretch did not adequately resolve evanescent modes.
   In the layered lossless TE control at k_air*L=0.21 and 256 layers, the error
   against its exact finite-PML solution was about 14% of the infinite-reference
   peak. Adding the positive real stretch reduced that to about 0.003%.
   The native strip campaign covers 540 field comparisons per profile: three
   wave scales, normal/45/80/89-degree and evanescent modes, conducting/lossless
   lower media, three resolutions, two thicknesses, and three scalar reductions.
   Grazing-wave hard-wall error remains distinct from discretization error;
   the nominal normal-incidence reflection target is not an all-angle bound.

3. The initial nested rectangular-ring mesh failed the 2D corner control.
   A deliberately strong shared side stretch in a homogeneous conducting medium
   gave about 18% relative field error near a bottom corner at 64, 128 and 256
   layers. Adding normal layers did not resolve Cartesian stretch transitions
   crossing tangential elements. The replacement mesh has ten Cartesian patches,
   aligned with both stretch onsets and y=0, and tensor-product corner resolution.
   At 64 layers the corrected mesh reduces the same maximum error to 0.4585%.
   At 256 layers it is 0.4531%; halving physical mesh sizes reduces it to 0.2705%.
   The 17,618 physical triangles are byte-for-byte coordinate-equivalent at
   64/128/256 PML layers. The remaining floor is therefore not removed by adding
   normal PML layers alone. A separate unstructured-corner trial was rejected:
   its low-frequency air error remained about 10% at 256 layers.

The current native cylindrical controls contain 792 baseline and 528
physical-refinement field/probe comparisons. The layered strip contains 540
comparisons for each of the old and current profiles. At 256 layers the worst
finite-PML strip error is 0.7403% of the infinite-reference peak, in a low-frequency
89-degree lossless TE case. Finite-wall grazing error is larger and is reported
separately; the normal-incidence reflection target is not a uniform angular bound.

## Measured cable results so far

The final Cartesian-patch sentinel completed all **72** case/frequency points
(four frequencies, both models, three layouts, three resistivities). At 1 MHz,
quasi-fw aerial G11 errors against unified are +0.188%, +0.5693%, and +0.0704%
for 1, 100, and 1000 ohm m, respectively.

At 1 ohm m and 166810.05372000593 Hz, baseline quasi-fw aerial G11 is
−3.4133623482e−10 S/m versus −3.5021212921e−10 analytically (−2.5344%).
Increasing PML layers from 128 to 256 gives −3.4955637848e−10 (−0.1872%).
Halving physical mesh sizes at 128 layers gives −3.4029025160e−10 (−2.8331%).
The doubled-thickness/256-layer run gives −3.4844280389e−10 (−0.5052%).
Changing triangle quadrature from 12 to 13 or 16 leaves the baseline result
unchanged to the displayed precision. Halving/doubling total attenuation gives
−2.8075%/−2.2683%, respectively. These controls identify remaining PML discretization
sensitivity; selecting an attenuation value to fit this point is not justified.

The 99-frequency baseline began at 00:46 on 2026-09-26. At 07:33 it had completed
all twelve aerial/buried cases (1,188 case/frequency points); mixed cases started.
The complete PML-refined sweep started at 07:11, with physical-refined and
doubled-thickness sweeps queued after it. Each phase has its own `runs.csv`.
The full field-map campaign is also running at the four sentinel frequencies.

The dense baseline exposes a genuine mutual-G zero crossing in buried 1 ohm m
soil near 794 kHz. At 794328.2347242822 Hz, quasi-fw G12 is −6.165781e−4 S/m,
unified is −3.217967e−4, and the independent finite-electrode scalar reference is
−2.547128e−4. A roughly 2.95e−4 absolute difference therefore gives a roughly 92%
relative error against the small unified value. Physical refinement of that
point and both neighbors is running. The zero crossing is shared by the
references and FEM; the error near it is not evidence of a branch-sign jump.

The h/2 control subsequently moved quasi-fw G12 at the crossing to
−5.226619e−4 S/m and quasi-tem to −3.901134e−4. Further h/4 and combined h/2,N256
controls are running. The first N256 low-frequency aerial result also changes
rho1 G12 at 0.1 Hz from +3.392161e−24 to −3.836142e−24 S/m, against
−5.568393e−24 analytically. It recovers the reference sign, but about 31% relative
error remains; this is not demonstrated relative convergence.

At 1 MHz, the six completed quasi-fw aerial/buried self-Z comparisons have
resistance errors between −0.31% and −0.12%, and reactance errors between
−0.54% and −0.43%. These terminal comparisons do not establish pointwise
resolution of conductor skin currents: the existing conductor target remains
geometry-based (disk radius/5), independently of the new exterior correction.

An independent finite-circle single-layer reference for the scalar PDE has now
completed all 891 layout/material/frequency combinations. Four-frequency controls
independently refine 64 to 128 contour panels and 256 to 512 spectral points.
For resolved entries above 1e−10 S/m, maximum relative changes are 9.501e−7 and
2.138e−12, respectively. Its interface kernel is derived from continuity of v
and kappa times its normal derivative, not from the package analytical kernel.
It subtracts the interior continuation current before forming terminal Y and
retains surface averaging for aerial rows. Its scalar-model differences from
unified remain plotted separately from FEM/reference differences.

Current full-spectrum raw data and plots are under
`/tmp/fem-pml-execution/final-baseline/`; scalar reference data are under
`/tmp/fem-pml-execution/scalar-reference-full/`. The aerial 1 ohm m low-frequency
mutual G is still sign-unstable at magnitudes around 1e−19 S/m and smaller.
The absolute 1e−12 criterion does not demonstrate that these values are resolved.

### Superseded ring-mesh results

The earlier ring-mesh candidate with diagonal scaling completed all 24 sentinel
points (four frequencies, two physics choices, three layouts) at 100 ohm m.
At 1 MHz its quasi-fw aerial G11/G12 were −0.290649/−0.191718 microS/m,
versus −0.288927/−0.189898 microS/m analytically. The real-plus-imaginary profile
on that same mesh gave −0.290627/−0.191695 microS/m. These are encouraging
candidate results, not acceptance of the final Cartesian-patch implementation.

For the ring candidate, quasi-fw buried G11 differed from the analytical value
by about 0.54–0.66% at the four sentinel frequencies. Mixed entries also remained
finite; one high-frequency mutual G error was about 1.26%, so physical-mesh
refinement remains necessary. The quasi-tem aerial scalar model had a much
larger analytical difference at high frequency; PML convergence must be judged
against that scalar PDE as well as displaying its analytical difference.

Raw matrices, reference matrices and plots are under
`/tmp/fem-pml-execution/sentinel-scaled/`. Plots show four sentinel points explicitly;
they are not the 99-frequency sweep. CSV files contain signed non-thresholded data.
The plotting script displays the signed-log linear threshold and also writes
linear and signed-error plots. No custom JSON summary is needed for the plots.

## Reproduction tools

Run from this checkout with `julia --project=test`. On this managed host the
working depot overlay is `/tmp/lcm-fem-depot:/home/amartins/.julia`.

- `run_fem_pml_validation.jl`: 99 frequencies, both physics choices, all three
  layouts and 1/100/1000 ohm m. `PML_PHASE` selects baseline, pml_refined,
  core_refined or thicker. `PML_OUTPUT` chooses the evidence directory.
  `PML_FREQUENCIES`, `PML_RHOS`, `PML_LAYOUTS`, `PML_PHYSICS`, `PML_MAPS` and
  `PML_WORKERS` select bounded controls. Completed and interrupted run paths are
  retained; reruns use the existing validated resume contract.
- `run_fem_pml_layered_modes.py`: independent closed-form finite-PML and
  infinite-half-space references for the native Fourier-mode controls.
  Requires numpy and `LINECABLEMODELS_GETDP` pointing to GetDP.
- `mesh_fem_pml_waves.jl`, followed by `run_fem_pml_waves.py`: native cylindrical
  K0 controls on the production 2D meshes, including physical corner probes.
  Set the same `PML_WAVE_OUTPUT` for both; Python additionally needs scipy.
- `plot_fem_pml_validation.py OUTPUT [VARIANT_OUTPUT ...]`: CSV comparisons and
  PNG/SVG plots, using numpy and matplotlib. It does not alter plotted values.
- `run_fem_scalar_reference.py`: independent finite-electrode scalar reference;
  set `PML_FREQUENCY_FILE` to a campaign's `frequencies.csv` for the 99-point grid.
  `PML_BEM_PANELS` and `PML_BEM_QUADRATURE` accept comma lists for refinement.
- `plot_fem_pml_fields.py RUN FREQ_INDEX BASIS OUTPUT`: physical close-ups,
  exterior magnitude/phase plots, and a vertical cut read directly from native
  `.pos` files. The cut retains separate element values across interfaces.
- `fem_pml_voltage_tail.py RUN FREQ_INDEX LAYOUT OUTPUT.csv`: scalar/vector
  contributions and resulting P/Y at physical and bottom-PML reference depths
  for the standard buried/mixed fixture. An independent affine gauge-cancellation
  check passed at five depths; the native buried map audit remains queued.
- `run_fem_pml_preservation.jl`: original rho0.1/radius controls, asymmetric
  and permittivity variants, and reversed mixed conductor order at the ten
  historical frequencies, in both models.
- `run_fem_domain_size.jl`: physical-domain comparisons with PML thickness and
  prescribed mesh-size targets held fixed; both models, all layouts/resistivities
  at decade endpoints and historical high-frequency points by default.

## Regression and plotting checks

The final Cartesian geometry passed 287 structural/option/resume checks and
260 native formulation, path, map and multi-frequency checks. Domain enlargement
passed 144 assertions that prescribed mesh targets and stretch remain fixed.
The public plotting API passed 22 signed-data/rendered-transform checks with
`clip=false` and with zero G cutoff. Toggling signed log off restores identity;
small nonzero values then flatten visually under a large linear dynamic range.
An actual desktop GL check could not open X display `:0` in the sandbox. The
subsequent desktop-access request was not executed because automatic approval
review failed authentication. Desktop GL verification remains outstanding.

## Practical aerial preset, 2026-09-26

The two-wire manual runner now selects `domain_skin_depths=24.0`,
`pml_layers=192`, and `mesh_size_factor=3.0`. The larger physical domain moves
the air PML away from the conductors; coarser physical elements keep the cost
manageable. This multiplier also coarsens conductor contours, so this is a
practical preset for conductance signs, not a percentage-converged mesh.
The runner retains two frequency workers and four solver threads.

Eight targeted quasi-fw aerial case/frequency checks completed with all **32/32**
raw `real(Y)` entries matching the analytical reference signs. The reference is
the unified formulation with Gamma=0. Heights and separation are both 1 m.

| Soil resistivity [ohm m] | Radius [m] | Frequency [Hz] | FEM G11 [S/m] | FEM G12 [S/m] |
| --- | --- | --- | --- | --- |
| 0.1 | 0.0425 | 0.1 | -1.339862e-24 | -3.570835e-25 |
| 0.1 | 0.0425 | 21.544347 | -8.774694e-20 | -4.206735e-20 |
| 0.1 | 0.0425 | 774.263683 | -2.759767e-16 | -2.167401e-16 |
| 0.1 | 0.0425 | 4641.588834 | -1.947818e-14 | -1.725784e-14 |
| 0.1 | 0.0425 | 1000000 | -1.696341e-8 | -1.641989e-8 |
| 0.1 | 0.085 | 0.1 | -1.985833e-24 | -3.612615e-25 |
| 100 | 0.0425 | 1000000 | -2.918366e-7 | -1.931141e-7 |
| 1000 | 0.0425 | 1000000 | -1.804422e-6 | -8.327514e-7 |

The 8-skin-depth/256-layer/unit-mesh candidate still gave positive mutual G at
0.1 Hz for rho=0.1. The 16-skin-depth/192-layer/mesh-factor-2 candidate passed
the standard radius but failed the 0.085 m radius. These rejected candidates
are retained with the selected preset's evidence under
`/tmp/fem-pml-execution/manual-sign-diagnosis/`.

Native solves with one solver thread took 164-625 seconds and reported peaks
of 3006-10540 MB. The two intermediate frequencies ran concurrently in 652
seconds including preparation. These measurements are not a benchmark of the
manual runner's four-thread setting. Expect hours for its full 70-point sweep.
That full sweep, other placements, and quasi-tem have not been revalidated
with this preset; this result does not complete the broader acceptance plan.

`practical-settings.md` records the selection and reproduction command;
`practical-signs.csv` contains all four entries and reference errors;
`practical-costs.csv` contains native costs and retained run names.
`practical-conductance.png` compares the previous manual run, analytical values,
and the five tested standard-radius rho=0.1 frequencies. New results appear as
markers only. The manual runner passed syntax and whitespace checks.

## Checks still required for completion

Completion of all four full cable sweeps, broad attenuation/domain/quadrature
controls, physical-depth voltage-tail checks, selected field maps and remaining
preservation controls is still required. Desktop GL verification remains
unavailable as described above. No full-suite or broadband pass is claimed.
