# Provenance of the historical thin-wire FEM comparison

The supplied `thin_wire.pro` leads to the user's historical results. A fresh
GetDP run on 2026-09-09 reproduces the archived 1 MHz electric and magnetic
outputs to numerical precision. The earlier custom PEC quasi-TEM audit solves
a different electric equation and is withdrawn as a reference for this test.

The source directory is
`/home/amartins/Documents/KUL/LineCableModels-temp/@INPROGRESS/thin-wire-mutual`.
No source files in that directory were modified.

## Source and output trail

1. `thin_wire.pro:1` includes `thin_wire_data.geo`, which defines geometry,
   materials, frequency, source selection and output directories.
2. `thin_wire.geo` includes the same data and creates the air/earth geometry.
   `jacobian_integration.pro:20` applies `VolSphShell` in the outer shell.
3. `thin_wire.pro:163` routes Electric/Y to `electrodynamic_extended_Y.pro`,
   Electric/P to `electrodynamic_extended_P.pro`, and Magnetic to
   `darwin_formulation.pro`. Hybrid is a separate branch.
4. `ZYcomparisons.xlsx`, sheet `ZY-equal_radii-0.1`, contains all eight historical
   rows: GetDP is column F and COMSOL is column G. Self/mutual Z occupy rows
   6–13/17–24; self/mutual Y occupy rows 28–35/39–46. These exactly contain the
   values supplied in the conversation, including the analytical comparisons.
5. `results/electrodynamics/Y.dat` preserves the unrounded historical 1 MHz
   source-1 self and mutual Y. `results/darwin/Z.dat` preserves two magnetic
   source experiments. Those are archived outputs, separate from the new runs.

The currently saved data selects **Hybrid**, not Electric. The saved
`thin_wire.msh` also contains 10 mm insulation annuli, inconsistent with the
current zero-insulation data. Therefore the replay explicitly selects the
branch and regenerates the bare geometry at each frequency. Reusing that stale
mesh with bare-domain definitions produced zero electric output because the
excluded annuli disconnected the electrodes; that trial is retained in
`saved-mesh-1MHz` and is not a physical result.

## Electric equation and extraction

Write `s = jω` and `σ̂ = σ + sε`. In
`electrodynamic_extended_Y.pro:33`, `k² = -sμσ̂`. Its weak form at lines 87–105 is

```math
\int_\Omega \hat\sigma\nabla v\cdot\nabla w
       +s\mu\hat\sigma^2 v w\,d\Omega
+\sum_p Q_p w_p=0.
```

The corresponding field equation is
`−div(σ̂ grad v) + sμσ̂²v = 0`. This is the scalar Helmholtz/diffusion model
implemented by the Electric branch; it is not a derivation of equivalence to
the manuscript's full Maxwell fields and vertical-path voltage.

The integration domain contains air, earth and their infinite shells. The
wire interiors are excluded. Grouped nodal degrees of freedom make each wire
an equipotential electrode. The source is 1 V, the other electrode is 0 V,
and both air-side and earth-side outer boundaries are at 0 V. Insulation is
absent in the replay. The conductor conductivity placeholders do not enter
this electric operator.

`electrodynamic_extended_Y.pro:173` returns `Y = -Q`, in S/m for the unit
voltage excitation. Each output row lists the electrode responses for one
source; the replay stores response rows and source columns in its matrices.

The P branch drives `Q = -1` on the source and zero on the other electrode,
then returns the electrode voltages. Its quantity named P is therefore
`P_current = Y⁻¹`, in Ω·m, rather than the manuscript's `Pe`, in m/F. This is
also what `p_to_y.py` implements. For the manuscript convention one would set
`Pe = s P_current`; no extra `s` multiplies `inv(P_current)` when computing Y.

## Difference from the earlier audit

`test/gauntlet/getdp/pec_boundary.pro:67` used

```math
\int_\Omega \hat\sigma\nabla\psi\cdot\nabla w
       +s\hat\sigma a w\,d\Omega=\sum_p I_p w_p.
```

Here the electric equation is driven by the magnetic solution `a`; it does
not contain the Electric branch's `sμσ̂²ψ` term. It also imposes the electric
reference only on the earth-side far boundary. Its output `inv(ψ/I)` is not
the same boundary map as the direct Electric electrode-reaction experiment.
The production `quasi-tem.pro` likewise uses coupled `a/ur/phi` equations.
Consequently neither earlier curve reproduces this historical Electric test.
Mesh refinement or multiplication by `jω` cannot make these operators identical.

The historical Magnetic branch is also worth distinguishing: it assigns
`σ = 10¹² S/m` to the active wire and zero to the inactive wire, integrates
the wire interiors, prescribes active current 1 and inactive current 0, and
returns `Z = -U`. Replaying it establishes provenance; it does not establish
an exact exterior solution with both wires represented as PEC boundaries.

## Fresh execution and numerical checks

Eight meshes and both source columns of Y and Z were run at
0.1, 1, 10, 100, 1,000, 10,000, 100,000 and 1,000,000 Hz. At 1 MHz both
current-driven electric columns were also solved: 34 GetDP solves in total.
Source `.pro` and `.geo` files were copied without changing equations.
The programs were ONELAB Gmsh `4.14.0-git-67db5bd93` and GetDP
`3.6.0-git-1cf7fa06`, with one solver thread per process.

Geometry: radii 0.0425 m; centres `(−0.5,−1)` and `(0.5,−1)` m; no insulation;
earth resistivity 0.1 Ω·m; relative permittivity and permeability 1.
The original mesh controls were retained: conductor size `r/10`, interface
size 0.0255 m, physical radius `max(5 m, skin depth)`, shell radius `1.25R`.

At 1 MHz the new source-1 results are:

| Quantity | Fresh GetDP result |
|---|---:|
| Z self (Ω/m) | 0.752382315537412 + j1.569339425587136 |
| Z mutual (Ω/m) | 0.000160454659251 + j0.001069205027862 |
| Y self (S/m) | 40.90908104799158 + j19.61318039135062 |
| Y mutual (S/m) | −0.020709917700088 − j0.019107515924300 |

These reproduce the unrounded archived source-1 results. The spreadsheet's Z
entries have slightly different rounding/values from its raw archive; the
replay matches the raw archive. At 100 kHz the new Y also agrees with the
spreadsheet to its displayed precision. Lower-frequency Y entries are close
but do not reproduce every historical digit; the largest complex relative
difference in mutual Y is approximately 0.59%, at 1 kHz. The exact historical
meshes/settings for those rows have not been recovered. This is a reproduction
check with current saved source, not a mesh-convergence study or a COMSOL run.

Electric reciprocity holds within `3.3e-16` relative matrix norm across the
eight frequencies. At 1 MHz `‖P_current Y − I‖ = 7.8e-15` and the relative
difference between direct Y and `inv(P_current)` is `5.8e-15` or less.

The replay, source hashes, commands, solver logs, matrices and spreadsheet
extraction are preserved under `.linecablemodels/qa/thin-wire-trail/`:

- `replay.py`: driver; completed frequency results are reused on repeat runs.
- `source-manifest.json`: source location and SHA-256 hashes.
- `regenerated/`: copied original equations, new meshes and per-source outputs.
- `replay.json`: full 2×2 Y and Z matrices, plus 1 MHz `P_current`.
- `comparison.tsv`: new self/mutual values alongside archived GetDP and COMSOL.
- `historical-sheet4.json`: the source spreadsheet cells, without recomputation.
- `admittance-provenance.png`/`.pdf`/`.svg`: self and mutual Y comparison against
  the archived GetDP/COMSOL values and the earlier quasi-TEM audit.

Use this Electric branch as the historical comparison baseline. A separate
mathematical equivalence check is still needed before calling it a validation
of the manuscript's full field/path formulation.

## Proposal and Xue compared with the replay

The initial provenance plot omitted the analytical results. The corrected
comparison is produced by `test/gauntlet/plot_thin_wire_fem_comparison.jl`:

```sh
julia --startup-file=no --compiled-modules=existing --project=. test/gauntlet/plot_thin_wire_fem_comparison.jl
```

It reads the independently evaluated analytical earth matrices from
`.linecablemodels/qa/earth-matrices-manual/matrices.json` and the fresh FEM
`replay.json`. It does not feed analytical results into the FEM. The analytical
curves contain 281 quadrature samples; the FEM and archived COMSOL markers
contain eight frequencies. These curves include only earth Z and Y.

Blue is the proposed field average with the isolated primary-current factor Cf
used in the user's table. Magenta is the proposal with the full manuscript
current map L. The dotted dark curve is Xue. Both proposal normalizations are
shown explicitly, because their mutual entries differ appreciably at 1 MHz.

Relative complex differences from the fresh original GetDP at 1 MHz are:

| Entry | Proposed field average, Cf | Proposed full manuscript, L | Xue |
|---|---:|---:|---:|
| Z self | 0.2183% | 0.2183% | 12.4920% |
| Z mutual | 9.2942% | 0.8584% | 24.2629% |
| Y self | 0.2178% | 0.2178% | 13.0688% |
| Y mutual | 9.5976% | 1.1622% | 1.1627% |

These are discrepancies against the stated FEM operators and meshes, not
converged estimates of analytical error. The magnetic branch's active/inactive
conductor treatment described above still applies.

Outputs under `.linecablemodels/qa/thin-wire-trail/` are:

- `proposal-xue-fem-ZY.{png,pdf,svg}`: real and imaginary self/mutual Z and Y.
- `proposal-xue-fem-Z.{png,pdf,svg}` and `proposal-xue-fem-Y.{png,pdf,svg}`:
  separate figures for easier inspection.
- `proposal-xue-fem-relative.{png,pdf,svg}`: relative complex differences at the
  eight FEM frequencies, exposing differences hidden by overlapping curves.
- `proposal-xue-fem.csv`: every analytical and FEM matrix entry, with relative
  complex differences, for both Z and Y and all three analytical curves.
