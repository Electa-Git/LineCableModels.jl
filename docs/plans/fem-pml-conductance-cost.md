# Recover aerial conductance accuracy at lower exterior cost

Status: completed on 2026-09-29, including numerical qualification, sixteen
figures in SVG and PNG, and delivery checks. The user authorized end-to-end autonomous
execution with one live log and no per-trial approvals.
Evidence and the append-only log are under
`.linecablemodels/fem/pml-conductance-cost/`.

Results: `(144,144,96)` passed all 363 managed conductance sign comparisons
in 79 frequency solves / 167 source columns. One independently meshed detached
solve passed signs; same-mesh comparison against the saved managed-worker solve
passed unchanged parity bounds with identical Z/P. Independent-mesh G differed
by up to 7.08e-26 S/m and did not satisfy the same-mesh bound; both results are
retained separately. See the [delivery assessment](../../test/manual/calculations/fem_pml_conductance_cost.md)
for errors, timing scope and the prescribed controls. No production engine
code was changed by this execution.

## Objective and scope

Find a cheaper prescribed exterior discretization that preserves the observed
quasi-fw conductance signs across the user's ordinary studies. The immediate
accuracy objective is the user's sign requirement; conductance magnitude and
all separate R/X/G/B errors must still be reported. The 192-layer result is a
comparison, not an exact solution. Aim for at least 30% lower native solve time
than 192/192/192 on matching cases; report the measured saving without promising
that target in advance.

Keep conductor geometry/grading, physical domain, materials, excitation,
voltage references, native voltage integration and equations fixed initially.
Use Julia and native Gmsh/GetDP. Preserve managed and detached execution.
Qualification and any selection of settings belong entirely in the manual
harness. Production will have no sign correction, clipping, scientific
acceptance, retry, fallback or automatic refinement.

## New evidence obtained while planning

At radius 0.085 m, soil 0.1 ohm m and 21.544346900318832 Hz, the measured change
`real(P192-P96)` in ohm m is

```text
[-0.015708534  -0.015706576
 -0.015706576  -0.015704670]
```

Removing the best equal-entry offset leaves only 8.70e-5 of this matrix's
Frobenius norm. In the orthonormal equal-current/opposite-current basis, its
diagonal changes are -0.0314132 and -2.58e-8 ohm m: approximately a factor of
1.22 million apart. This is a strong lead toward shared exterior/common-mode
error, not proof of a particular boundary mechanism. It does not authorize
subtracting the offset from results. Individual scalar-potential outputs are
gauge-dependent; assess the complete scalar-plus-reference-plus-vector voltage.

The earlier source-expression controls exclude the PML expression optimization
as the cause of the sign reversal. Increasing bottom layers alone from 48 to 96
does not resolve the aerial case. Existing data and native meshes are reusable.
See [the reproduction and controls](../../test/manual/calculations/fem_pml_low_frequency.md).

## Phase 1: three discriminating solves

Use the failing 0.085 m aerial fixture at 21.544346900318832 Hz, with both source
columns. Retain the complete original problem/frequency definitions when
exporting geometry; solve only the selected frequency. Compare to saved 96 and
192 baselines rather than repeating them.

| New control | What changes | Question answered |
|---|---|---|
| Quadrature | On the saved 96 mesh, change triangle integration from 12 to the supported 13 points | Does integrating the rapidly varying stretched tensors more accurately remove a substantial part of the common-mode error without more unknowns? |
| Side resolution | Set counts to (192,96,96); retain grading and all other targets | Is side resolution responsible? |
| Top resolution | Set counts to (96,192,96); retain grading and all other targets | Is top resolution responsible? |

Report all P and Y entries, their separate real/imaginary parts, common and
differential changes, DOFs, assembly/factorization time and peak memory. Native
scalar/reference/path contributions may be printed separately in copied test
inputs to localize changes; the complete voltage remains authoritative.
Request field maps only for the most informative control if coefficient data
leave an ambiguity. Physical/PML interfaces and conductor nodes are checked
against the baseline so an unintended local remesh is not mistaken for a PML
effect. Shared interface constraints and side stretching remain conforming.

This is three new factorizations/six source columns. A change of sign alone is
insufficient: compare the distance to the analytical component and the 192
result. One 13-point check can identify sensitivity, but cannot prove quadrature
convergence. If quadrature dominates, confirm its benefit at the second probe
before treating it as the cheaper route.

## Phase 2: at most four additional solves

Use the Phase 1 result to pursue one explanation. If an asymmetric count choice
already meets the cost/accuracy objective, spend this budget confirming it,
not searching gratuitously. Otherwise retain 96 intervals and test two node
distributions, initially bracketing the current exponent 7.37535 with 6 and 9
in the implicated side/top direction(s). Keep the bottom distribution fixed.
Increasing the exponent concentrates nodes near the physical/PML interface;
decreasing it reallocates nodes toward the exterior. The PML stretching law,
thickness and absorption strength stay fixed, isolating discretization.

Before solving, inspect native curve coordinates, adjacent step ratios and
changes of stretched coefficients across elements. Use the same node count;
do not invent a new meshing framework or optimize G against its analytical value.
Use the existing `pml_grading` control and native Gmsh transfinite curves.

Each of the two distributions is checked at 0.1 Hz and 21.544346900318832 Hz:
four solves/eight columns. Keep both self and mutual values. No case-specific
offset, sign clamp, analytical replacement or hidden frequency switch is allowed.

**Decision after no more than seven new diagnostic solves:** either nominate
one cheaper setting with supporting evidence or reject this approach. If the
error does not respond to quadrature or node redistribution, do not launch a
larger sweep. Use the common-mode evidence to propose a separate exterior
operator/stretch-profile design with its mathematical and native-export
requirements made explicit. The seven-solve budget does not promise a fix.

Execution decision (2026-09-29): quadrature and both grading controls failed;
side-only and top-only refinement helped at 21.54 Hz, but top-only failed at
0.1 Hz. The sixth solve completed that rejection. The measured changes were
consistent with approximately inverse-square convergence in side/top counts,
so the seventh solve tested `(144,144,96)` at 0.1 Hz. It passed. Under the user's
instruction to continue autonomously, one additional targeted solve confirmed
21.54 Hz before the final sweep: eight diagnostic solves in total. This is a
recorded extension of the initial bound, not an open-ended mesh search.

Both 144 probes have the required signs, but their G12 magnitude errors remain
about 84% and 46%, respectively. At 21.54 Hz the native solve took 36.01 s,
versus 67.01 s for the matching current-source 192 control. The mesh audit found
identical conductor coordinates at 0.1 Hz and identical physical/PML interfaces
at both probes. At 21.54 Hz, independent export meshing moved 17 interior nodes
in conductor 1 (maximum nearest-node difference 1.35 mm); conductor 2 matched.
Therefore that probe is not claimed to isolate PML changes perfectly. Final
qualification uses ordinary managed computation and checks detached execution.

## Phase 3: one candidate, actual user cases

Only a candidate that survives both low-frequency probes proceeds. Qualification
uses the same fixed setting across fixtures; no per-case tuning or calibration.

- Run the ordinary two-wire study unchanged: four resistivities at radius
  0.0425 m and three radii at 0.1 ohm m, with the existing ten frequencies:
  70 frequency solves. Reuse compatible saved results where available.
- Check the existing three-conductor fixture in buried and mixed placements
  at 0.1 Hz, 21.544346900318832 Hz and 1 MHz: six frequency solves. Keep every
  energization/observation pair and explicit terminal order.
- Retain one 1 MHz preservation point for each screen, tubular-shell and sector
  fixture: three frequency solves. Compare with saved conductor qualification
  results and retain their existing scientific limitations.

This final selection is at most 79 managed frequency solves / 167 source columns,
plus one independently meshed detached two-wire solve / two columns for export
preservation: 80 solves / 169 columns total, for one candidate only.
Inspect each layout as it completes; stop this qualification
batch at a demonstrated sign regression. No new full 192-layer sweep is planned.
Use existing 192 results only when physical inputs and source differences are
accounted for. Any missing matched reference is identified before launch.

Judge signed G directly against matching analytical references where applicable,
and against independent saved/refinement evidence for explicit cable shapes.
Show errors of R, X, G and B separately; complex norms cannot establish a tiny
conductance's accuracy. Preserve the user's existing scientific tolerances;
introduce no new percentage criterion or near-zero clipping. Wrong signs on
the qualification fixtures fail the candidate. Failure means the cost/accuracy
objective remains unmet, even if a smaller mesh runs faster.

Measure native stage times at matching thread counts. Separate Julia first-use
compilation from warmed execution. Repeat a representative matched timing only
if a claimed saving is unclear. Estimate the final batch duration from the
diagnostic runs before launching it.

## Delivery and execution protocol

If existing controls suffice, deliver a documented prescribed configuration and
focused regression fixtures. Add no production abstraction. If qualification
demonstrates a necessary code change, implement only that change afterward,
through the existing mesh/source owners, and check managed/detached parity on
the same problem. The one detached bundle is independently meshed and solved
using the exported geometry and shared native formulation, at the failing
0.085 m aerial fixture's 21.544346900318832 Hz sample.

One implementation agent, serial diagnostic solves, ordinary native verbosity,
a single live log, per-case completion records and a saved input manifest.
After launch, verify startup and return the log path/command. Resume completed
cases after interruption; do not spend turns polling or rerun completed grids.
No Python, background agents, new monitoring framework or production scientific
validation machinery is needed. Existing solver frequency workers may be used
only for the final approved ordinary study at its established resource settings.

## Technical basis

- Steven G. Johnson, [Notes on Perfectly Matched Layers, sections 6 and 7.1](https://ocw.mit.edu/courses/18-303-linear-partial-differential-equations-analysis-and-numerics-fall-2014/0ad128a4b3d9dbb860e83a59a47b1b01_MIT18_303F14_pml.pdf): distinguishes evanescent behavior and discretization reflections from continuous matching.
- [COMSOL PML implementation](https://doc.comsol.com/6.4/doc/com.comsol.help.comsol/comsol_ref_definitions.21.137.html): discusses allocating resolution for mixtures of propagating and evanescent components. Its stretching-curvature control is not numerically identical to our mesh-grading exponent; no parameter value is copied from it.
- Druskin, Guddati and Hagstrom, [On generalized discrete PML optimized for propagative and evanescent waves](https://arxiv.org/abs/1210.7862): supports considering the discrete exterior response if ordinary grading fails. Its optimal-grid construction is a separate design, not an unqualified drop-in implementation for this coupled layered problem.
