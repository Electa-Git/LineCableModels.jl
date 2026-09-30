# Fixed-system corner diagnostic — 2026-09-30

The four prescribed executions completed. Three native MUMPS corrections barely
changed the structured result, but substantially changed the corner candidate.
For the candidate, the first-source residual worsened and the second-source
residual improved. The corrected candidate acquired positive G11 and G21.
This demonstrates sensitivity to the algebraic solution procedure on an unchanged
assembled problem. It does not identify the exactly solved discrete output,
establish convergence, or separate all remaining discretization and algebraic error.
The tested refinement setting is not a supported production correction.

## Fixed scope and retained evidence

- Two existing meshes: `baseline` and `metric-1.35`, aerial, 0.1 Hz, Gamma=0.
- Two configurations per mesh; both sequential source columns. No remeshing.
- Original GetDP 3.5.0 executable, PETSc 3.14.4 complex, MUMPS 5.3.3.
- Saved `.pre` files loaded using `-cal`: no regeneration of the gauge tree or DOF ordering.
- Original solver controls retained, including diagonal scaling and restoration,
  and `GenerateRHSGroup` / `SolveAgain` factorization reuse for the second source.
- The sole correction setting was `ICNTL(10)=-3`: exactly three native refinement
  steps. Original effective `ICNTL(10)=0`. Both used identical diagnostic controls
  `ICNTL(4)=3`, `ICNTL(11)=2` and `-ksp_view`.
- Effective `CNTL(2)=1.49012e-8` was unchanged. It is inactive for the prescribed
  negative ICNTL(10); termination was the fixed three-step count. Native
  `INFOG(15)` confirms 0 or 3 corrections separately for both field source solves.
  No KSP tolerance was substituted. The global flag also reached the existing
  small native P-inversion solves; independent high-precision inversion checks
  distinguish those rounding effects from the observed field-output changes.

The [live log](../../../../.linecablemodels/fem/pml-corner-algebraic/live.log)
contains all four runs and the sparse arithmetic assessment. Results and commands
are under `.linecablemodels/fem/pml-corner-algebraic/{baseline,metric-1.35}/{original,refine3}`.
The [complete numeric assessment](../../../../.linecablemodels/fem/pml-corner-algebraic/assessment.toml)
includes every complex P and Y entry, per-column diagnostics, scaling drift and
inverse-difference contributions. Original retained experiment files and feature
sources were preserved: 105 recorded files have unchanged hashes.

The diagnostic scripts are [check_algebra.jl](check_algebra.jl) and
[assess_algebra.jl](assess_algebra.jl). The latter performs arithmetic only.
No further solver setting or fixture is scheduled; an inconclusive accuracy
assessment is a completed outcome of this experiment.

## Original-coordinate system checks

Native system dumps saved the sparse matrix before the first scaling operation.
Binary `.res` records preserve b1, x1, b2, x2 without decimal serialization loss;
the temporary RHS-to-solution copy was restored before each solve. All comparisons
use the untouched original matrix A0 and each saved pre-scaling RHS.

For each mesh, the following match exactly between original and refine3:

- original binary matrix coefficients and sparse ordering;
- saved preprocessing, constraints and DOF ordering;
- both binary RHS vectors;
- the working matrix immediately before the second source solve.

The restored working matrix is not bitwise identical to A0. The largest coefficient
change divided by the largest original coefficient is at most 1.99e-16 for the
structured mesh and 3.59e-16 for the candidate after both solves. These are
roundoff-scale changes in that particular norm, not a forward-error bound.
The independent residual calculation uses A0 rather than assuming restoration is exact.

## Algebraic diagnostics

Residuals were accumulated by sparse row operations in 128-bit BigFloat arithmetic
on the exact saved binary64 coefficients and vectors. There was no high-precision
assembly, solve, dense factorization, adjoint or error estimator. Norms below are
in the original assembled coordinates.

| Mesh | Configuration | Source | Native relative L2 residual | Independent relative L2 residual | Componentwise backward error |
|---|---|---:|---:|---:|---:|
| Structured | Original | 1 | 2.152e-10 | 2.106e-10 | 0.949588 |
| Structured | Original | 2 | 2.401e-10 | 2.164e-10 | 0.965846 |
| Structured | Three corrections | 1 | 9.139e-11 | 7.206e-11 | 0.885604 |
| Structured | Three corrections | 2 | 9.925e-11 | 1.122e-10 | 0.788628 |
| Candidate | Original | 1 | 0.1774794 | 0.1774792 | 0.899103 |
| Candidate | Original | 2 | 0.2838333 | 0.2838333 | 0.939994 |
| Candidate | Three corrections | 1 | 0.3908461 | 0.3908460 | 0.922521 |
| Candidate | Three corrections | 2 | 0.0107598 | 0.0107597 | 0.922675 |

The requested componentwise diagnostic is

\[
\eta_c=\max_i\frac{|b_i-(A_0x)_i|}{\sum_j |(A_0)_{ij}|\,|x_j|+|b_i|}.
\]

Here complex absolute value is the Euclidean modulus. Zero/zero is defined as
zero; nonzero/zero as infinity. No rows were excluded and no denominator threshold
was added. All observed denominators were nonzero. The large maxima are attained
at zero-RHS rows with tiny terms. For example, structured/original source 1 has
numerator 1.785e-18 and denominator 1.880e-18 at its worst row; after correction,
the worst row changes and has numerator 6.143e-28 and denominator 6.937e-28.
Candidate worst-row denominators range from 1.454e-21 to 6.853e-20.
These raw ratios should not be read as percent errors in conductance. MUMPS's
native W1/W2 diagnostics use a different treatment of small-denominator rows;
they are separately retained in the log, not substituted for the requested formula.

## Observed P and Y changes

| Mesh | Original G12 [S/m] | After three corrections [S/m] |
|---|---:|---:|
| Structured | -1.178474438623045e-24 | -1.178474322412431e-24 |
| Candidate | -6.752131874604136e-22 | -2.119778043706841e-21 |

The structured G12 changes by 1.1621e-31 S/m, approximately 9.86e-8 relative.
Candidate G12 changes by -1.444564856e-21 S/m. Candidate G11 changes from
-6.7467e-22 to +7.5992e-21 S/m; G21 changes from -6.7427e-22 to +7.6039e-21 S/m.
These are raw values; nothing was clipped or replaced.

The exact inverse-difference attribution, evaluated with 256-bit arithmetic on the
saved decimal P and Y tables, gives the following contributions to the candidate's
observed change in G12:

\[
\Delta Y=-Y_{\rm original}\,\Delta P\,Y_{\rm refine3}.
\]

| Changed P entry | Contribution to delta G12 [S/m] |
|---|---:|
| P11 | -2.230516141e-21 |
| P21 | +4.706347741e-22 |
| P12 | +3.996396827e-22 |
| P22 | -8.432317256e-23 |
| Sum | -1.444564856e-21 |

The first source column contributes -1.759881366e-21 S/m; the second contributes
+3.153165102e-22 S/m. Thus the output shift is not attributable to the improving
second-source residual alone. Closure error in the real identity is 3.34e-37 S/m,
consistent with saved-table rounding. This is attribution of an observed change,
not an accuracy certificate. The complete complex P values are tabulated below.
All four terminal P matrices have condition number approximately 1.535. Independent
256-bit inversion reproduces every native real Y entry to within 4.57e-16 relative;
the final small-matrix inversion does not account for the observed conductance changes.

## Historical reproduction limit and conclusion

The original-setting reruns did not exactly reproduce the older saved tiny real
outputs. Structured G12 was previously -2.342250689e-25 S/m, versus
-1.178474439e-24 in this diagnostic; the largest P difference is a real P12 shift
of -0.01017687845 ohm m, while its printed imaginary entries are unchanged.
The candidate was previously -7.316754025e-22 S/m, versus -6.752131875e-22 here.
The historical assembled matrices and solution vectors were not retained, so
coefficient and solution identity against those older executions cannot be checked.
This experiment does not resolve the origin of that historical reproduction gap.
Its within-mesh original/refine3 comparison does have verified identical systems.

The analytical G12 retained for this case is -5.953725939e-25 S/m. Proximity to
that value, or to the structured result, is not used as evidence of algebraic
accuracy or as a stopping rule.

For the fixed comparisons, changes in the solution procedure alone materially
change the candidate's extracted output. The large original-coordinate residual
is confirmed independently and is not merely a scaled logging artifact. However,
three corrections do not uniformly improve the candidate, do not yield a known
P_h, and do not isolate the exact discrete-operator difference between meshes.
In the decomposition

\[
\widehat P_c-\widehat P_r=(P_{h,c}-P_{h,r})+(E_c-E_r),
\]

the separate terms remain unresolved. Neither algebraic error nor discretization
error is excluded. No mesh, solver setting, domain size or production default is
promoted. Production and existing scientific tolerances remain unchanged.

Native elapsed times including diagnostic matrix dumps were 65.21/71.81 seconds
for structured original/refine3 and 46.28/52.96 seconds for candidate
original/refine3. These are diagnostic executions, not warmed performance claims.

References: [PETSc diagonal scaling](https://petsc.org/release/manualpages/KSP/KSPSetDiagonalScale/),
[PETSc MUMPS controls](https://petsc.org/release/manualpages/Mat/MATSOLVERMUMPS/),
[KSPPREONLY tolerance semantics](https://petsc.org/release/manualpages/KSP/KSPPREONLY/),
and [MUMPS user guide](https://mumps-solver.org/doc/userguide_5.7.3.pdf).
The installed MUMPS 5.3.3 execution logs supply the effective settings and step counts.

## Complete P comparison

All entries below have units ohm m; rows are observations and columns excitations.
Real and imaginary parts are copied from the native tables without clipping.

| Mesh | Entry | Original real | Corrected real | Original imaginary | Corrected imaginary |
|---|---|---:|---:|---:|---:|
| baseline | 11 | -0.016382875057965884 | -0.016382873037788858 | -108676851694.54036 | -108676851694.54163 |
| baseline | 12 | -0.019213480736333558 | -0.019213478722820757 | -22920945748.428646 | -22920945748.428688 |
| baseline | 21 | -0.0090189738371274641 | -0.0090189718170439193 | -22920945748.428642 | -22920945748.428837 |
| baseline | 22 | -0.016389228045807391 | -0.016389226032371196 | -108625710505.41879 | -108625710505.41879 |
| metric-1.35 | 11 | -11.68470057956192 | 102.29766197495063 | -108676851569.20616 | -108676851844.4269 |
| metric-1.35 | 12 | -11.687466525165636 | -7.3802533171672327 | -22920945623.094723 | -22920945661.563747 |
| metric-1.35 | 21 | -11.676783663520212 | 102.29999977559507 | -22920945623.10062 | -22920945898.308437 |
| metric-1.35 | 22 | -11.684089264845195 | -7.3770899493161268 | -108625710380.09131 | -108625710418.5584 |
