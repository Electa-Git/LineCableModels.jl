# Finite-Gamma sensitivity: metal-drive cancellation

Investigation on 2026-09-30, before releasing localization in both media.
All original bundles, meshes and results remain unchanged under
`.linecablemodels/fem/local-interface-mesh/`. These experiments and acceptance
comparisons belong to the manual harness, not the production engine.

## Reproduction

Two copper wires, radius 0.0425 m, centres (0,-1) and (1,+1) m; soil resistivity
0.1 ohm m; frequency 1 MHz; prescribed complex Gamma = 0.99 gamma_earth.
Use the retained `mixed-f1.0e6-gamma0.99/{baseline,localized,localized-decay}`
bundles. The three meshes have 95963, 95964 and 95965 nodes. Their PML, terminal
contour and native voltage-path coordinate hashes are identical. The physical
wave-size targets are mathematically identical throughout this particular
domain; see `constant-target-bounds.toml`.

The original equations give maximum component changes of 57.50% and 66.74% for
compact and wider localization. The affected quantities belong mainly to the
aerial-source mutual column. This is not caused by admittance inversion: the
potential coefficient already changes before inversion.

## Fixed-mesh solver controls

`finite_gamma_solver.jl` copies the frozen baseline mesh and native sources.
The following changes use exactly the same assembled formulation and mesh:

| Control | Largest R/X/G/B change from original | Observation |
| --- | ---: | --- |
| Rebuild/factor each source | 0 | Factorization reuse is not the cause |
| MUMPS AMF ordering | 4.104% | Sensitive to elimination ordering |
| GMRES with LU | 24.271% | Small recurrent residual does not resolve the issue |
| MUMPS scaling and default refinement stopping criterion | 0 | Zero refinement steps performed |
| Refinement stopping criterion 1e-16 | 29.341% | Residual improves but coefficients remain sensitive |

The original source-2 residual norm is 1.5493e-5 for RHS norm sqrt(2).
Tighter native refinement reduces it to 1.6297e-6. MUMPS condition estimates
reported by the full error analysis are system/RHS-dependent diagnostics, not
an independent certification of any weak terminal observable. The scaled
backward residual alone is insufficient evidence of accurate coefficients.
Relevant native controls are documented in
[PETSc's MUMPS interface](https://petsc.org/release/manualpages/Mat/MATSOLVERMUMPS/)
and the [MUMPS user guide](https://mumps-solver.org/doc/userguide_5.8.2.pdf).

## Exact change of metal unknown

With s = j omega and kappa = sigma + s epsilon, the existing metal branch is

    Ez = -(s a + u - Gamma^2 V)

The scalar potential is one constant V on each terminal. Thus the terms
`kappa u` and `-Gamma^2 kappa V` have exactly the same metal support and constant
basis. Assembling these separate large terms makes the small driven current
depend on their numerical cancellation. Define instead

    w = u - Gamma^2 V
    Ez = -(s a + w)                 in the finite metal
    u = w + Gamma^2 V              for the longitudinal-drive output

This is an invertible linear substitution. It adds no approximation, new
boundary condition, material assumption, mesh parameter, unknown, or solve.
The exterior Gamma-squared terms, scalar continuity, transverse vector
equation, PML and native voltage-path integration remain unchanged. The metal
v terms cancel before assembly in both the magnetic and total-current rows.
At Gamma=0 it is the original unknown exactly.

`finite_gamma_terminal_shift.jl` applies this substitution only to copied
qualification bundles. `qualify_terminal_shift.jl` compares old and corrected
results separately; it never overwrites the original reference or silently
counts a corrected-reference comparison as original-reference preservation.

## Isolated result

On the frozen high-frequency mixed meshes, the exact substitution gives:

| Comparison, with substitution on both sides | Largest component change |
| --- | ---: |
| Original vs compact mesh | 0.0003461% |
| Original vs wider mesh | 0.0014414% |
| Original mesh, default vs AMF ordering | 0.0001717% |

All component signs agree. The source-2 residual on the original mesh falls
to 1.2076e-10. This directly isolates metal-drive cancellation as the dominant
cause of the reported remeshing sensitivity, without changing meshing or
voltage extraction. Native elapsed time is 33.58 s for the corrected original
mesh and 33.62 s for the corrected compact mesh in this control; no extra
factorizations or iterations are requested.

The original weak X12 was -0.00245317 ohm/m; the corrected value is
-0.000204241 ohm/m; analytical is -0.000241417 ohm/m. This correction therefore
does not preserve the numerically unstable original coefficient within 2%.
The 2% localization comparison must hold the equations fixed on both meshes.
Original references are retained explicitly to expose that distinction.
Other analytical discrepancies remain, including the pre-existing sign of
weak X21. Agreement between corrected meshes is not a claim of full analytical
accuracy or universal absence of sign disagreements.

All six finite-Gamma placement/frequency pairs passed before production edits:

| Placement | 0.1 Hz maximum change | 1 MHz maximum change |
| --- | ---: | ---: |
| Mixed | 0.234699% | 0.0003461% |
| Air | 0.227507% | 0.0004273% |
| Soil | 0.159896% | 0.0002537% |

No comparison changed a component sign. The Gamma=0 control gave bitwise
identical Z, Y and P on each retained mesh, with and without the substitution.
Together with the earlier nine Gamma=0 localization comparisons this completes
the prescribed 15-case numerical qualification. Production now uses the exact
metal substitution and both-media localization. Public managed/detached
execution verification also completed. The qualification
records are `terminal-shift-comparison.toml` and `terminal-shift-components.csv`
in each case directory. Single serial log: `local-interface-mesh/live.log`.

## Public construction-route control

The public 1 MHz mixed-case mesh initially differed from the corrected native
qualification mesh by 2.0654% in weak X21 (about 2.6 micro-ohm/m). The two
execution routes agreed on an identical supplied mesh; the discrepancy arose
between independently constructed meshes. It was not counted as a 2% pass.

`managed_reference.jl` then used the captured pre-edit Julia mesh builder in an
isolated qualification module, with the same corrected equations. No production
method was replaced. The original managed mesh already differs from the
original detached mesh by 2.0655%. Managed before/after localization differs
by only 0.0001547%, without new signs. This separates a pre-existing
construction-route discretization difference from the optimization. It does
not establish 2% interchangeability of independently generated meshes.

The high-frequency public preservation check therefore uses that matched
managed reference; the original native comparison and its failure are retained.
The limit remains 2%. `managed-reference-terminal-shift/comparison.toml` records
both differences, and the captured builder is `managed-reference-configure.jl`
under the evidence root. The detached qualification still covers all 15 points.

Public managed/detached/resume checks passed at mixed 0.1 Hz for Gamma=0 and
finite Gamma, and at mixed 1 MHz for finite Gamma (42 assertions). Six maintained
test items passed, including prescribed-Gamma manufactured Maxwell fields,
cylindrical reference, native complex extraction, finite-Gamma integrated
metal currents, map output and interface controls. Production screen/tube/sector
mesh checks preserved all prescribed conductor controls and PML/contour/path
coordinate hashes. No production validation or refinement loop was added.

The first public 0.1 Hz computation took 72.82 s, with 38.25 s reported compilation;
the next finite-Gamma computation took 42.13 s, with 0.028 s compilation. These
are different physical cases, not a speedup comparison. The final 1 MHz check
reused its existing solve: its 17.37 s were almost entirely fresh-process Julia
compilation and must not be reported as a FEM solve time. The matched Gamma=0
native pair took 30.34 -> 28.52 s; the earlier repeat was 32.43 -> 32.60 s.
Mesh/DOF and memory savings are established; a general runtime saving is not.
