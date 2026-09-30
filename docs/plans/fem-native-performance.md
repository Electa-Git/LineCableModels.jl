# Native FEM performance execution — 2026-09-30

**Superseded delivery policy:** the researcher subsequently instructed that
numerical choices be implemented for evaluation regardless of screening
thresholds. Harmonic terms are consolidated; physical quadrature and PML
quadrangles are selectable through existing computation/export options, with
existing discretization defaults preserved. The implementation and plot review
are documented in `test/manual/calculations/fem_native_performance/README.md`.
The historical measurements and plan below remain as provenance, not release gates.

Authorized scope: optimize the existing quasi-fw implementation using native
Gmsh/GetDP facilities. No formulation rewrite, changed physical domain, PML
stretch, conductor grading, voltage paths, clipping, or scientific acceptance
inside production. Keep the user's manual runner and previous evidence intact.
Both Julia-managed execution and detached ONELAB must remain supported.

## Fixed sequence

1. Freeze current sources and copy retained three pilot bundles/meshes into
   `.linecablemodels/fem/native-performance-20260930`. Record hashes. Run the
   current sources on these meshes for matched references.
2. On copies only, separately qualify (a) harmonic conductivity/displacement
   assembly consolidation, (b) native region-selected three-point quadrature
   in constant-material linear physical triangles, and (c) structured PML
   quadrangles with 9 and 16 Gauss–Legendre points. Keep the physical mesh and
   all shared boundaries unchanged. No other mesh candidates are authorized.
3. Compare every raw R/X/G/B entry, with at most 2% relative change per nonzero
   entry, exact zeros preserved and no new sign reversals. Report analytical
   values alongside these comparisons. Failed candidates do not expand into
   fixture sweeps and do not enter feature code.
4. Repeat compatible reference/survivor runs twice after the initial run.
   Report meshing, assembly, solve, total native elapsed and peak RSS, plus
   Julia first-use compilation separately. Require a 20% median total saving
   across the pilot workload before delivering a new PML element choice.
   Equivalent assembly changes must preserve results and show measured benefit.
5. Qualify surviving numerical changes on the retained 15 air/buried/mixed
   frequency–Gamma cases, then the three-conductor and screen/tube/sector
   low/high endpoints. Reuse compatible saved references; never regenerate
   completed compatible work. Combine successful changes and check interaction.
   Expose `physical_volume_quadrature=3` separately from the existing
   `volume_quadrature=12` PML rule. Both remain explicit prescribed controls;
   users can retain a higher physical rule without changing PML resolution.
6. Add optional native MUMPS ordering and PETSc preallocation controls without
   changing defaults. Verify managed command and detached UI propagation and
   native execution. Do not claim that a particular ordering is universally best.
7. Qualify binary MSH write/read preservation and implement it if beneficial.
   Audit existing mesh reuse; avoid adding a second cache or changing the
   user's explicit remesh request. Finish focused regression tests and docs.

The manual Julia harness owns all comparisons and decisions. Execution is
serial and resumable. One live log is `native-performance-20260930/live.log`;
native solver output is flushed there. No approvals between trials. A failed
candidate ends that candidate, not the other independent optimizations.

References: [GetDP integration objects](https://getdp.info/doc/texinfo/getdp.html#Integration),
[Gmsh structured/recombined meshes](https://gmsh.info/doc/texinfo/gmsh.html#t6),
[PETSc MUMPS controls](https://petsc.org/release/manualpages/Mat/MATSOLVERMUMPS/).
Installed GetDP 3.5.0 source confirms zero-based `Integration Criterion`
selection and tensor-product Gauss–Legendre quadrangle rules.

## Qualification outcome

The PML quadrangles failed the pilot preservation limit. The combined harmonic
assembly/physical quadrature change passed all 15 bare-wire cases but failed
the screen and low-frequency tube checks. Each assembly change also failed
independently on the low-frequency screen. Consequently step 5 does **not**
authorize promotion: no physical-quadrature option, consolidated weak terms or
quadrangle PML enters feature code. The measured 5.3% combined pilot speedup
belongs to a rejected candidate, not the delivered implementation.

Binary MSH and explicit native solver controls are the independently qualified
deliverables. See the [execution report](../../test/manual/calculations/fem_native_performance/results.md)
for costs, preservation evidence and verification limits.
