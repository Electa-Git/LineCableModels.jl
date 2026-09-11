# Manual coupled quasi-full experiment

`ext/LineCableModelsGmshExt/getdp/model.pro` selects `quasi-full.pro` through
the ONELAB constant `Physics=1`. The public backend now also exposes this model
as `options=(physics=:quasi_fw,)`; its default remains `:quasi_tem`.
The manual experiment below keeps the explicit PEC conductor branch.

From the repository root, run in the Julia REPL:

```julia
include("dev/run_quasi_full.jl")
```

The script builds two bare wires with radius 0.0425 m, separation 1 m and depth
1 m, in earth with resistivity 0.1 Ω m and relative permittivity/permeability 1.
The material object has resistivity 10⁻¹² Ω m, but the script explicitly passes
`-setnumber PerfectConductors 1`: both field domains exclude metal interiors,
and the axial problem uses exact PEC contour current constraints. There are no
insulation regions or internal impedance contributions in this manual run.

The script constructs `system`, `problem` and `formulation`, uses the existing
material adapter and Gmsh mesher, then invokes `model.pro -setnumber Physics 1` through
the GetDP CLI. It performs no production `compute` call. It prints each command,
native solver progress and the complete Z, Pe and Y matrices. A new run directory
under `.linecablemodels/quasi-full/runs/` retains the meshes, material data,
solver-source snapshot, voltage paths and native output columns.

The following variables remain available in the REPL, indexed by response,
source and frequency:

| Variable | Meaning | Units |
|---|---|---|
| `Zqf` | Earth/exterior series impedance with exact PEC contours | Ω/m |
| `Mqf` | Inverse shunt admittance; GetDP's raw `P` output | Ω m |
| `Peqf` | Potential coefficients, `jω Mqf` | m/F |
| `Yqf` | Shunt admittance, `inv(Mqf)` | S/m |

Edit the inputs at the top of the script to change frequency or material.
Set `plot_field_maps=true` for native Gmsh field maps. The path helper currently
delegates to the backend's mesh integration, using the specified circular
electrode endpoints. The public backend derives endpoints from terminal contours
for general geometry. Paths may cross equipotential metal, where the transverse
field is zero; leaving the field mesh is rejected. The `.pro` accepts general
mesh-coordinate integration points and oriented line weights.

The `.pro` defaults to `PerfectConductors=0`, retaining the existing finite-metal
`a/u` branch. The script's explicit PEC setting is deliberate: at 0.1 Hz and
metal resistivity 10⁻¹² Ω m, that retained branch returns negative self resistance
on this mesh. A fresh call to the unchanged production `model.pro` reproduces
the same value, −1.47888×10⁻⁷ Ω/m. This is a failure in the inherited magnetic
calculation for those inputs, separate from the transverse-field repair; its
underlying cause is not established here. Exact PEC gives positive resistance
without approximating a perfect conductor by an extreme conductivity.

## Equations and extraction

Use `exp(jωt − Γz)`, κ = σ + jωε and ν = 1/μ. The solved variables are
`a = A_z`, `bt = A_t/Γ`, `v = φ/Γ`, and the existing axial current-port unknowns
`u`, `U`, `I` in the finite-metal branch. The PEC branch instead groups the
constant axial-potential trace `A` and its associated current `I` on each
electrode. In the surrounding media the equations are

```math
-\nabla_t\cdot(\nu\nabla_t a)+j\omega\kappa a=0,
```

```math
C^*(\nu C\mathbf b)+j\omega\kappa\mathbf b
+\kappa\nabla_t v-\nu\nabla_t a=0,
```

```math
-\nabla_t\cdot[\kappa(j\omega\mathbf b+\nabla_t v)]
+j\omega\kappa a=0.
```

Here `C b = ∂x b_y − ∂y b_x` and `C* h = (∂y h, −∂x h)`. The axial block
retains the working `a/u` equations inside the metal when enabled. For PEC,
`Z = jω A/I` is extracted from the exterior magnetic solve. The normalized electric
terminal source is the same global axial current `I`; there is no second,
independently prescribed electric excitation.

This is the analytically normalized Γ → 0 formulation. It retains the first
order transverse fields and omits second order feedback into the axial block.
It does not solve for a propagation eigenvalue or implement finite-Γ modes.

The transverse field uses `Form1` / `BF_Edge`. A tree gauge is completed on
**all** boundaries with prescribed transverse-potential circulation: the outer
boundary and all electrode contours. Omitting the electrode contours from the
tree's starting boundary can remove physical loop degrees of freedom. The
native regression changes both node and element ordering to obtain a different
tree and checks that physical voltage remains invariant while scalar traces
change.

The voltage is evaluated along a reference-to-electrode path:

```math
\frac{V_i}{\Gamma}=v_i-v_{\rm ref}
+j\omega\int_{\rm ref}^{i}\mathbf b\cdot d\boldsymbol\ell,
\qquad M_{ij}=\frac{V_i^{(j)}}{\Gamma I_j}.
```

The manual paths run from earth infinity to the bottom of each electrode.
Their finite-domain portion is vertical. The shell portion is the pullback of
that physical vertical ray through GetDP's `VolSphShell` transformation,
approximated by 128 straight segments. Each segment is cut at triangle
boundaries; one midpoint integrates the tangential component of a lowest-order
edge field exactly on each straight piece. `Jacobian Plain` evaluates the
pulled-back one-form, preserving its circulation through the mapped shell.
The native postoperation accumulates this integral before writing `P.tsv`.
`Pscalar.tsv` separately records the gauge-dependent scalar contribution.

The tree mechanism and edge space are documented in the
[GetDP manual](https://getdp.info/doc/texinfo/getdp.html#Edge-finite-element-space-with-gauge-condition).
The shell mapping is defined in
[GetDP 3.5.0 Get_Geometry.cpp](https://github.com/getdp-project/getdp/blob/getdp_3_5_0/src/kernel/Get_Geometry.cpp).
The preceding operator diagnosis is in
[fem-coupled-electric-extraction.md](fem-coupled-electric-extraction.md).

## Native validation, 2026-09-11

The corrected manual script completed 0.1, 1, 10, 100, 1,000, 10,000, 100,000
and 1,000,000 Hz: eight GetDP processes and sixteen current-excitation columns,
with one factorization per frequency. These are fresh solves, with no fixture
or analytical matrix supplying any FEM output.

The 39 checks in `test/extensions/fem_quasi_full.jl` pass. They cover exact
edge-field circulation, shared-edge counting, path orientation, invalid paths,
native gauge invariance on reordered meshes, factorization reuse, both conductor
treatments and finite matrix output. The eight-frequency PEC run also has
positive-definite real parts of Z and Y. Additional native controls show:

- Uniform-earth coupled and independent scalar Helmholtz matrices agree to
  6.95×10⁻¹⁵ relative error with exact PEC contours.
- Increasing shell-path resolution from 128 to 512 segments changes M by at
  most 3.21×10⁻¹² relative error at 0.1 Hz, 100 kHz and 1 MHz.
- At 1 MHz, Y11 = 40.9005067429 + j19.6572058671 S/m and
  Y21 = −0.0205614319750 − j0.0191033680397 S/m.

Full matrices and native evidence are under
`.linecablemodels/qa/quasi-full/validated-matrices.tsv`,
`validated-run.txt`, `validation-pec.log`, `regression-tests.log` and
`manual-pec.log`. Early exploratory runs with an incomplete tree boundary
are superseded by the run identified in `validated-run.txt`.

These checks establish the implemented formulation and extraction identities;
they do not establish mesh convergence for every matrix entry or geometry.
