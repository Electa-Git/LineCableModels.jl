# Prescribed complex Γ in the coupled FEM formulation

Implemented and checked on 2026-09-29 in the FEM investigation checkout.
The separate manufactured-field experiment was executed before extending the
production GetDP equations. The old `electrodynamic_full_wave.pro` was inspected
for conventions; the implementation retains the current terminal spaces,
native voltage paths, material domains and extraction owner.

## Input and exact normalization

```julia
fem = Formulation(:LineCableModelsFEM;
    options=(physics=:quasi_fw, Γ=0.01+0.02im))
```

Γ has units 1/m. A finite scalar applies at every frequency; a finite vector
follows `problem.frequencies`, including repeated frequencies. Zero is the
default. The physics selector defaults to `quasi_fw`, its supported choice.
Γ belongs to formulation options, and changes invalidate run reuse.

With exp(jωt−Γz), write Az=a, At=Γb, φ=Γv, s=jω, κ=σ+sε and ν=1/μ.
For nonzero Γ this is an exact change of variables. The exterior equations are

```math
-\nabla_t\cdot[\nu(\nabla_t a+\Gamma^2 b)]+s\kappa a-\Gamma^2\kappa v=0,
```

```math
C^*(\nu Cb)+(s\kappa-\Gamma^2\nu)b+\kappa\nabla_t v-\nu\nabla_t a=0,
```

```math
-\nabla_t\cdot[\kappa(sb+\nabla_t v)]+\kappa(sa-\Gamma^2v)=0.
```

Here Cb=∂x by−∂y bx and C*h=(∂y h,−∂x h). The native tree gauge removes
potential redundancy; continuity supplies the nodal part of transverse Ampère.
The transformed PML tensors require the explicit rotations retained in the
GetDP weak form. Γ² is the complex product Γ·Γ.

The exact-zero branch omits only the Γ² terms. Both branches use the same
normalized unknowns, so there is no division by Γ or small-value threshold.

The physical exterior fields are

```math
E_t=-\Gamma(sb+\nabla_t v),\qquad E_z=-sa+\Gamma^2v,
```

```math
B_t=-\hat z\times(\nabla_t a+\Gamma^2b),\qquad B_z=\Gamma Cb.
```

Finite metal retains the existing axial a/u model and transverse equipotential
terminal approximation. Its axial field and total-current constraint include
Γ²v. This extension does not add a transverse field solve inside finite metal.
The exterior Maxwell equations are complete under these conductor assumptions.

## Native extraction and maps

The existing reference-to-terminal measurement path gives

```math
P_{ij}=\frac{v_i-v_{R_i}+s\int_{\ell_i}b\cdot d\ell}{I_j},\quad
K_{ij}=-U_i/I_j,\quad Z_{ij}=K_{ij}+\Gamma^2P_{ij},\quad Y=P^{-1}.
```

For a PEC contour, K=(sA−Γ²V)/I. The Γ²P correction is applied in GetDP after
the physical voltage is assembled, before matrix reductions. P has units Ω m;
Z has units Ω/m. The corresponding analytical coefficient is Pe=sP.

Maps include the complete B, Ez and axial current, plus transverse physical
loss at finite Γ. The regular transverse maps retain Et/Γ and Jt/Γ, with
explicit labels; multiply their complex values by Γ to recover physical fields.
PML map values are analytic continuation.

## PML and export

The transverse scales are qm=√(sμmκm−Γ²), using the decaying root and the
positive-imaginary root on the imaginary axis. Physical wave-size controls,
PML strengths and the prescribed physical strip distribution consume these
scales for each frequency/Γ pair.

The cubic stretch is S=1+(1−jη)b(d/t)³. One common η is used in both media to
preserve interface matching. It starts at one; if Im(qm)<0, it is bounded by
Re(qm)/(2|Im(qm)|), ensuring positive decay along the complex ray. Side
strength uses both media. The original Γ=0 strengths and meshes are preserved.

For finite Γ, strengths use rates
rm=Re[qm(1−jη)] Im(qm0)/(Re(qm0)+Im(qm0)), where qm0 is the Γ=0 root.
This retains the original normal-wave attenuation convention. At an exact
transverse cutoff, the base stretch stays finite; exponential absorption of
the zero transverse mode is impossible, and algebraic truncation requires
its own domain/refinement check.

Managed jobs receive ΓRe, ΓIm and the stretch slope explicitly. Detached
ONELAB exports store their frequency-aligned arrays and select all coefficients
together. Re-export after changing Γ to keep the equations, PML and mesh
consistent. Directly changing a solver coefficient on an old mesh does not
constitute a new mesh/PML qualification.

## Independent references and observed results

Maintained tests are in
[`fem_prescribed_gamma.jl`](../../extensions/fem_prescribed_gamma.jl), with
fixtures in [`prescribed_gamma`](../../fixtures/data/fem/prescribed_gamma/).
All five new test items passed, totaling 199 assertions across focused runs.

| Check | Observed result |
| --- | --- |
| Input ownership, frequency transport, PML ray and physical strips | 36 assertions passed |
| Finite-metal solve, repeated-frequency Γ vector, resume and detached native matrices | 9 assertions passed |
| Manufactured Maxwell fields, independent fixture and production equation block | 130 assertions passed |
| Finite cylindrical exterior, exact Bessel reference | 12 assertions passed |
| Infinite cylindrical exterior with production Cartesian PML | 12 assertions passed |

The manufactured problem exercises Γ=0, ±(2+j), 10⁻⁶(2+j), 2j and 2, both in
isotropic material and with complex anisotropic coordinate tensors. Gauge-
invariant fields converge at their expected finite-element orders. It checks
complex-square signs, Γ sign symmetry and the regular zero limit.

For a PEC cylinder of radius 0.2, μ=1.25, κ=0.7+0.68j and ω=1.7, the finite
outer radius is 2. The refined h=0.08 errors against the exact annular Bessel
solution are:

| Γ | Relative Z error | Relative P error |
| --- | ---: | ---: |
| 0 | 0.04585% | 0.04585% |
| 0.99√(sμκ) | 0.04459% | 0.04459% |
| 0.6+0.3j | 0.05824% | 0.05824% |

For the infinite-domain reference,
a/I=μ K0(qr)/(2πrqK1(qr)), Z=sa/I and P=a/(μκI).
The Cartesian PML starts at |x| or |y|=1 and ends at 3. Its entire variable
tensor profile and the voltage path through it are exercised. Refined errors
are 0.04198%, 0.56868% and 0.05426% for Γ=0, 0.99√(sμκ) and 2+2j, respectively,
for both Z and P. The last case needs a reduced ray slope. Halving h reduces
the errors by factors of 3.1–3.8. This economical fixture uses a 1% refined
acceptance bound; the configured continuous reflection is not a mesh-error bound.

A saved production two-conductor mesh was also solved with the original and
extended equations: Γ=0 Z and P were exactly identical. With 10⁻⁶γearth,
relative changes were 1.80×10⁻¹² in Z and 1.05×10⁻¹⁴ in P; ±Γ gave identical
matrices. This isolates the operator on a fixed mesh. Its finite-Γ numbers
with the original PML are not a qualification of the new open-domain problem.
A separate finite-Γ map check verified Bz=Γμ(Hz/Γ) and Jz=κEz in physical air.

Existing regression checks also passed: physical-PML strip preservation
(338 assertions), native complex extraction (28), physics/maps/resume (35),
detached numerical parity/relocation (68), resume ownership (85), option
contracts (174), and the export ownership/geometry checks (53).

Two existing assertions in `fem_export.jl:173` still fail because a frequency
scan and a separate solve do not produce byte-identical table text. The same
two failures were reproduced using the saved original field equations and
original PML ray with the current geometry/export path (42 pass, 2 fail).
That comparison isolates this equation extension; it is not a pristine-checkout
baseline. The strict assertions were left unchanged. The observed matrix
differences are at numerical roundoff/solver-accuracy scale (Z differences
around 10⁻¹⁰ Ω/m). The wider suite is therefore not reported as entirely green.

Commands used for the maintained checks include:

```sh
julia --project=test --startup-file=no test/runtests.jl fem_prescribed_gamma
julia --project=test --startup-file=no test/runtests.jl 'physical PML preserves'
```

The new items were also run individually after corrections. Native ONELAB scan
checks require local Unix-socket access. Transient logs and isolated experiments
are under `/tmp/lcm-gamma-study/`; maintained tests reconstruct their own data.

These checks establish the equations, extraction, regular zero limit and PML
continuation. They do not establish arbitrary-geometry agreement with the
analytical circumferential-averaging model, nor qualify every finite-Γ
low-frequency conductance sign with the default mesh controls.
