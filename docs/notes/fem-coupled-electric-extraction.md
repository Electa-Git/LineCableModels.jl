# Why the original coupled FEM gives a different shunt admittance

2026-09-11. This investigation leaves the working scalar electrodynamic backend
unchanged. It diagnoses the earlier `A_z/u/phi` equations retained with native
run `run-XtjQUt`.

**The original field ansatz omits the transverse Maxwell–Ampère equations.**
Current continuity does not supply those missing vector equations. In a uniform
surrounding medium the omission can be harmless for this extraction. At an
air–earth interface it permits an additional harmonic potential contribution,
so `Phi/(Gamma I)` does not recover the scalar Helmholtz potential coefficients.
The error survives the small-Gamma limit used to extract the electric response.

The earlier explanation that the extraction was merely "mesh-sensitive" was
insufficient. Mesh error exists, but the controlled tests below separate it from
this continuum-operator mismatch. Neither division by a small number nor the
wire's interior skin-depth discretization explains the main discrepancy.

## Equations actually assembled

Use exp(jωt − Γz), κ = σ + jωε, ν = 1/μ, and p = φ/Γ. Outside the wires,
the original formulation tends, as Γ → 0, to

```math
-\nabla_t\cdot(\nu\nabla_t A_z)+j\omega\kappa A_z=0,
\qquad
-\nabla_t\cdot(\kappa\nabla_t p)+j\omega\kappa A_z=0.
```

The first equation is axial Ampère; the second is current continuity after
normalization. The source conductor supplies the same unit-current boundary
reaction to both equations. The finite-conductivity version supplies the second
reaction through its integrated axial current; eliminating that volume source
in favor of the imposed terminal current produces the same numerical result.

Subtracting the equations gives

```math
\nabla_t\cdot\mathbf r=0,
\qquad \mathbf r=\kappa\nabla_t p-\nu\nabla_t A_z.
```

But the original ansatz has only A = A_z ẑ and E_t = −Γ∇t p. It therefore gives

```math
\mathbf H_t=\nu(\partial_y A_z,-\partial_x A_z),
\qquad
(\nabla\times\mathbf H-\kappa\mathbf E)_t
=\Gamma\bigl(\kappa\nabla_t p-\nu\nabla_t A_z\bigr)=\Gamma\mathbf r.
```

**Transverse Ampère requires r = 0. The implemented scalar equation only
requires div r = 0.** A vector field can have zero divergence without being zero.
The missing condition is consequential even when every assembled equation is
solved accurately.

This restriction was already present in the supplied `fullwave.pro`. Its
`BF_PerpendicularEdge` basis represents only `(0,0,A_z)`; that is the basis
documented in the [GetDP manual](https://getdp.info/doc/texinfo/getdp.html).
There are no transverse vector-potential unknowns or corresponding vector test
functions. The full A–φ formulation requires the vector Maxwell–Ampère equation
alongside continuity; see equations 10–11 of
[Ciuprina and Sabriego (2024)](https://link.springer.com/article/10.1186/s13362-024-00165-6).
The component reduction and residual above are derived directly from the
retained GetDP equations, rather than inferred from that paper's benchmarks.

## Why the interface exposes the omission

In a homogeneous medium with constant μ and κ and matching outer references,

```math
p=\frac{A_z}{\mu\kappa}
```

satisfies both equations and both unit-current terminal constraints. Substituting
it yields the scalar electrodynamic equation

```math
\nabla_t^2p-q^2p=0,\qquad q^2=j\omega\mu\kappa.
```

That identity also holds algebraically on a common finite-element mesh when
both fields have the same PEC terminal traces. It does not need a refined mesh.

Across an interface, the original weak forms impose continuity of A_z and p,
and of ν∂n A_z and κ∂n p, respectively. However, A_z/(μκ) is discontinuous when
μκ changes. The original coupled solution must then have the form

```math
p_m=\frac{A_{z,m}}{\mu_m\kappa_m}+h_m,
\qquad \nabla_t^2h_m=0
```

within each homogeneous medium. The harmonic terms enforce the interface and
terminal conditions. They do not solve the scalar Helmholtz equation:

```math
(\nabla_t^2-q_m^2)p_m=-q_m^2h_m.
```

Equivalently the original p satisfies
`∇t²(∇t² − q²)p = 0` in each source-free homogeneous region, allowing a Laplace
component in addition to the desired Helmholtz component. For a transverse
Fourier mode λ, the extra component decays with `exp(−|λ|d)`; the Helmholtz
component decays with `exp(−sqrt(λ²+q²)d)`. Their range and phase differ. This is
particularly significant for a small mutual admittance.

The incompatibility can also be seen directly at the interface. If r = 0 on
both sides, tangential differentiation would require

```math
\partial_xp=\frac{\partial_x A_z}{\mu_0\kappa_0}
=\frac{\partial_x A_z}{\mu_g\kappa_g}.
```

Both p and A_z have continuous traces, but μκ differs between air and earth and
the localized wire field varies along the interface. These equalities cannot
generally hold. Additional transverse field degrees of freedom are required.

## Why a tiny Gamma does not cure it

A consistent expansion generally contains transverse vector-potential terms
at order Γ, say A_t = Γb + higher-order terms. Then

```math
\frac{\mathbf E_t}{\Gamma}
=-\nabla_t p-j\omega\mathbf b+\text{higher-order terms}.
```

The original ansatz sets b = 0. Making Γ smaller reduces both the retained and
omitted transverse fields by the same factor. Extracting P by dividing by Γ
preserves the leading-order error. The corresponding feedback into the axial
magnetic problem enters at order Γ², explaining why the series impedance can
look accurate while the shunt extraction disagrees.

This is a consistency issue in the field reduction. The tests do not implicate
the arithmetic operation of dividing by Γ or the later inversion of P.

## Controlled native results

The experiments reuse the exact same two-wire geometry and meshes. In the
homogeneous control, both upper and lower media have the earth's conductivity,
permittivity and permeability, including the appropriate conductivity
integration domains. PEC controls exclude metal interiors from both field
integrals and impose constant magnetic and electric terminal traces. Every
entry below was computed by native GetDP during this investigation.

At 100 kHz, on the original coarse mesh:

| Surrounding medium; exact PEC wires | Coupled Y21 [S/m] | Scalar Y21 [S/m] |
|---|---:|---:|
| Uniform earth | 0.129028668129 + j1.117124233529 | 0.129028668128 + j1.117124233529 |
| Air over earth | 0.364158478297 + j1.660230995941 | 0.171390183016 + j1.064980280145 |

On the refined production meshes, still with exact PEC wires:

| Frequency | Coupled Y21 [S/m] | Scalar Y21 [S/m] |
|---|---:|---:|
| 100 kHz, 29,465 nodes | 0.356772418273 + j1.652994471935 | 0.171903721782 + j1.059775191905 |
| 1 MHz, 93,391 nodes | −0.028424037228 − j0.020856469229 | −0.020563230005 − j0.019102149439 |

The original finite-copper result at 100 kHz on that refined mesh is
0.356773258728 + j1.652993509139 S/m. Replacing the metal by exact PEC thus leaves
the large mismatch intact. Eliminating the conductor-volume load by the known
terminal current changes results only around 10⁻⁸ S/m in these controls.

The homogeneous PEC identities agree to roundoff on both resolutions, even
though each common discretized result can still differ from the exact circular
continuum solution. The interface mismatch persists under refinement. These
are operator comparisons on shared meshes, not claims of a fully converged
continuum value for the original layered operator.

The first-order transverse Ampère residual was also evaluated directly from
the solved PEC fields. Define, separately in each finite physical medium,

```math
\eta_m=\frac{\|\kappa\nabla_t p-\nu\nabla_t A_z\|_{L^2(\Omega_m)}}
{\|\nu\nabla_t A_z\|_{L^2(\Omega_m)}}.
```

For source 1 at 100 kHz:

| Medium configuration | Mesh | η in lower medium | η in upper medium |
|---|---|---:|---:|
| Uniform earth | Coarse | 1.01×10⁻¹⁴ | 8.61×10⁻¹⁴ |
| Uniform earth | Refined | 1.36×10⁻¹⁴ | 1.21×10⁻¹³ |
| Air over earth | Coarse | 0.04237 | 1.000008 |
| Air over earth | Refined | 0.04212 | 1.000008 |

These are residuals divided by Γ, appropriate to the order extracted for P.
They are not percentages of the total Maxwell residual at exactly Γ = 0.
The upper-medium denominator is small; the ratio identifies a missing balance,
not a claim of 100% error in a terminal matrix. No field symmetrization or fitted
correction was applied.

## Mesh-free interface identity

An independent Fourier calculation gives the same extra harmonic component.
For equal permeability on both sides, define
`a_m = sqrt(λ²+jωμκ_m)`. For a line source of unit axial current at depth d_s,
the magnetic interface coefficient is

```math
B=\frac{\mu}{a_g+a_0}e^{-a_gd_s}.
```

Matching the original p and κ∂n p across the interface gives, at target depth d_t,

```math
\widehat h_g(d_t)=
B\frac{\kappa_g-\kappa_0}{\mu\kappa_g(\kappa_g+\kappa_0)}e^{-|\lambda|d_t}.
```

This vanishes when the media coincide. Otherwise it survives at Γ = 0 and is
not generally invariant under exchanging d_s and d_t. The direct scalar Green
function has the reciprocal reflection term

```math
R_p=\frac{\kappa_g a_g-\kappa_0a_0}{\kappa_g a_g+\kappa_0a_0},
\qquad
\widehat p_{\rm scalar}
=\frac{e^{-a_g|d_t-d_s|}+R_p e^{-a_g(d_t+d_s)}}{2\kappa_g a_g}.
```

The diagnostic script verifies the interface traces, fluxes, homogeneous
identity, and scalar reciprocity directly, without quadrature or a mesh. This
is a line-source operator counterexample; it is not substituted for the
finite-radius FEM matrices above. Unequal-depth nonreciprocity in this
counterexample is distinct from the numerical asymmetry of the equal-depth
two-wire mesh.

## Evidence and scope

Reproducible scripts and immutable source/output directories are under
`.linecablemodels/qa/fem-electric-diagnosis/root-cause/`:

- `audit.jl`, `*-coarse-v2/`, `*-fine-v2/`: material, PEC and volume-load controls;
- `residual.jl`, `ampere-residual-*/`: native field residuals;
- `fourier_identity.jl`: exact Fourier coefficients and interface checks;
- `initial-control-correction.md`: records an invalid exploratory homogeneous
  control and why only the corrected `v2` controls are used as evidence.

The working scalar electrodynamic backend remains in place. This investigation
identifies the missing equations in its predecessor; it does not implement a
complete full-vector finite-Gamma FEM model. Such a model would also need
consistent gauge, terminal-voltage definitions, interface conditions and
outer-boundary treatment.

A separate, normalized first-order implementation is now available for manual
use as `quasi-full.pro`; see
[the manual experiment and native validation](fem-quasi-full-manual.md).
It restores the transverse vector equations and path-voltage extraction while
leaving the registered backend unchanged.
