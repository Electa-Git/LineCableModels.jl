# Finite element solver

```@meta
CurrentModule = LineCableModels
```

The Gmsh extension computes the series impedances and shunt admittances of a problem by means
of the finite-element method. The cross-section of the system is built and meshed by Gmsh,
and each frequency is solved by GetDP. The matrices of the terminals are then reduced in the
same way as in the coaxial computation.

## Loading the extension

The extension is loaded by `using Gmsh`, and the Julia package Gmsh.jl provides the Gmsh
library with it. The field equations are solved by GetDP 3.5.0 in complex arithmetic, an
external executable that the extension finds or installs as described below. The
finite-element computation is tested on Linux.

## GetDP executable

The computation option `getdp_executable` selects the executable by a
[`GetdpBinary`](@ref Main.GmshExt.GetdpBinary) choice or by a path:

| Choice | Executable |
| --- | --- |
| `:discover`, the default | the `getdp` found on `PATH` when its output of `getdp -info` reports version 3.5 and complex arithmetic, and otherwise the executable of `:download` |
| `:download` | the GetDP 3.5.0 complex build installed as an artifact of the package |
| a path, given as a string | the executable at that path, used as is |

Any other symbol throws an `ArgumentError`. With `:download`, or with `:discover` when no
suitable GetDP is found on `PATH`, the extension installs the lazy artifact `getdp` of the
`Artifacts.toml` of the package at its first use. Pkg then downloads the archive of the
platform and verifies its SHA-256 hash and its tree hash. The archive is unpacked under the
`artifacts/` folder of the Julia depot. Only this first use requires network access: later computations
run the installed copy, which stays in the depot until `Pkg.gc()` removes it. Since
`:discover` looks at `PATH` before every computation, a suitable GetDP installed later is
used in place of the artifact.

| Platform | Archive |
| --- | --- |
| Linux on x86-64, glibc | `getdp-3.5.0-Linux64c.tgz` of getdp.info |
| macOS on x86-64 and on Apple silicon | `getdp-3.5.0-MacOSXc.tgz` of getdp.info, an x86-64 build that Apple silicon runs through Rosetta 2 |
| Windows on x86-64 | `getdp-3.5.0-Windows64c.tar.gz`, the files of the getdp.info build repacked as a tarball in the release `getdp-3.5.0` of the package |

On another platform, or on a machine without network access, GetDP 3.5 with complex
arithmetic is installed separately and given on `PATH` or by its path.

Before each computation, the extension runs the selected executable with `-info` and
requires an output that names GetDP and a version number. Otherwise, it throws an error. The
path of the executable must be shorter than 256 bytes, a limit of the native path
arguments. The file `logs/getdp-info.log` of the run directory records the executable, the
choice and the output of `-info`.

## Finite-element formulation

The [`FEM`](@ref) formulation selects the material laws of a finite-element computation. Its
keywords are:

| Keyword | Default | Meaning |
| --- | --- | --- |
| `insulation_admittance` | `:lossless` | admittivity law for insulator layers |
| `semicon_admittance` | `:lossless` | admittivity law for semicon layers |
| `resistivity` | `:linear` | conductor resistivity law |
| `earth_properties` | `:constant` | air and earth properties law |
| `Γ` | `0` | longitudinal wavenumber in 1/m, one value or one per frequency |
| `bundle` | `true` | merge the conductors that share a positive connection index |
| `kron` | `true` | remove the conductors of connection index zero |
| `transposition` | `false` | apply the ideal cyclic transposition |

Since the trait [`backend`](@ref) of [`FEM`](@ref) is [`GmshBackend`](@ref), the call
`compute(problem, FEM())` runs the extension.

## Field equations

The finite element model solves the fields of the cross-section with the dependence
``e^{\mathrm{j}\omega t - \Gamma z}``, where ``z`` is the axial coordinate, ``x`` and ``y``
are the horizontal and vertical coordinates of the cross-section, and the interface between
air and earth lies at ``y = 0``. Each region has the complex conductivity
``\widehat\sigma = \sigma + \mathrm{j}\omega\varepsilon`` and the reluctivity
``\nu = 1/\mu``. A conductor takes the admittivity
``\kappa = 1/\rho(\theta) + \mathrm{j}\omega\varepsilon_0\varepsilon_r``, with ``\rho(\theta)``
from the resistivity law at the operating temperature, and a dielectric the admittivity of its
law. Both then have ``\sigma = \operatorname{Re}\kappa`` and
``\varepsilon = \operatorname{Im}\kappa/\omega``. The air and the earth take
``\sigma = 1/\rho``, ``\varepsilon = \varepsilon_0\varepsilon_r`` and ``\mu = \mu_0\mu_r`` from
the law of the earth properties.

The fields are expressed by the normalized potentials ``a``, ``\mathbf b`` and ``v``, with
``A_z = a``, ``\mathbf A_t = \Gamma\mathbf b`` and ``\phi = \Gamma v``, in the framework of
the potential formulation of [Ciuprina2024](@cite). The normalization retains every
``\Gamma^2`` term and has a regular limit for ``\Gamma = 0``:

```math
\mathbf E_t = -\Gamma\left(\nabla_t v + \mathrm{j}\omega\mathbf b\right), \qquad
E_z = \Gamma^2 v - \mathrm{j}\omega a, \qquad
\mathbf B_t = \nabla_t \times (a\hat{\mathbf z}) - \Gamma^2\hat{\mathbf z} \times \mathbf b.
```

Outside the metal, the axial component of Ampère's law, its transverse component divided by
``\Gamma``, and the continuity of the current divided by ``\Gamma`` give

```math
\begin{aligned}
\hat{\mathbf z} \cdot \nabla_t \times (\nu\mathbf B_t) &= \widehat\sigma E_z, \\
\nabla_t \times (\nu\nabla_t \times \mathbf b)
    + (\mathrm{j}\omega\widehat\sigma - \Gamma^2\nu)\mathbf b
    + \widehat\sigma\nabla_t v - \nu\nabla_t a &= 0, \\
\nabla_t \cdot \left[\widehat\sigma\left(\nabla_t v + \mathrm{j}\omega\mathbf b\right)\right]
    + \widehat\sigma\left(\Gamma^2 v - \mathrm{j}\omega a\right) &= 0.
\end{aligned}
```

In the metal of terminal ``k``, the axial field is ``E_z = -\mathrm{j}\omega a - w_k``, with
one unknown ``w_k`` per terminal, and the axial equation holds with the total current

```math
\int_{\Omega_k} \widehat\sigma E_z \, \mathrm{d}\Omega = I_k.
```

The transverse fields are solved outside the metal only. The contour of each terminal is an
equipotential electrode at the floating potential ``V_k``, with ``\mathbf b \times \mathbf n
= 0``, and its outward leakage current per unit length is ``\Gamma I_k``, which follows from
the continuity of the axial current. The physical axial drive of the terminal is
``w_k + \Gamma^2 V_k``. The outer boundary of the perfectly matched layer has ``a = 0``,
``\mathbf b \times \mathbf n = 0`` and ``v = 0``, and a tree gauge, started on the conductor
contours and on the outer boundary, removes the gradient part of ``\mathbf b``. GetDP solves
the Galerkin form of these equations with first-order nodal functions for ``a`` and ``v``
and edge functions for ``\mathbf b``, by the direct solver MUMPS. Each column ``j`` of the
solution imposes ``I_j = 1`` A on terminal ``j`` and zero on the other terminals, in both
the axial current and the leakage.

## Line parameters from the fields

The voltage of terminal ``i`` divided by ``\Gamma`` is the integral of the normalized
transverse field along a vertical line ``\ell_i``, from the interface for a terminal in the
air, or from the outer boundary of the bottom layer for a terminal in the earth, to the metal
of terminal ``i``. With the excitation of column ``j``,

```math
P_{ij} = \frac{1}{I_j} \int_{\ell_i} \left(\nabla_t v + \mathrm{j}\omega\mathbf b\right)
    \cdot \mathrm{d}\boldsymbol\ell, \qquad
K_{ij} = -\frac{w_i + \Gamma^2 V_i}{I_j}, \qquad
Z_{ij} = K_{ij} + \Gamma^2 P_{ij}, \qquad
\mathbf Y = \mathbf P^{-1},
```

where ``w_i`` and ``V_i`` are the axial drive and the floating potential of terminal ``i`` in
column ``j``, ``P`` is in Ω·m, ``K`` and ``Z`` in Ω/m, and ``Y`` in S/m. The potential
coefficients of the coaxial computation correspond to ``\mathrm{j}\omega P``. The
starting points of the lines reproduce the voltage references of `:martins2026`. The reference
is the interface for a conductor in the air and infinite depth for a conductor in the earth.
Each line is shifted horizontally by ``\min(10^{-5}, 0.01 d_i)`` m from the anchor of
the terminal, where ``d_i`` is the smallest metal dimension of the terminal, and the metal
interiors add nothing to the integral. The matrices are then reduced by bundling, Kron
elimination and ideal transposition, as in the coaxial computation.

## Domain and perfectly matched layers

For each medium ``m``, the base root ``q_{0,m} = \sqrt{\mathrm{j}\omega\mu_m\widehat\sigma_m}
= a_{0,m} + \mathrm{j} b_{0,m}`` has a nonnegative real part, and the transverse root
``q_m = a_m + \mathrm{j} b_m`` of
``q_m^2 = \mathrm{j}\omega\mu_m\widehat\sigma_m - \Gamma^2`` has ``b_m \geq 0`` when
``\operatorname{Im}\Gamma < b_{0,m}`` and ``a_m \geq 0`` otherwise. The layout radius
``R_l`` is the largest of 5 m, of the distance from ``(x_c, 0)`` to each cable plus its
radius, and of twice the distance between each pair of cables, where ``x_c`` is the mean
horizontal position of the cables. The sizing length of the earth ``L`` and the half-width
``D`` of the physical domain ``|x - x_c| \leq D``, ``|y| \leq D`` are

```math
L_0 = \min\left(\frac{1}{a_{0,e}}, \frac{2\pi}{|q_{0,e}|}\right), \qquad
L = \min\left(L_0, \sqrt{\frac{2\rho_\mathrm{cap}}{\omega\mu_e}}\right), \qquad
D = \max\left(R_l, m_D L \frac{|q_{0,e}|}{\max(|q_{0,e}|, |q_e|)}\right),
```

where the index ``e`` denotes the earth, ``1/a_{0,e}`` is omitted for ``a_{0,e} = 0``,
``\rho_\mathrm{cap} = 10^5`` Ω·m, and ``m_D`` is `domain_size_factor`. The side, top and
bottom layers have the thicknesses ``L_d = f_d D``, with ``f_d`` the overrides
`PmlSideThicknessFactor`, `PmlTopThicknessFactor` and `PmlBottomThicknessFactor`.

Each layer stretches its normal coordinate by the staged profile

```math
s_d(u) = 1 + A_d\left(u^3 - \mathrm{j}\frac{9}{4}\eta_d u^8\right),
```

where ``u`` is the normalized depth in the layer, ``(|x - x_c| - D)/L_d`` on the sides and
``(|y| - D)/L_d`` at the top and bottom. The real part starts cubically, and the absorption
of propagating waves is delayed to the eighth degree. The stretches enter the materials as
``\nu\mathbf T^{-1}`` and ``\widehat\sigma\mathbf T`` in the transverse terms and as
``\widehat\sigma s_x s_y`` in the axial terms, with
``\mathbf T = \operatorname{diag}(s_y/s_x, s_x/s_y, s_x s_y)``, and the side stretch is the
same in both media, which keeps the interface matched. The layers of the sides use both
media, the top layer the air and the bottom layer the earth, denoted by the sets ``M_d``.
The imaginary stretch and the strength of each layer are

```math
\begin{aligned}
\eta_d &= \min\left(\max_{m \in M_d} \nu_m, \min_{m \in M_d} \beta_m\right), \qquad
\nu_m = \operatorname{clamp}\left(\frac{b_m - a_m}{b_m}, 0, 1\right), \qquad
\beta_m = \min\left(1, \frac{2 a_m}{9 |b_m|}\right), \\
A_d &= \frac{4 T}{L_d \min_{m \in M_d} \widehat d_{m,d}}, \qquad
\widehat d_{m,d} = \max\left(\frac{a_m + \eta_d b_m}{C_m}, 0.1 b_{0,m}\right), \qquad
C_m = \frac{a_{0,m} + b_{0,m}}{b_{0,m}}, \qquad
T = -\frac{\ln R}{2},
\end{aligned}
```

where ``\nu_m = 0`` for ``b_m \leq 0``, ``\beta_m = 1`` for ``b_m \geq 0``, and ``R`` is
`pml_reflection`. When the air field is quasi-static over 100 half-widths,
``100 |q_{0,a}| D \leq 0.1``, and the rate of a layer is below ``T / (100 D)``, the layer
takes ``\eta_d = 0`` and the rate ``\max(T / (100 D), a_m / C_m, 0.1 b_{0,m})``, which limits
the added real stretch to 100 half-widths. With ``X_d = L_d(1 + A_d/4 -
\mathrm{j}\eta_d A_d/4)``, the net exponent
``\operatorname{Re}[q_m(D + X_d)]`` must be positive in every medium of each layer, except
at the exact cutoff ``q_m = 0``. Otherwise, GetDP stops with an error.

## Mesh

The number of intervals across layer ``d`` and their progression ``e^{g_d/N_d}`` are

```math
N_d = \max\left(N_{\min,d},
    \left\lceil \frac{\max_{m \in M_d} \mathrm{PPW}_m \Phi_{m,d}}{2\pi} \right\rceil\right),
\qquad
\Phi_{m,d} = |q_m||X_d| \min\left(1, \frac{T}{E_{m,d}}\right), \qquad
E_{m,d} = \operatorname{Re}(q_m X_d),
```

where ``\Phi_{m,d} = |q_m||X_d|`` for ``E_{m,d} \leq 0``,
``\mathrm{PPW}_m = 10 \operatorname{clamp}\left(\sqrt{0.1|q_m|/a_m}, 1, 3\right)/m_h`` for
``a_m > 0`` and ``30/m_h`` otherwise, ``m_h`` is `mesh_size_factor`, ``N_{\min,d}`` are the
overrides `PmlSideLayers`, `PmlTopLayers` and `PmlBottomLayers`, and ``g_d`` are
`PmlSideGrading`, `PmlTopGrading` and `PmlBottomGrading`. The air takes part only when its
field is excited, which happens when a cable lies in the air or when the exterior of a buried
cable lies within ``6/a_e`` of the interface, and always for ``a_e = 0``.

The physical regions are meshed with triangles and the layers with quadrangles. With the
resolution radius ``R_r = \max(R_l, L)``, the bulk target ``h_b``, the remote target
``h_r``, the wave target ``h_{w,m}`` and the decay radius ``r_{d,m}`` of each medium are

```math
h_b = \frac{m_h R_r}{20}, \qquad
h_r = \min\left(f_\mathrm{ext} h_b, h_b + 0.2 m_h \max(0, D - 2 R_r)\right), \qquad
h_{w,m} = \min\left(h_b, \frac{m_h}{8 |q_m|}\right), \qquad
r_{d,m} = \min\left(2 R_r, \frac{6}{a_m}\right),
```

where ``f_\mathrm{ext}`` is `ExteriorMeshSizeFactor` and ``6/a_m`` is omitted for
``a_m \leq 0``. In a lossy earth, the exterior size grows with the distance ``d`` from the
cables as ``h(d) = \min(h_{r,e}, h_{w,e} e^{a_e d / 2})``. A mesh boundary layer on the earth side
of the interface has the first size ``h_{w,e}``, the growth ratio 1.4 and the thickness
``\min(r_{d,e}, \frac{1}{2}\min_p(|y_p| - r_p - h_p))`` over the buried cables ``p``, with
``h_p`` the target of the neighbourhood of the cable. The layer is used when
``h_{w,e} < h_{r,e}`` and its thickness is at least ``2 h_{w,e}``.

A circular arc of radius ``R`` takes at least ``\max(n_c, 2\pi R/h)`` segments per turn,
with ``h`` the target of its region and
``n_c = 12 \cdot 2^{\max(0, \lceil \log_2(\sqrt{2\pi^2/(3\tau)}/12) \rceil)}``, which bounds
the relative error of the area of a circle by ``\tau``, the override
`ConductorGeometryTolerance`. Each region has the target ``h`` equal to ``m_h`` times
``r/5`` for a solid circle of radius ``r``, ``r`` for a round strand, ``(r_o - r_i)/2`` for an
annulus, ``\min(w, h_\square)/5`` for a rectangle of width ``w`` and height ``h_\square``,
and ``\min(w, h_\square)/2`` for a rectangular strand. A thin passive annulus is meshed as a
strip of ``\max(2, \lceil 2/m_h \rceil)`` radial intervals, with the target of the region at
least ``2\pi r_o/n_c``.

In a circular or annular conductor of width ``w``, the radius of a circle or the thickness of
an annulus, the skin depth ``\delta`` and the phase rate ``\beta`` follow from
``\sqrt{\mathrm{j}\omega\mu\widehat\sigma} = 1/\delta + \mathrm{j}\beta``. A graded layer of
elements along the contour of the conductor, of depth ``e`` and ``N`` intervals growing by
the ratio ``g``, starts at the first size ``h_1``:

```math
e = \min(n_s\delta, f w), \qquad
N = \left\lceil \frac{\ln\left(1 + (g - 1) e / h_1^\ast\right)}{\ln g} \right\rceil, \qquad
h_1 = \frac{e (g - 1)}{g^N - 1}, \qquad
h_1^\ast = \min\left(\frac{m_h\delta}{n_\delta}, \frac{m_h w}{n_t}\right),
```

where ``n_s``, ``n_\delta``, ``g`` and, for an annulus, ``n_t`` are the overrides
`ConductorSkinDepths`, `ConductorSkinDepthElements`, `ConductorMeshGrowth` and
`ConductorThicknessElements`, ``n_t = 5`` for a circle, and ``f`` is 0.5 for a circle and
0.45 for an annulus. The bent rectangular strands keep the target of their region. The interior of the conductor has the target
``\min(h, m_h w / n_t)``, further limited to ``2\pi m_h / (12\beta)`` when
``w / \delta \leq n_s``.

## Computation options

The keywords of `compute(problem, FEM(); kw...)` are the computation options of
[`GmshBackend`](@ref):

| Option | Default | Meaning |
| --- | --- | --- |
| `mesh_size_factor` | `1.0` | factor ``m_h`` of every native resolution target (larger values coarsen the mesh) |
| `domain_size_factor` | `2.0` | multiplier ``m_D`` of the earth sizing length for the half-width of the domain |
| `pml_reflection` | `1e-3` | nominal attenuation target ``R`` of the perfectly matched layer, in the open interval (0, 1) |
| `workers` | `clamp(Sys.CPU_THREADS ÷ 4, 1, 8)` | number of Gmsh and GetDP processes at once |
| `solver_threads` | `1` | number of BLAS and OpenMP threads of each GetDP process |
| `keep_run_directory` | `false` | whether to keep the run directory |
| `trace` | `false` | whether to keep the matrices `Z` and `P` by terminal in the result |
| `getdp_executable` | `:discover` | the GetDP binary, a [`GetdpBinary`](@ref Main.GmshExt.GetdpBinary) choice or a path |
| `gmsh_verbosity` | `2` | Gmsh's `-v` level, from 0 to 5 |
| `getdp_verbosity` | `2` | GetDP's `-v` level, from 0 to 5 |
| `overrides` | each at its default | native expert controls by native name |

The factors must be finite and positive, and `workers` and `solver_threads` must be at
least 1. Each native control must lie inside its range, otherwise `validate` throws an
`ArgumentError`. The native controls accepted by `overrides` are:

```@docs
Main.GmshExt.OVERRIDES
```

```@eval
using LineCableModels, Gmsh, Markdown
GmshExt = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
rows = ["| `$(name)` | `$(repr(entry.default))` | `$(entry.range)` |"
        for (name, entry) in pairs(GmshExt.OVERRIDES)]
Markdown.parse(join(["| Control | Default | Range |", "| :-- | :-- | :-- |", rows...], "\n"))
```

## Run directory

Each computation writes its inputs and outputs under the current directory, in
`.linecablemodels/fem/<system id>/<yyyymmdd-HHMMSS>-<4 random characters>/`. The run
directory contains the native sources and the geometry under `input/`, the meshes under
`mesh/`, the matrices `results/Z.tsv` and `results/P.tsv`, and the logs under `logs/`, where
`logs/getdp-info.log` records the executable used in the computation. With the default
`keep_run_directory = false`, the run directory is removed after the result is built. With
`keep_run_directory = true`, the directory is kept and its path is logged. When a native
process fails, the run directory is kept, and the error message gives the frequency and the
log file with its last twenty lines.

## Example

The reference case `two_bare_wires` consists of two bare copper wires, buried at 1 m depth
with a horizontal spacing of 1 m, in a homogeneous earth of 0.1 Ω·m. The example below
computes this system at three frequencies of its range.

Build the system and compute its line parameters:

```@example fem
using LineCableModels, Gmsh
copper = Material(:conductor, 1.7241e-8, 1, 0.999994, 20, 0.00393) # [Ω·m], [°C], [1/°C]
wire = build(CableDesign, "bare_wire",
    Terminal(:core, Region(:core_metal, Circle(0.0425), copper))) # [m]
system = build(CableSystem, [wire, wire], [Placement(0.0, -1.0), Placement(1.0, -1.0)];
    connections = [(core = 1,), (core = 2,)], id = "two_bare_wires", line_length = 1.0)
frequencies = (10.0 .^ range(-1, 7; length = 101))[[1, 51, 101]] # [Hz]
problem = LineParametersProblem(system, homogeneous(rho = 0.1); frequencies) # [Ω·m]
result = compute(problem, FEM())
result.Z[1, 1, :] # [Ω/m]
```

Self admittance of the first wire at the three frequencies:

```@example fem
result.Y[1, 1, :] # [S/m]
```

## Accuracy

The test suite compares the finite-element results with the analytical ones for the six
reference cases, at each frequency, by means of the relative errors

```math
e_Z = \frac{\max_{ij} |\Delta Z_{ij}|}{\max_i |Z_{ii}|}, \qquad
e_Y = \frac{\max_{ij} |\Delta Y_{ij}|}{\max_i |Y_{ii}|}.
```

The accuracy level of the finite-element computation is 1 %. The maximum errors over the
eleven finite-element frequencies of each case are:

| Case | ``e_Z`` | ``e_Y`` |
| --- | --- | --- |
| `two_bare_wires` | 0.568 % | 0.798 % |
| `two_insulated_wires` | 0.407 % | 0.0447 % |
| `coaxial_circular_strands` | 4.92 % | 10.5 % |
| `coaxial_rectangular_strands` | 5.22 % | 10.4 % |
| `overhead_bare_wires` | 0.610 % | 0.817 % |
| `mixed_bare_wires` | 0.609 % | 0.617 % |

The differences in the coaxial cases result from the different representations of the cable
in the two models:

- ``e_Y``: the analytical model represents the armor wires as a tube. In the finite-element
  model, the wires touch the bedding of 3 mm only at their crowns, with a filler of relative
  permittivity 1 in the wedges between them.
- ``e_Z`` below about 25 Hz: the lay correction of the armor is included in the equivalent
  conductor of the analytical model and absent from the two-dimensional model.
- ``e_Z`` above about 40 kHz: the lay correction of the core, for the same reason.
- ``e_Z`` from 25 Hz to 40 kHz: the magnetic armor wires compared with the equivalent tube.

For the four wire cases, the finite-element results of this package agree with those of the
reference computation within 1e-6 of the largest diagonal entry.

## Continuous integration

The workflow `.github/workflows/FEM.yml` runs the complete test suite, including the
finite-element tests of the reference cases, on every pull request that is not a draft.

## Uncertain inputs

For any number type of the problem, the backend computes the nominal problem, in which each
number is replaced by the `Float64` value of its nominal value, and returns a
`LineParameters{Float64}` without uncertainty. The propagation of uncertainties is performed
by the coaxial computation, as described in [Uncertainty quantification](uncertainty.md).
