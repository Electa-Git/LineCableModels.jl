# Finite elements with Gmsh and GetDP

```@meta
CurrentModule = LineCableModels
```

The Gmsh extension computes the series impedances and shunt admittances of a problem with the
finite-element method. Gmsh builds and meshes the cross-section, and GetDP solves each
frequency. The extension then reduces the matrices of the terminals as the coaxial
computation does.

## Loading

`using Gmsh` loads the extension. The computation runs GetDP 3.5.0 in complex arithmetic,
which [`GetdpBinary`](@ref Main.GmshExt.GetdpBinary) selects:

- `:discover`, the default, runs the GetDP found on `PATH` when its `-info` reports version
  3.5 and complex arithmetic, and the downloaded one otherwise.
- `:download` runs getdp.info's GetDP 3.5.0 complex build, downloaded by its pinned hash.
- A path given as a string runs as is.

The download happens once per Julia depot, at the first use, and needs the network then. The
package's artifact covers Linux on x86-64, macOS on x86-64, macOS on Apple silicon through
Rosetta 2, and Windows on x86-64. The finite-element chain on Windows and on macOS through
Rosetta 2 is untested.

## The formulation

[`FEM`](@ref) selects the material laws of a finite-element computation. Its keywords:

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

The trait [`backend`](@ref) of [`FEM`](@ref) is [`GmshBackend`](@ref), so
`compute(problem, FEM())` runs the extension.

## Computation options

The keywords of `compute(problem, FEM(); kw...)` are the computation options of
[`GmshBackend`](@ref):

| Option | Default | Meaning |
| --- | --- | --- |
| `mesh_size_factor` | `1.0` | factor of every native resolution target (larger values coarsen the mesh) |
| `domain_size_factor` | `2.0` | multiplier of the earth length scale for the half-width of the domain |
| `pml_reflection` | `1e-3` | nominal attenuation target of the perfectly matched layer, in the open interval (0, 1) |
| `workers` | `clamp(Sys.CPU_THREADS ÷ 4, 1, 8)` | number of Gmsh and GetDP processes at once |
| `solver_threads` | `1` | number of BLAS and OpenMP threads of each GetDP process |
| `keep_run_directory` | `false` | whether to keep the run directory |
| `trace` | `false` | whether to keep the matrices `Z` and `P` by terminal in the result |
| `getdp_executable` | `:discover` | the GetDP binary, a [`GetdpBinary`](@ref Main.GmshExt.GetdpBinary) choice or a path |
| `gmsh_verbosity` | `2` | Gmsh's `-v` level, from 0 to 5 |
| `getdp_verbosity` | `2` | GetDP's `-v` level, from 0 to 5 |
| `overrides` | each at its default | native expert controls by native name |

The factors must be finite and positive, `workers` and `solver_threads` at least 1, and a
native control inside its range, or `validate` throws an `ArgumentError`. The native
controls of `overrides`, with their defaults and ranges:

```@docs
Main.GmshExt.OVERRIDES
```

```@example fem
using LineCableModels, Gmsh
using Markdown
GmshExt = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
Markdown.MD(Markdown.Table(
    [[["Control"], ["Default"], ["Range"]];
     [[Markdown.Code(string(name)), Markdown.Code(repr(entry.default)),
          Markdown.Code(string(entry.range))]
      for (name, entry) in pairs(GmshExt.OVERRIDES)]],
    [:l, :l, :l]))
```

## The run directory

Each computation writes its inputs and outputs under the current directory, in
`.linecablemodels/fem/<system id>/<yyyymmdd-HHMMSS>-<4 random characters>/`. Among them are
the native sources and the geometry under `input/`, the meshes under `mesh/`, the matrices
`results/Z.tsv` and `results/P.tsv`, and the logs under `logs/`, with `logs/getdp-info.log`
recording the executable that ran. With `keep_run_directory = false`, the default, the run
directory is removed after the result is built. With `true` it stays, and its path is logged.
A failed native process keeps the run directory, and its error message gives the frequency and
the log with its last twenty lines.

## An example

The reference case `two_bare_wires` has two bare copper wires, buried 1 m deep and 1 m apart in
earth of 0.1 Ω·m. The example computes it at three frequencies of its scan:

```@example fem
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

```@example fem
result.Y[1, 1, :] # [S/m]
```

## Accuracy

The test tree compares the finite-element results with the analytical ones of the six
reference cases, per frequency, as

```math
e_Z = \frac{\max_{ij} |\Delta Z_{ij}|}{\max_i |Z_{ii}|}, \qquad
e_Y = \frac{\max_{ij} |\Delta Y_{ij}|}{\max_i |Y_{ii}|}.
```

The accuracy level of the finite-element computation is 1 %. The maxima over the eleven
finite-element frequencies of each case:

| Case | ``e_Z`` | ``e_Y`` |
| --- | --- | --- |
| `two_bare_wires` | 0.568 % | 0.798 % |
| `two_insulated_wires` | 0.407 % | 0.0447 % |
| `coaxial_circular_strands` | 4.92 % | 10.5 % |
| `coaxial_rectangular_strands` | 5.22 % | 10.4 % |
| `overhead_bare_wires` | 0.610 % | 0.817 % |
| `mixed_bare_wires` | 0.609 % | 0.617 % |

The coaxial cases differ by the differences between the two models:

- ``e_Y``: the analytical model treats the armor wires as a tube. The wires meet the 3 mm
  bedding only at their crowns, with fill of relative permittivity 1 in the wedges between
  them.
- ``e_Z`` below about 25 Hz: the lay correction of the armor, which the analytical model
  includes in its equivalent conductor and a two-dimensional model leaves out.
- ``e_Z`` above about 40 kHz: the lay correction of the core, in the same way.
- ``e_Z`` from 25 Hz to 40 kHz: the magnetic armor wires against the tube.

For the four wire cases, the finite-element rows of this package agree with those of the
source computation within 1e-6 of the largest diagonal entry.

## Continuous integration

The workflow `.github/workflows/FEM.yml` runs the whole test suite with the finite-element
tests of the reference cases on every pull request that is not a draft.

## Uncertain inputs

For any number type of the problem, the backend computes the nominal problem, each number
replaced by the `Float64` value of its nominal value, and returns a `LineParameters{Float64}`
without uncertainty. The coaxial computation propagates the uncertainty, as
[Uncertainty](uncertainty.md) describes.
