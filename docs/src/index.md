# LineCableModels.jl

[`LineCableModels.jl`](https://github.com/Electa-Git/LineCableModels.jl) computes the series
impedances and shunt admittances of cables per unit length. The computation takes a problem,
which describes the cables, the earth, the operating temperature and the frequencies, and a
formulation, which selects the formulas. The backend that executes it is the coaxial
analytical computation of the package or the finite-element computation of the Gmsh
extension. The Measurements extension propagates the uncertainty of the inputs.

## Documentation

- [Tutorials](tutorials.md) build four cable systems and compute their line parameters.
- [Cable data model](data-model.md) describes materials, shapes, parts, designs and systems.
- [Line parameters](computing.md) covers problems, formulations, computations and results.
- [Finite elements](fem.md) covers the Gmsh extension and its options.
- [Uncertainty](uncertainty.md) covers the propagation of input uncertainty.
- [Theory](theory/contents.md) gives the matrix formulation and the formulas.
- [API reference](reference.md) lists the public names.
- [Developers](developers/architecture.md) records the architecture, the rules, the ways of
  working and the history of the package.

## Features

- Cables declared by their parts: materials, solid and stranded conductors, round wires and
  rectangular strands on helices, layers, fills, terminals and systems of several cables.
- A coaxial formulation with Schelkunoff's surface impedances, Ametani's insulation term,
  lossless admittivities, the series chain of the potential coefficients and the unified earth
  return with a prescribed longitudinal wavenumber Γ.
- Bundle and Kron reductions and the ideal transposition of the conductor matrices.
- A finite-element backend with Gmsh and GetDP.
- Linear propagation of input uncertainty with Measurements.
- The test tree checks the analytical and finite-element computations against six frozen
  reference systems.
- `Float32` inputs computed in `Float64`.

## Installation

Install the registered package from Julia's package manager:

```julia-repl
pkg> add LineCableModels
```

Then load the package:

```julia
using LineCableModels
```

The extensions load with their packages: `using Measurements` for the propagation of
uncertainty, and `using Gmsh` for the finite-element backend. The finite-element backend runs
the GetDP found on `PATH` when it is suitable, and downloads GetDP 3.5.0 on its first use
otherwise.

## Work in progress

The features of the previous releases are refactored into this structure one by one, never
merged from the old code:

- the earth-return formulas of Carson, Pollaczek, Wise, Saad, Wedepohl and Wilcox, Gary,
  Lucca, Ametani 2009 and Xue
- the lossy dielectric law, pipe-type cables and ducts, sector and Milliken conductors, tapes,
  rectangular armor wires and compaction
- layered, frequency-dependent and equivalent homogeneous earth
- cable constants, modal decomposition and finite-length propagation
- grids and parametric studies, Monte Carlo evaluation, observables and plotting through
  Makie, reports and tables
- exports to ATPDraw, PSCAD and TRALIN, material and cable libraries, saving and loading
- the resume of a run, mesh reuse, the ONELAB export, field maps and the GMRES solver in the
  finite-element backend

## User statistics

![Top Julia package-server regions observed for LineCableModels.jl](assets/user-statistics.svg)

The map is generated in CI from Julia's public package-server request logs. It
shows the top server regions by the sum of `request_addrs` for requests marked
as user traffic. These regional aggregates are not a count of distinct people
and do not represent country-level telemetry.

## License

The source code is licensed under the
[BSD 3-Clause License](https://github.com/Electa-Git/LineCableModels.jl/blob/main/LICENSE).

---
```@raw html
<p align="left">Documentation generated using <a target="_blank" href="https://github.com/JuliaDocs/Documenter.jl">Documenter.jl</a> and <a target="_blank" href="https://github.com/fredrikekre/Literate.jl">Literate.jl</a>.</p>
```
