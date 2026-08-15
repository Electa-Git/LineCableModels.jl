# LineCableModels.jl

<img src="docs/src/assets/logo.svg" align="left" width="150" alt="LineCableModels.jl logo" />

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://electa-git.github.io/LineCableModels.jl/dev/)
[![Build Status](https://github.com/Electa-Git/LineCableModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Electa-Git/LineCableModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](LICENSE)
[![codecov](https://codecov.io/gh/Electa-Git/LineCableModels.jl/graph/badge.svg?token=6H12DDBZ0T)](https://codecov.io/gh/Electa-Git/LineCableModels.jl)

LineCableModels.jl computes electrical parameters for underground and overhead
power cables and supports uncertainty propagation through cable geometry and
material data.

## Features

- Solid, tubular, stranded, screened, and armored cable models.
- Frequency-dependent series impedance and shunt admittance calculations.
- Earth-return, modal transformation, ATPDraw, and PSCAD integration.
- Material and cable libraries with JSON import and export.
- Optional Makie plotting through explicit backend extensions.

## Installation

After version 0.2.0 is accepted into the Julia General registry:

```julia-repl
pkg> add LineCableModels
```

Until registration is complete, install the release commit from GitHub:

```julia-repl
pkg> add https://github.com/Electa-Git/LineCableModels.jl
```

Core usage has no plotting dependency:

```julia
using LineCableModels
```

## Optional plotting

Load one Makie backend explicitly before calling `preview`, `plot`, or
`set_backend!`:

```julia
using LineCableModels
using CairoMakie

set_backend!(:cairo)

plots = plot(line_parameters) # Vector{UIPlot}, one page per quantity
export_svg(first(plots); path = "series_resistance.svg")
```

`GLMakie` and `WGLMakie` are supported in the same way. LineCableModels
never imports or selects a backend dynamically. `preview` and Monte Carlo
distribution plots return one `UIPlot`; line-parameter plots always return a
`Vector{UIPlot}`. The Export SVG control preserves the current declarative plot
state and requires CairoMakie to have been loaded explicitly.

## Result access

`CableConstants` stores canonical per-metre R/L/C values. `LineParameters`
stores its frequency domain and either a `:per_length` or `:total` basis:

```julia
basis(line_parameters)
R(line_parameters, 1, 1)       # complete frequency response
Z(line_parameters, 1, 1, 2:5) # selected frequency samples
abs.(Z(line_parameters, 1, 1))
```

Monte Carlo results use first-class `CableConstantsMC` and `LineParametersMC`
containers. Use `statistics`, `mean`, `std`, `quantile`, `samples`, `trial`,
`distribution`, and `surrogate`; joint `trial`/`rand` calls require retained
samples.

## Retired FEM and sector support

FEM/GetDP integration and sector-shaped cable support were removed before the
v0.2 release. Their final snapshot is branch `legacy/fem-sector` at commit
`b75dd2723f90a83ec090b20605ea42af57f4a9c3`. To use that historical version in
the current Julia project:

```sh
julia --project=. -e "using Pkg; Pkg.add(Pkg.PackageSpec(url=ARGS[1], rev=ARGS[2]))" https://github.com/Electa-Git/LineCableModels.jl.git b75dd2723f90a83ec090b20605ea42af57f4a9c3
```

See the [documentation](https://electa-git.github.io/LineCableModels.jl/) and
[examples](examples) for supported workflows.

## License and citation

LineCableModels.jl is distributed under the [BSD 3-Clause License](LICENSE).
Citation metadata is provided in [CITATION.cff](CITATION.cff).

## Acknowledgements

This work is supported by the Etch Competence Hub of EnergyVille, financed by
the Flemish Government. The primary developer is Amauri Martins
([@amaurigmartins](https://github.com/amaurigmartins)).

<p align="left">
  <br><img src="assets/img/ETCH_LOGO_RGB_NEG.svg" width="150" alt="Etch logo">
  <br><img src="assets/img/ENERGYVILLE-LOGO.svg" width="150" alt="EnergyVille logo">
  <br><img src="assets/img/kul_logo.svg" width="150" alt="KU Leuven logo">
</p>
