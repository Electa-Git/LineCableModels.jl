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
- Earth-return calculations, modal transformations, and ATPDraw/PSCAD data
  exchange.
- Material and cable libraries with JSON import and export.
- One typed `Grid`/`Gridspace` grammar for deterministic and uncertain designs.
- Ordinary, combinatorial, and conditional Monte Carlo execution through
  `compute`.
- Optional Makie plotting through Julia package extensions.

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

The physical model is built directly from materials, primitives, and cable
parts:

```julia
copper = Material(; kind=:conductor, rho=1.7241e-8)
xlpe = Material(; kind=:insulator, rho=1.97e14, eps_r=2.5)

root = Stack(
    Group(:core, Conductor.Solid(:core_metal, copper; r=10e-3)),
    Insulator.Shell(:insulation, xlpe; t=8e-3),
)
design = build(CableDesign, "example", root)

constants = CableConstants(design)
```

Scalar `build` calls return completed domain objects. An explicit `Grid` at a
construction boundary returns a `Gridspace` whose points invoke the same
scalar `build` method.

## Optional plotting

Load one Makie backend explicitly before calling `preview` or `plot`:

```julia
using LineCableModels
using CairoMakie

plot_pages = plot(line_parameters) # R, X, G, B matrix-dashboard pages
export_svg(plot_pages[1]; path = "series_resistance.svg")

self_impedance = plot(line_parameters, @observe Z[1, 1, :])
# `self_impedance` is one UIPlot with separate R and X axes.
```

`GLMakie` and `WGLMakie` are supported in the same way. LineCableModels
never imports a backend dynamically. The optional `backend=:cairo`, `:gl`, or
`:wgl` keyword is explicit call-local sugar. Observation requests determine
physical quantity/coordinate axes; `layout` only groups those facets into
figures. Singular recipes return a `UIPlot`; calls that produce several figures
return a vector. Each handle exposes native Makie objects that remain
caller-owned. `figure_title` names the whole figure, `panel_titles` names its
logical axes, and `series_labels` names overlaid result containers in legends.
Matrix coordinates belong in semantic axis titles, never legend entries.
Legends accept outer docks or
`legend_position=:inside` with an anchor such as `:rt`; `figurelegend!`,
`panellegend!`, `figuretitle!`, and `paneltitle!` can change those native blocks
after construction. The Export SVG control saves their current live state and
requires CairoMakie to have been loaded explicitly.

## Result access

`CableConstants` stores R/L/C values per metre. `LineParameters`
stores its frequency domain and either a `:pul` or `:total` basis:

```julia
basis(line_parameters)
@observe line_parameters R[1, 1, :]       # complete frequency response
@observe line_parameters Z[1, 1, 2:5]     # selected frequency samples
@observe line_parameters (Z, abs)[1, 2, :]

label(R)                         # "Series resistance"
symbol(Z, angle)                 # "∠Z"
label(display_unit(R, :pul))     # "Ω/km"
```

Complete parameter traversals return `ParametricResult{T}`, including a space
with cardinality one. `Formulation`, `CableConstantsFormulation`, and
`ModalTransformationFormulation` accept explicit `Grid` fields; combinatorial
calculation forms the Cartesian product of problem and formulation points and
retains both axes for `result(run, problem_index, formulation_index)` lookup.
Conditional Monte Carlo propagation returns `MonteCarloResult{T}`. Loading
PolyChaos.jl enables validated non-intrusive propagation through
`PolynomialChaos` and `PolynomialChaosResult{T}`. Use `statistics`,
`expansions`, and `validation` for its fitted products; PCE results do not
fabricate empirical samples or histograms. Result order is
problem-index-fastest within formulation order. Unresolved traversal state is
not copied into completed results. Parametric, linear-error, and Monte Carlo
results are ordinary finite collections: indexing and iteration return stored
core results, and Base
`first`, `last`, `only`, `collect`, `map`, and `zip` retain their standard
meanings.
`DataFrame(monte_carlo_result)` renders marginal summaries. After loading a
Makie package, `Makie.hist`, `Makie.stairs`, `Makie.ecdfplot`, `Makie.lines`,
and `Makie.qqplot` display retained distribution information through the
LineCableModels shell.

Higher-order calculations keep supplemental computation output separate from
their scientific products. `details(result)` returns the empty named tuple by
default. Construct `Combinatorial`, `LinearError`, or `MonteCarlo` with
`options=(retain_details=true,)` only when the core computation owner has
registered a `computation_details` method and those records are needed.
Monte Carlo uses strict failure propagation by default. Physically unsupported
draws can be rejected explicitly with
`options=(retain_details=true, on_error=:retry, max_failures=100)`. Only
`DomainError` is retryable; retained details report every rejected argument
tuple and its error summary. The resulting statistics are conditional on a
successful realisation.

Physics and numerical-method choices belong to `Formulation`. Execution choices
are passed as a named tuple. For a materialised line system, pass
`options=(output_basis=:total,)` to `compute` to scale both Z and Y by the line
length. Composite calculations select their operation explicitly, for example
`compute(ParametricProblem(space), Combinatorial(Formulation()))`.

Human-facing XLSX output is a ReportBuilder operation:

```julia
using XLSX

artifact = report(
    XLSXReportDefinition(file_name="line_parameters.xlsx"),
    line_parameters,
)
artifact.output
```

XLSX is an optional dependency; loading it activates the workbook writer. The
established `export_data(:xlsx, line_parameters; ...)` call delegates to the
same report and returns its output path. Relative paths and the default
`ZY_export.xlsx` resolve from the caller's current working directory.

## Retired FEM and sector support

FEM/GetDP support and sector-shaped cable support were removed before the
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
