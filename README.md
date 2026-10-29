# LineCableModels.jl

<img src="docs/src/assets/logo.svg" align="left" width="150" alt="LineCableModels.jl logo" />

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://electa-git.github.io/LineCableModels.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://electa-git.github.io/LineCableModels.jl/dev/)
[![Build Status](https://github.com/Electa-Git/LineCableModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Electa-Git/LineCableModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![codecov](https://codecov.io/gh/Electa-Git/LineCableModels.jl/graph/badge.svg?token=6H12DDBZ0T)](https://codecov.io/gh/Electa-Git/LineCableModels.jl)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](LICENSE)

LineCableModels.jl computes the series impedances and shunt admittances of cables per unit
length. The computation takes a problem, which describes the cables, the earth, the operating
temperature and the frequencies, and a formulation, which selects the formulas. The backend
that executes it is the coaxial analytical computation of the package or the finite-element
computation of the Gmsh extension. The Measurements extension propagates the uncertainty of
the inputs.

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

## Example

A pair of insulated copper wires, 1 m apart and 1 m deep in earth of 100 Ω·m:

```julia
using LineCableModels

copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
polyethylene = Material(:insulator, 1.0e14, 2.3)
wire = @cable "insulated_wire" begin
    @terminal :core begin
        core(copper; r = 4.0e-3) # [m]
    end
    insulation(polyethylene; t = 1.0e-3) # [m]
end
system = build(CableSystem, [wire, wire], [Placement(0.0, -1.0), Placement(1.0, -1.0)];
    connections = [(core = 1,), (core = 2,)], id = "two_wires", line_length = 1000.0)
problem = LineParametersProblem(system, homogeneous(rho = 100.0);
    frequencies = [50.0, 1e3, 1e6]) # [Hz]
result = compute(problem, CoaxialFormulation())
result.Z[:, :, 1] # [Ω/m] at 50 Hz
result.terminals
```

## Finite elements

After `using Gmsh`, `compute(problem, FEM())` computes the same problem with Gmsh and GetDP.
The [finite-element guide](https://electa-git.github.io/LineCableModels.jl/dev/fem/) gives
its options, its run directory and its accuracy on the reference systems.

## Uncertainty

After `using Measurements`, an input such as `Circle(4.0e-3 ± 1.0e-5)` gives `Z` and `Y` with
uncertain real and imaginary parts. The
[uncertainty guide](https://electa-git.github.io/LineCableModels.jl/dev/uncertainty/) gives
what propagates and how the release checks it.

## Documentation

The documentation is at <https://electa-git.github.io/LineCableModels.jl>.

## License

The source code is licensed under the [BSD 3-Clause License](LICENSE).
