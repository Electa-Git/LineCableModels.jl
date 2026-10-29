# LineCableModels.jl

<img src="docs/src/assets/logo.svg" align="left" width="150" alt="LineCableModels.jl logo" />

[![Coverage](https://codecov.io/github/Electa-Git/LineCableModels.jl/branch/main/graph/badge.svg)](https://app.codecov.io/github/Electa-Git/LineCableModels.jl)
[![Aqua QA](https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://electa-git.github.io/LineCableModels.jl/stable/)
[![CI](https://github.com/Electa-Git/LineCableModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Electa-Git/LineCableModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](LICENSE)

`LineCableModels.jl` is a Julia package for computing the series impedances and shunt admittances per unit length of arbitrary
arrangements of underground and overhead power cables, with uncertainty propagation. A computation is defined by a problem, which describes the cables, the earth, the operating temperature and the frequencies of interest, and by a formulation, which selects the formulas for each physical contribution. The computation is carried out either by the coaxial analytical backend of the package or by the finite-element backend of the Gmsh extension.

## Main features

- **High-fidelity cable modeling:** Detailed representation of solid and stranded conductors,
  with round wires or rectangular strands laid on helices, semicon screens, insulation layers,
  sheaths, beddings, armors of round wires, jackets and fillers. Any number of cables can be
  placed in the same cross-section to form a cable system.
- **Equivalent tubular representation:** Each terminal of a cable is represented by an
  equivalent tubular conductor that preserves the DC resistance and the geometric mean radius
  of its regions, including the overlength of the round wires and rectangular strands laid on
  helices. The relative permeability of the insulation layers is corrected to account for the
  solenoid effect introduced by the helical path of the current in stranded conductors.
- **Internal impedances and admittances:** Internal impedances of multilayered coaxial cables
  based on the surface impedances of Schelkunoff and on the general formulation of Ametani,
  with conductor resistivities corrected to the operating temperature.
- **Earth-return impedances and admittances:** Earth-return impedances and potential
  coefficients of underground, overhead and mixed systems in homogeneous earth, computed by a
  unified formulation with a prescribed longitudinal wavenumber.
- **Phase-domain matrices:** Series impedance and shunt admittance matrices of polyphase
  systems with any number of conductors per phase, with bundling of bonded conductors, Kron
  elimination of grounded conductors and ideal transposition.
- **Finite element modeling:** Computation of the same problems on a mesh of the actual
  cross-section of the cables, by means of Gmsh and GetDP.
- **Uncertainty quantification:** Linear propagation of the uncertainties of materials,
  dimensions, positions, earth properties and temperature to the line parameters, by means of
  the [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl) package.

## Installation

Add the registered package to the Julia environment:

```julia-repl
pkg> add LineCableModels
```

Then, in your Julia code, import the package:

```julia
using LineCableModels
```

The extensions are loaded together with their packages, namely `using Measurements` for the
propagation of uncertainties and `using Gmsh` for the finite-element backend. The
finite-element backend uses the GetDP executable found on `PATH` when it is suitable, and
otherwise downloads GetDP 3.5.0 at its first use.

## Example

The following example computes the line parameters of two insulated copper wires, buried at
1 m depth with a horizontal spacing of 1 m, in a homogeneous earth of 100 Ω·m:

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

## Finite element solver

After `using Gmsh`, the same problem is computed with Gmsh and GetDP by
`compute(problem, FEM())`. The
[finite-element guide](https://electa-git.github.io/LineCableModels.jl/dev/fem/) describes the
options of the backend, its run directory and its accuracy on the reference systems.

## Uncertainty quantification

After `using Measurements`, any input of the problem can be given as an uncertain value, such
as `Circle(4.0e-3 ± 1.0e-5)`, and the resulting `Z` and `Y` have uncertain real and imaginary
parts. The [uncertainty guide](https://electa-git.github.io/LineCableModels.jl/dev/uncertainty/)
describes the propagation and the tests that verify it.

## Documentation

See the [full documentation](https://electa-git.github.io/LineCableModels.jl) for detailed
usage instructions, technical background and examples.

## Contributing

Contributions are welcome. [CONTRIBUTING.md](CONTRIBUTING.md) gives the principles, the
checks and the form of a pull request, and the developer section of the documentation
describes the architecture, the rules and the conventions.

## License

The source code is provided under the [BSD 3-Clause License](LICENSE).

## Acknowledgements

This work is supported by the Etch Competence Hub of EnergyVille, financed by the Flemish Government. The primary developer is Amauri Martins ([@amaurigmartins](https://github.com/amaurigmartins)).

<p align = "left">
  <p><br><img src="docs/src/assets/img/ETCH_LOGO_RGB_COLOR.svg" width="150" alt="Etch logo"></p>
  <p><img src="docs/src/assets/img/ENERGYVILLE-LOGO.svg" width="150" alt="EV logo"></p>
  <p><img src="docs/src/assets/img/kul_logo.svg" width="150" alt="KUL logo"></p>
</p>
