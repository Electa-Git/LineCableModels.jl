# LineCableModels.jl

<img src="docs/src/assets/logo.svg" align="left" width="150" alt="LineCableModels.jl logo" />

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://electa-git.github.io/LineCableModels.jl/dev/)
[![Build Status](https://github.com/Electa-Git/LineCableModels.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Electa-Git/LineCableModels.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![codecov](https://codecov.io/gh/Electa-Git/LineCableModels.jl/graph/badge.svg?token=6H12DDBZ0T)](https://codecov.io/gh/Electa-Git/LineCableModels.jl)

`LineCableModels.jl` is a Julia package for computing the electrical parameters of arbitrary arrangements of underground and overhead power cables, with built-in uncertainty quantification. It is designed as a general-purpose and scalable toolbox to calculate transmission line parameters and to construct models for steady-state analysis and electromagnetic transient (EMT) simulations.
  
## Main features

- **Comprehensive cable modeling:** Detailed representation of conductors (solid, tubular, stranded), insulation layers, screens, armoring, and semicons.
- **Line and cable constants:** Accurate DC and AC parameters (R, L, C, G) with correction factors for temperature, stranding, and helical effects.
- **Propagation characteristics:** Rigorous electromagnetic   models for cable internal impedances and earth-return paths.
- **Multiple solvers:** Analytical formulations, finite element modeling, and interfaces to EMT programs, including ATPDraw and PSCAD.
- **Materials and cables library:** Store and reuse standardized material properties and cable designs across projects.

## Documentation

See the [full documentation](https://electa-git.github.io/LineCableModels.jl/) for detailed usage instructions, technical background, and examples.

## Usage

Clone the package and add to the Julia environment:

```julia-repl
pkg> add https://github.com/Electa-Git/LineCableModels.jl.git
```

If you are using the finite-element solver, it is recommended to run the build script to retrieve the binaries needed by the [GetDP.jl](https://github.com/Electa-Git/GetDP.jl) front-end:

```julia-repl
pkg> build LineCableModels
```

Then, in your Julia code, import the package:

```julia
using LineCableModels
```

For application examples, please refer to the [tutorials section](https://electa-git.github.io/LineCableModels.jl/) and the [examples folder](examples).

## License

The source code is provided under the [BSD 3-Clause License](LICENSE).

## Acknowledgements

This work is supported by the Etch Competence Hub of EnergyVille, financed by the Flemish Government. The primary developer is Amauri Martins ([@amaurigmartins](https://github.com/amaurigmartins)).

<p align = "left">
  <p><br><img src="assets/img/ETCH_LOGO_RGB_NEG.svg" width="150" alt="Etch logo"></p>
  <p><img src="assets/img/ENERGYVILLE-LOGO.svg" width="150" alt="EV logo"></p>
  <p><img src="assets/img/kul_logo.svg" width="150" alt="KUL logo"></p>
</p>
