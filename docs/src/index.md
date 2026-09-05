# LineCableModels.jl

[`LineCableModels.jl`](https://github.com/Electa-Git/LineCableModels.jl)
calculates frequency-domain electrical parameters for underground and overhead
cable systems. The models include conductor skin effect, dielectric loss,
earth return, frequency-dependent earth properties, and declared uncertainty
in geometry and material data.

## Documentation

- [Tutorials](tutorials.md) introduce cable construction and calculation.
- [Modelling and results](usage.md) covers calculations, result access,
  uncertainty, tables, and plots.
- [Gridspace and uncertainty](gridspace.md) specifies finite variation and
  uncertainty realisation.
- [API reference](reference.md) lists the line and cable calculation API.
- [Conveniences](conveniences.md) covers estimates, scalar formulas, and VDE
  designation parsing.
- [Benchmarks](gauntlet.md) publishes the Gauntlet validation results and plots.
- [Developers](developers.md) records grammar invariants, CI checks, extension
  APIs, and project conventions.

## Features

- Construct deterministic and uncertain designs with the typed
  `Grid`/`Gridspace` grammar and evaluate them with `compute`.
- Calculate base cable parameters for solid, tubular, and stranded cores,
  semiconductors, screens, armors, sheaths, tapes, and water-blocking materials.
- Apply temperature, stranding, and twisting corrections to DC resistance
  [app14198982](@cite), GMR [6521501](@cite), and base inductance
  [yang2008gmr](@cite).
- Calculate dielectric loss and equivalent insulation resistance
  [916943](@cite), including the solenoid contribution of twisted strands to
  insulation permeability [5743045](@cite).
- Assemble phase-domain Z/Y matrices for polyphase systems with any number of
  conductors per phase, with optional Measurements-based direct propagation or
  conditional Monte Carlo analysis.
- Replace stranded assemblies with equivalent tubular conductors for EMT
  calculations and export cable data to ATPDraw and PSCAD.
- Calculate internal impedance for solid, tubular, and multilayer coaxial
  single-core cables with the formulation in [4113884](@cite) or the
  approximations documented by
  [PSCAD](https://www.pscad.com/webhelp/EMTDC/Transmission_Lines/Deriving_System_Y_and_Z_Matrices.htm).
- Calculate earth-return impedance and admittance for underground conductors in
  homogeneous soil from the electric-Hertz-vector solution of the Helmholtz
  equation, up to 10 MHz [5437464](@cite).

## Installation

Install the registered package from Julia's package manager:

```julia-repl
pkg> add LineCableModels
```

Then load the core package:

```julia
using LineCableModels
```

Plotting is optional. Load one backend explicitly before calling `preview` or
`plot`:

```julia
using LineCableModels
using CairoMakie
```

## User statistics

![Top Julia package-server regions observed for LineCableModels.jl](assets/user-statistics.svg)

The map is generated in CI from Julia's public package-server request logs. It
shows the top server regions by the sum of `request_addrs` for requests marked
as user traffic. These regional aggregates are not a count of distinct people
and do not represent country-level telemetry.


## License

The source code is licensed under the
[BSD 3-Clause License](https://github.com/Electa-Git/LineCableModels.jl/LICENSE).

---
```@raw html
<p align="left">Documentation generated using <a target="_blank" href="https://github.com/JuliaDocs/Documenter.jl">Documenter.jl</a> and <a target="_blank" href="https://github.com/fredrikekre/Literate.jl">Literate.jl</a>.</p>
```
