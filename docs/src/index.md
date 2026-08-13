# LineCableModels.jl

[`LineCableModels.jl`](https://github.com/Electa-Git/LineCableModels.jl) is a specialized Julia package designed to compute the electrical parameters of coaxial arbitrarily-layered underground/overhead cables with uncertainty quantification. It focuses on calculating line and cable impedances and admittances in the frequency-domain, accounting for skin effect, insulation properties, and earth-return impedances with frequency-dependent soil models.

## Documentation outline

```@contents
Pages = [
    "index.md",
    "tutorials.md",
    "reference.md",
    "bib.md",
]
Depth = 1
```

## Features

- Calculates all base DC parameters of a given cable design (R, L, C and G), for solid, tubular or stranded cores, semiconductors, screens, armors, sheaths, tapes,  and water-blocking materials, with uncertainty propagation using the [Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl) package.
- Correction factors to account for temperature, stranding and twisting effects on the DC resistance [app14198982](@cite), GMR [6521501](@cite) and base inductance of stranded cores and wire screens [yang2008gmr](@cite).
- Explicit computation of dielectric losses and effective resistances for insulators and semiconductors [916943](@cite). Correction of the magnetic constant of insulation layers to account for the solenoid effect introduced by twisted strands [5743045](@cite).
- Computes phase-domain Z/Y matrices for poliphase systems with any number of conductors per phase, and sequence-domain components for three-phase systems, with uncertainty propagation.
- Improved equivalent tubular representation for EMT simulations and direct export to ATPDraw and PSCAD formats.
- Computes internal impedances of solid, tubular or coaxial multi-layered single-core (SC) cables, using rigorous [4113884](@cite) or equivalent approximate formulas available in [industry-standard EMT software](https://www.pscad.com/webhelp/EMTDC/Transmission_Lines/Deriving_System_Y_and_Z_Matrices.htm).
- Computes earth-return impedances and admittances of underground conductors in homogeneous soil, based on a rigorous solution of Helmholtz equation on the electric Hertzian vector, valid up to 10 MHz [5437464](@cite).

## Installation

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

## License

The source code is provided under the [BSD 3-Clause License](https://github.com/Electa-Git/LineCableModels.jl/LICENSE).

---
```@raw html
<p align="left">Documentation generated using <a target="_blank" href="https://github.com/JuliaDocs/Documenter.jl">Documenter.jl</a> and <a target="_blank" href="https://github.com/fredrikekre/Literate.jl">Literate.jl</a>.</p>
```