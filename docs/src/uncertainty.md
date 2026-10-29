# Uncertainty

```@meta
CurrentModule = LineCableModels
```

The Measurements extension propagates the uncertainty of the inputs to the line parameters by
linear propagation, the first-order method of
[Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl). `using Measurements`
loads it.

## The extension

The extension defines methods only, on the names of the package and of SpecialFunctions. The
Bessel recurrences it uses are those of the NIST Digital Library of Mathematical Functions
[NIST:DLMF](@cite).

```@docs
Main.MeasurementsExt
```

## What propagates

Any number of the data model and of the problem can be a `Measurement`. A material property
or a dimension of a design can be uncertain, and so can a placement in the system or a
property of the earth. The temperature can be uncertain too. The number type `T` of a design
is a `Measurement` when any of its numbers is, and the problem promotes its numbers to the
same type. The computation then gives `Z` and `Y` with uncertain real and imaginary parts.
Measurements keeps the correlations of a repeated source: two wires that share one uncertain
radius give correlated parameters, and the test tree checks that one source repeated in two
wires gives √2 times the uncertainty of two independent sources of the same value.

## What does not

- Γ and every other option are controls of the formulation or of the backend, not measured
  inputs. A `Measurement` value for an option throws a `MethodError` when the options record
  is built.
- The finite-element backend computes the nominal problem and returns a result without
  uncertainty, as [Finite elements](fem.md) describes.
- This release has no Monte Carlo propagation.

## An example

The insulated wire pair of the reference case `two_insulated_wires`, with an uncertain wire
radius and an uncertain resistivity of the copper:

```@example uncertainty
using LineCableModels, Measurements
copper = Material(:conductor, 1.7241e-8 ± 1.0e-10, 1, 0.999994, 20, 0.00393) # [Ω·m], [°C], [1/°C]
insulation = Material(:insulator, 1.97e14, 2.3) # [Ω·m]
metal = Terminal(:core, Region(:core_metal, Circle(0.0425 ± 0.0005), copper)) # [m]
design = build(CableDesign, "insulated_wire", metal,
    Region(:core_insulation, Offset(1.0e-3), insulation)) # [m]
system = build(CableSystem, [design, design], [Placement(0.0, -1.0), Placement(1.0, -1.0)];
    connections = [(core = 1,), (core = 2,)], id = "two_insulated_wires", line_length = 1.0)
frequencies = 10.0 .^ range(-1, 7; length = 101) # [Hz]
problem = LineParametersProblem(system, homogeneous(rho = 0.1); frequencies) # [Ω·m]
result = compute(problem, CoaxialFormulation())
k = argmin(abs.(result.frequencies .- 50.0))
result.frequencies[k] # [Hz]
```

```@example uncertainty
result.Z[1, 1, k] # [Ω/m]
```

```@example uncertainty
result.Y[1, 1, k] # [S/m]
```

```@example uncertainty
Measurements.uncertainty(real(result.Z[1, 1, k])) # [Ω/m]
```

## Tests of the propagation

The tests of `test/measurements/propagation.jl` compute the reference cases with uncertain
inputs and check three facts:

- The nominal parts of `Z` and `Y` equal the deterministic computation within 1e-14 of the
  largest diagonal entry.
- The derivatives with respect to each input equal central finite differences of the
  deterministic computation at the step 1e-5, within 1e-5.
- The uncertainty of independent sources adds as the root sum of squares, and one source
  repeated in two wires gives √2 times the uncertainty of two independent sources of the same
  value, both within 1e-12.

## Cost

Each result inside the quadrature of the earth-return integrals stores its list of
derivatives, one entry per uncertain source. The computation allocates memory for these lists
at each evaluation of an integrand. Measured on one workstation, on the reference case
`two_insulated_wires` at 11 frequencies, the deterministic computation takes 6 ms, with one
uncertain material property 83 ms, and with six uncertain sources 236 ms. The example above
computes 101 frequencies in about 1 s against 0.054 s without uncertainty.
