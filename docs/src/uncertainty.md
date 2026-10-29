# Uncertainty quantification

```@meta
CurrentModule = LineCableModels
```

The Measurements extension propagates the uncertainties of the inputs to the line parameters
by linear propagation, the first-order method of
[Measurements.jl](https://github.com/JuliaPhysics/Measurements.jl). The extension is loaded
by `using Measurements`.

## Measurements extension

The extension defines only methods, of the functions of the package and of SpecialFunctions.
The Bessel recurrences used in these methods are those of the NIST Digital Library of
Mathematical Functions [NIST:DLMF](@cite).

```@docs
Main.MeasurementsExt
```

## Uncertain inputs

Uncertain values, of type `Measurement`, are accepted for any number of the data model and of
the problem, such as a material property, a dimension of a cable design, a placement in the
system, a property of the earth or the operating temperature. The number type `T` of a design is a `Measurement` when any of
its numbers is uncertain, and the problem promotes its numbers to the same type. In this
case, the computation returns `Z` and `Y` with uncertain real and imaginary parts.

The correlations between results that depend on the same uncertain input are retained by
Measurements. Wires that share the same uncertain radius give correlated line
parameters, and the test suite verifies that one uncertain radius shared by two wires
gives √2 times the uncertainty obtained with two independent radii of the same value.

## Deterministic controls

- Γ and the other options are controls of the formulation or of the backend, and take
  deterministic values. A `Measurement` value given to an option throws a `MethodError` when
  the record of options is built.
- The finite-element backend computes the nominal problem and returns a result without
  uncertainty, as described in [Finite element solver](fem.md).

## Example

The following example computes the reference case `two_insulated_wires` with uncertain values
of the wire radius and of the copper resistivity.

Define the materials and the cable design with the uncertain inputs:

```@example uncertainty
using LineCableModels, Measurements
copper = Material(:conductor, 1.7241e-8 ± 1.0e-10, 1, 0.999994, 20, 0.00393) # [Ω·m], [°C], [1/°C]
insulation = Material(:insulator, 1.97e14, 2.3) # [Ω·m]
metal = Terminal(:core, Region(:core_metal, Circle(0.0425 ± 0.0005), copper)) # [m]
design = build(CableDesign, "insulated_wire", metal,
    Region(:core_insulation, Offset(1.0e-3), insulation)) # [m]
nothing # hide
```

Build the system of two wires and compute its line parameters:

```@example uncertainty
system = build(CableSystem, [design, design], [Placement(0.0, -1.0), Placement(1.0, -1.0)];
    connections = [(core = 1,), (core = 2,)], id = "two_insulated_wires", line_length = 1.0)
frequencies = 10.0 .^ range(-1, 7; length = 101) # [Hz]
problem = LineParametersProblem(system, homogeneous(rho = 0.1); frequencies) # [Ω·m]
result = compute(problem, CoaxialFormulation())
k = argmin(abs.(result.frequencies .- 50.0))
Measurements.value(result.frequencies[k]) # [Hz]
```

Real and imaginary parts of the self impedance of the first wire at this frequency:

```@example uncertainty
(R = real(result.Z[1, 1, k]), X = imag(result.Z[1, 1, k])) # [Ω/m]
```

Real and imaginary parts of the self admittance of the first wire at this frequency:

```@example uncertainty
(G = real(result.Y[1, 1, k]), B = imag(result.Y[1, 1, k])) # [S/m]
```

Uncertainty of the self resistance:

```@example uncertainty
Measurements.uncertainty(real(result.Z[1, 1, k])) # [Ω/m]
```

## Tests of the propagation

The tests in `test/measurements/propagation.jl` compute the reference cases with uncertain
inputs and verify that:

- The nominal values of `Z` and `Y` equal the deterministic results within 1e-14 of the
  largest diagonal entry.
- The derivatives with respect to each input equal the central finite differences of the
  deterministic computation with a step of 1e-5, within 1e-5.
- The uncertainties of independent inputs add as the root sum of squares, and one uncertain
  input shared by two wires gives √2 times the uncertainty obtained with two independent
  inputs of the same value, both within 1e-12.

## Computational cost

Each result evaluated inside the quadrature of the earth-return integrals stores its list of
derivatives, with one entry per uncertain input, and the computation allocates memory for
these lists at each evaluation of an integrand. On one workstation, for the reference case
`two_insulated_wires` at 11 frequencies, the deterministic computation takes 6 ms, against
83 ms with one uncertain material property and 236 ms with six uncertain inputs. The example
above computes 101 frequencies in about 1 s, against 0.054 s without uncertainty.
