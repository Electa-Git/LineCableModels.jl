# # Single-core cable with a 1000 mm² solid aluminum conductor
#=
This example demonstrates how to build the model of a single-core cable and compute its
series impedance and shunt admittance per unit length. The objective is to represent a cable
with a solid aluminum conductor of 1000 mm² and an XLPE insulation 8.3 mm thick, buried at
1 m depth in a homogeneous earth of 100 Ω·m, over the frequency range from 0.1 Hz to 10 MHz.

This example covers:

1. Defining the materials of the cable.
2. Building the cable design with the [`@cable`](@ref) and [`@terminal`](@ref) blocks.
3. Examining the equivalent coaxial representation of the design.
4. Defining the cable system and the line parameters problem.
5. Computing the line parameters and extracting R, L and C at the power frequency.
6. Examining the frequency dependence of the series resistance.
=#

# Load the package:
using LineCableModels

#=
## Materials

A [`Material`](@ref) is defined by its class (`:conductor`, `:insulator` or `:semicon`),
followed by the electrical resistivity, relative permittivity, relative permeability,
reference temperature and temperature coefficient of the resistivity. Omitted arguments take
their default values.
=#

# Aluminum conductor and XLPE insulation:
aluminum = Material(:conductor, 2.8264e-8, 1.0, 1.000022, 20.0, 0.00429) # [Ω·m], [°C], [1/°C]
xlpe = Material(:insulator, 1.97e14, 2.5) # [Ω·m]

#=
## Cable design

The [`@cable`](@ref) block lists the parts of the cable from the inside out. The
[`@terminal`](@ref) block assigns the conductive regions declared inside it to one terminal,
in this case the core. The [`core`](@ref) builder declares a solid circular conductor by its
radius, and the [`insulation`](@ref) builder declares a layer by its thickness, placed over
the outline of the part beneath it.
=#

# Radius of a solid core with a cross-section of 1000 mm²:
core_radius = sqrt(1000.0e-6 / π) # [m]

# Build the cable design and inspect its outer radius:
design = @cable "solid_1000mm2" begin
    @terminal :core begin
        core(aluminum; r = core_radius, tag = :core_metal)
    end
    insulation(xlpe; t = 8.3e-3, tag = :core_insulation) # insulation_thickness [m]
end
r_ex(design.outer) # [m]

#=
## Equivalent coaxial representation

Coaxial formulations represent each terminal of the cable by an equivalent tubular
conductor, with the dielectric layers located between consecutive conductors. The function
`LineCableModels.flatten` returns this representation of the design, which in this case
consists of one solid conductor and one insulation layer.
=#

# Equivalent conductor, with radii in m and resistivity in Ω·m:
blueprint = LineCableModels.flatten(design)
[(terminal = conductor.terminal, r_in = conductor.r_in, r_ex = conductor.r_ex,
     rho = conductor.material.rho) for conductor in blueprint.conductors]

# Insulation layer of the core:
[(conductor = layer.conductor, r_in = layer.r_in, r_ex = layer.r_ex,
     eps_r = layer.material.eps_r) for layer in blueprint.dielectrics]

#=
## Cable system and problem definition

A [`CableSystem`](@ref) places one or more cable designs in the cross-section, with the
vertical coordinate `y` negative below the earth surface. The connections assign a conductor
index to each terminal, here index 1 to the core. The [`LineParametersProblem`](@ref)
completes the definition with the earth model, the conductor temperature and the frequencies
of interest.
=#

# Homogeneous earth model:
earth = homogeneous(rho = 100.0, eps_r = 10.0) # [Ω·m]

# Cable buried at 1 m depth, in a line 1 km long:
system = build(CableSystem, [design], [Placement(0.0, -1.0)]; # [m]
    connections = [(core = 1,)], id = "solid_1000mm2_single", line_length = 1000.0) # [m]

# Frequency range from 0.1 Hz to 10 MHz, logarithmically spaced:
frequencies = 10.0 .^ range(-1, 7; length = 101) # [Hz]

# Line parameters problem at 20 °C:
problem = LineParametersProblem(system, earth; temperature = 20.0, frequencies) # [°C]
length(problem.frequencies)

#=
## Line parameters at power frequency

The [`CoaxialFormulation`](@ref) selects one formula for each physical contribution to the
series impedance and shunt admittance, by default the formulas implemented in this release.
The [`compute`](@ref) function returns a [`LineParameters`](@ref) object containing the
matrices `Z`, in Ω/m, and `Y`, in S/m, for each frequency, together with the label of each
row.
=#

# Compute the line parameters and inspect the row labels:
result = compute(problem, CoaxialFormulation())
result.terminals

# Frequency closest to 50 Hz, since the samples do not include 50 Hz exactly:
k = argmin(abs.(result.frequencies .- 50.0))
result.frequencies[k] # [Hz]

# Series impedance at this frequency:
result.Z[:, :, k] # [Ω/m]

# Shunt admittance at this frequency:
result.Y[:, :, k] # [S/m]

#=
The series resistance `R` and inductance `L` per unit length are obtained from the real and
imaginary parts of `Z`, respectively, and the capacitance `C` from the imaginary part of
`Y`, with ``\omega = 2 \pi f``. Since `Z` is the impedance of the loop formed by the core
and the earth, `R` and `L` include the earth-return contribution. The capacitance `C`
results from the potential coefficients of the insulation and of the surrounding earth.
=#

# Per-unit-length R, L and C at the power frequency:
ω = 2π * result.frequencies[k] # [rad/s]
(R = real(result.Z[1, 1, k]), L = imag(result.Z[1, 1, k]) / ω,
    C = imag(result.Y[1, 1, k]) / ω) # [Ω/m], [H/m], [F/m]

#=
## Frequency dependence of the series resistance

The DC resistance of the core at 20 °C is the ratio of its resistivity to its cross-section.
At 0.1 Hz, the loop resistance is close to this value. As the frequency increases, the
current in the core concentrates toward its outer surface due to the skin effect, and the
return current in the earth concentrates closer to the cable, both increasing the loop
resistance.
=#

# DC resistance of the core and loop resistance at 0.1 Hz, 50 Hz and 10 MHz:
(R_dc = aluminum.rho / (π * core_radius^2), R_first = real(result.Z[1, 1, 1]),
    R_50Hz = real(result.Z[1, 1, k]), R_last = real(result.Z[1, 1, end])) # [Ω/m]

# Size of the impedance array and number of frequencies:
size(result.Z), length(result.frequencies)

#=
## Conclusion

This example has demonstrated how to:

1. Define the materials and the parts of a single-core cable.
2. Obtain the equivalent coaxial representation of the cable design.
3. Define a cable system and a line parameters problem in homogeneous earth.
4. Compute the per-unit-length line parameters and extract R, L and C at the power
   frequency.
5. Observe the increase of the series resistance with frequency.

The example of the 320 kV bipole separates the matrices of the cables from the earth-return
contribution, which isolates the internal impedance of the conductors.
=#
