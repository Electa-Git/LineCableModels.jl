# # Positive-sequence parameters of 18/30 kV cables in trefoil
#=
This example computes the positive-sequence series impedance and shunt admittance per unit
length of a three-phase line of 18/30 kV single-core cables. The frequencies are the
harmonics of 50 Hz up to the 50th, at 2500 Hz. Each cable has an aluminum core of 630 mm²
and an XLPE insulation 8 mm thick between two semiconducting screens. A copper screen under
a polyethylene jacket covers the insulation. The cables touch each other in a trefoil
formation, with the lower two at 1 m depth. The copper screens connect to the earth at both
ends of the line, and the line has an ideal transposition of its phases.

This example covers:

1. Building the cable design with the [`@design`](@ref) and [`@terminal`](@ref) blocks.
2. Placing three cables in a trefoil formation with [`trefoil`](@ref).
3. Eliminating the grounded screens and transposing the phases in the formulation.
4. Computing the positive-sequence parameters at the harmonic frequencies.
=#

# Load the package:
using LineCableModels

#=
## Materials

The core is aluminum and the screen copper. The XLPE insulation and the polyethylene jacket
have a resistivity of 10¹⁴ Ω·m and a relative permittivity of 2.3. The semiconducting screens
have a resistivity of 1000 Ω·m and a relative permittivity of 1000.
=#

# Conductors, semiconducting screens and insulators:
aluminum = Material(:conductor, 2.8264e-8, 1.0, 1.0, 20.0, 0.00403) # [Ω·m], [°C], [1/°C]
copper = Material(:conductor, 1.7241e-8, 1.0, 1.0, 20.0, 0.00393) # [Ω·m], [°C], [1/°C]
semicon = Material(:semicon, 1000.0, 1000.0) # [Ω·m]
xlpe = Material(:insulator, 1.0e14, 2.3) # [Ω·m]
polyethylene = Material(:insulator, 1.0e14, 2.3) # [Ω·m]

#=
## Cable design

The core is a solid conductor of radius 14.16 mm, the radius of a circle of 630 mm². The
copper screen is a tube 0.3 mm thick. Each [`@terminal`](@ref) block declares one conductor of
the cable: the core, and the copper screen under the name `:sheath`.
=#

# Build the cable design and inspect its outer radius:
design = @design "xlpe_18_30kv_1x630al" begin
    @terminal :core begin
        core(aluminum; r = 14.16e-3, tag = :core) # [m]
    end
    screen(semicon; t = 0.8e-3, tag = :semicon_inner) # [m]
    insulation(xlpe; t = 8.0e-3, tag = :insulation) # [m]
    screen(semicon; t = 0.8e-3, tag = :semicon_outer) # [m]
    @terminal :sheath begin
        sheath(copper; t = 0.3e-3, tag = :screen) # [m]
    end
    jacket(polyethylene; t = 2.5e-3, tag = :jacket) # [m]
end
r_ex(design.outer) # [m]

#=
## Trefoil formation

The [`trefoil`](@ref) function places three cables on the vertices of an equilateral
triangle of side `spacing`, one at the top and two at the base, about the center of the
formation `at`. The cables touch each other when the spacing equals their outer diameter.
The center lies at a third of the height of the triangle above its base, at
``s / (2\sqrt{3})`` for the spacing ``s``. With the center at that height above 1 m depth,
the two lower cables lie at 1 m depth.
=#

# Touching cables, the lower two at 1 m depth:
spacing = 2r_ex(design.outer) # [m]
placements = trefoil(design; spacing, at = Placement(0.0, -1.0 + spacing / (2sqrt(3)))) # [m]

#=
## Cable system and problem definition

The connections assign the conductor indices 1, 2 and 3 to the cores and the index zero to
the screens. With the index zero, a conductor stays at the potential of the earth along the
whole line, as a screen grounded at both ends does. The formulation eliminates these conductors
from the matrices by Kron reduction, and the currents induced in the screens remain in the
impedance of the cores. The line parameters are per unit length, so the line length does not
change them.
=#

# Cores as conductors 1 to 3, screens grounded:
system = build(CableSystem, fill(design, 3), placements;
    connections = [(core = k, sheath = 0) for k in 1:3], id = "cable_30kv_630mm2_trefoil",
    line_length = 1.0) # [m]

# Homogeneous earth and the harmonics of 50 Hz up to the 50th:
earth = homogeneous(rho = 100.0, eps_r = 10.0) # [Ω·m]
frequencies = 50.0:50.0:2500.0 # [Hz]
problem = LineParametersProblem(system, earth; temperature = 20.0, frequencies) # [°C]
length(problem.frequencies)

#=
## Positive-sequence parameters

With `transposition = true`, the [`CoaxialFormulation`](@ref) applies the ideal cyclic
transposition of the three phases. Each matrix of the result then has one value ``z_s`` on
its diagonal and one value ``z_m`` off it, and the positive-sequence impedance is
``z_1 = z_s - z_m``. The positive-sequence admittance ``y_1 = y_s - y_m`` follows in the same
way from the admittance matrix.
=#

# Compute the line parameters and inspect the row labels:
result = compute(problem, CoaxialFormulation(transposition = true))
result.terminals

# Positive-sequence impedance and admittance at every frequency:
z1 = result.Z[1, 1, :] .- result.Z[1, 2, :] # [Ω/m]
y1 = result.Y[1, 1, :] .- result.Y[1, 2, :] # [S/m]

# Resistance, inductance and capacitance at 50 Hz:
ω = 2π * result.frequencies[1] # [rad/s]
(R1 = real(z1[1]), L1 = imag(z1[1]) / ω, C1 = imag(y1[1]) / ω) # [Ω/m], [H/m], [F/m]

#=
## Resistance at the harmonics

The positive-sequence resistance increases with the harmonic order ``h``, through the skin
effect in the core and the currents induced in the grounded screens. The output below gives
the ratio of the resistance at the harmonics 5, 7, 11, 13, 25 and 50 to the resistance at
50 Hz. Beside it stands ``\sqrt{h}``, the growth of the resistance of a conductor whose skin
depth is small against its radius.
=#

# Resistance at selected harmonics relative to 50 Hz:
[(h = h, ratio = real(z1[h]) / real(z1[1]), sqrt_h = sqrt(h))
 for h in (5, 7, 11, 13, 25, 50)]

#=
## Conclusion

This example has demonstrated how to:

1. Place three cables in a trefoil formation.
2. Remove the grounded screens with the connection index zero.
3. Apply the ideal transposition in the formulation and get the positive-sequence
   parameters.
4. Compute the line parameters at the harmonics of the power frequency.
=#
