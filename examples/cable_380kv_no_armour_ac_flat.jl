# # 380 kV cables in flat formation
#=
This example demonstrates how to model a three-phase underground system composed of 380 kV
single-core cables without armor. Each cable consists of a stranded copper core, inner and
outer semicon screens, XLPE insulation, a water-blocking layer, a lead sheath, an inner
polyethylene layer and an outer jacket. The cables are buried at 1 m depth in a flat
formation, with a horizontal spacing of 1 m, in a homogeneous earth of 1000 Ω·m.

This example covers:

1. Representing a compacted stranded core with rectangular strands.
2. Examining the equivalent tubular conductor of the stranded core and the solenoid
   correction of the insulation permeability.
3. Computing the line parameters of the six conductors of the system.
4. Eliminating the sheaths of a system with sheaths grounded at both ends.
5. Applying the ideal transposition to the phase conductors.
=#

# Load the package:
using LineCableModels

#=
## Materials

The semicon screens are defined as materials of the `:semicon` class, with resistivities of
1000 Ω·m and 500 Ω·m for the inner and outer screens, respectively, and a relative
permittivity of 1000. The XLPE and the outer jacket include the relative permeability given
by the cable data.
=#

# Conductor, semicon and insulation materials:
core_metal = Material(:conductor, 2.3853e-8, 1.0, 1.01663, 20.0, 0.00393) # [Ω·m]
semicon_inner = Material(:semicon, 1_000.0, 1_000.0) # [Ω·m]
semicon_outer = Material(:semicon, 500.0, 1_000.0) # [Ω·m]
xlpe = Material(:insulator, 1.0e14, 2.77781, 1.59891) # [Ω·m]
lead = Material(:conductor, 2.14e-7, 1.0, 0.999983, 20.0, 0.004) # [Ω·m], [°C], [1/°C]
pe = Material(:insulator, 1.97e14, 2.51862) # [Ω·m]
pp = Material(:insulator, 1.0e15, 2.8, 1.1263) # [Ω·m]

#=
## Stranded core with rectangular strands

The core consists of 127 copper wires with a diameter of 3.66 mm, arranged as a central
wire and six layers of 6, 12, 18, 24, 30 and 36 wires, all with a lay ratio of 15. In the
actual cable, the wires are compacted to fill a circle with the total area of the 127 wires.
This compacted core is represented with rectangular strands: layer `k` occupies a ring
between two radii, and each of its `6k` rectangular strands has the area of one wire. The
[`CircularPattern`](@ref) bends each rectangle into an annular sector of the same area, so
that the central wire and the sectors tile the full cross-section of the conductor.
=#

# Strand radius, strand area and outer radius of each layer:
strand_radius = 3.66e-3 / 2 # strand_diameter / 2 [m]
strand_area = π * strand_radius^2 # [m²]
rings = [strand_radius * sqrt(1 + 3k * (k + 1)) for k in 0:6] # [m]

# Layers of rectangular strands, laid on helices around the central wire:
courses = [Lay(Helix(LayRatio(15.0)), # core_lay_ratio
               Repeat(CircularPattern(6k; r = Offset(0)),
                   Region(Symbol(:core_strands_, k),
                       Rectangle(strand_area / (rings[k + 1] - rings[k]),
                           rings[k + 1] - rings[k]),
                       core_metal)))
           for k in 1:6]
rings[end] # [m]

#=
## Cable design

The central wire is declared with the [`core`](@ref) builder, followed by the strand
layers in the same [`@terminal`](@ref) block. The [`screen`](@ref) builder declares the
semicon screens and the water-blocking layer, the [`sheath`](@ref) builder declares the lead
sheath, and the [`jacket`](@ref) builder declares the inner polyethylene and the outer
jacket.
=#

# Build the cable design and inspect its outer radius:
design = @cable "cable_380kv_no_armour" begin
    @terminal :core begin
        core(core_metal; r = strand_radius, tag = :core_center)
        courses...
    end
    screen(semicon_inner; t = 2.0e-3, tag = :core_semicon_inner) # inner_semicon_thickness [m]
    insulation(xlpe; t = 20.0e-3, tag = :core_insulation) # insulation_thickness [m]
    screen(semicon_outer; t = 1.5e-3, tag = :core_semicon_outer) # outer_semicon_thickness [m]
    screen(semicon_outer; t = 0.3e-3, tag = :core_water_blocking) # water_blocking_thickness [m]
    @terminal :sheath begin
        sheath(lead; t = 3.3e-3, tag = :sheath_lead) # lead_sheath_thickness [m]
    end
    jacket(pe; t = 3.0e-3, tag = :sheath_inner_pe) # inner_pe_thickness [m]
    jacket(pp; t = 10.0e-3, tag = :sheath_outer_jacket) # jacket_thickness [m]
end
r_ex(design.outer) # [m]

#=
## Lay parameters of the core

Each strand layer corresponds to one row in the table of lays of the design, containing the
lay radius of the rectangular strands, the pitch of the helix and the overlength, defined as the length
of a strand per unit length of cable. For a constant lay ratio, the overlength is the same
for all layers.
=#

# Lay radius, pitch and overlength of each layer:
[(radius = row.radii[1], pitch = row.pitches[1], overlength = row.overlengths[1])
 for row in design.lays] # [m], [m], [dimensionless]

#=
## Equivalent coaxial representation

The function `LineCableModels.flatten` reduces each terminal to one equivalent tubular
conductor. The central wire and the six layers conduct in parallel. The central wire is
straight, whereas the strands of each layer are longer than the cable by their overlength,
which increases their resistance per unit length of cable. The equivalent tube of the core
then has a resistivity higher than that of copper, and its relative permeability
reproduces the geometric mean radius of the central wire and the rectangular strands.
=#

# Equivalent conductors, with radii in m and resistivity in Ω·m:
blueprint = LineCableModels.flatten(design)
[(terminal = conductor.terminal, r_in = conductor.r_in, r_ex = conductor.r_ex,
     rho = conductor.material.rho, mu_r = conductor.material.mu_r)
 for conductor in blueprint.conductors] # [m], [m], [Ω·m], [dimensionless]

# Ratio of the core equivalent resistivity to the copper resistivity:
blueprint.conductors[1].material.rho / core_metal.rho

#=
The dielectric layers between the core and the sheath form the insulation interval of the
core. The semicon screens are included in this interval as dielectrics with a relative
permittivity of 1000. Due to this high permittivity, their contribution to the potential
coefficient of the interval is small compared to that of the XLPE.
=#

# Dielectric layers of the core and sheath intervals:
[(conductor = layer.conductor, r_in = layer.r_in, r_ex = layer.r_ex,
     eps_r = layer.material.eps_r, mu_r = layer.material.mu_r)
 for layer in blueprint.dielectrics] # [m], [m], [dimensionless], [dimensionless]

#=
The current in a stranded core follows the helices of its rectangular strands and also
circulates around the cable axis, like the current in a solenoid. This effect introduces
additional magnetic energy in the insulation interval. The relative permeability of each
interval, `mu_r_ins`, combines this solenoid correction with the mean permeability of the
layers. Since the sheath has no lay, the value of its interval is the mean permeability of
its layers.
=#

# Relative permeability of the core and sheath intervals:
blueprint.mu_r_ins

#=
## Three-phase system in flat formation

The cables are placed at `x = −1, 0, 1` m. Each terminal is assigned its own conductor
index, with core and sheath of cable `i` connected to conductors `2i − 1` and `2i`,
respectively.
=#

# Homogeneous earth model and frequency range:
earth = homogeneous(rho = 1_000.0, eps_r = 10.0) # [Ω·m]
frequencies = 10.0 .^ range(-1, 7; length = 101) # [Hz]

# Cable positions and cable system for a line 20 km long:
placements = [Placement(x, -1.0) for x in (-1.0, 0.0, 1.0)] # [m]
system = build(CableSystem, fill(design, 3), placements;
    connections = [(core = 2i - 1, sheath = 2i) for i in 1:3],
    id = "cable_380kv_no_armour_ac_flat", line_length = 20_000.0) # [m]

# Line parameters problem at 70 °C:
problem = LineParametersProblem(system, earth; temperature = 70.0, frequencies) # [°C]

#=
## Line parameters at power frequency

The bundle and Kron reductions are enabled by default. Since no two terminals share an
index and no terminal has the index zero, the result retains the six conductors of the
system. The matrices at the frequency closest to 50 Hz are presented below, with `Z` in Ω/km
and `Y` in µS/km, rounded.
=#

# Compute the line parameters and inspect the row labels:
result = compute(problem, CoaxialFormulation())
result.terminals

# Frequency closest to 50 Hz:
k = argmin(abs.(result.frequencies .- 50.0))
result.frequencies[k] # [Hz]

# Series resistance matrix:
round.(real(result.Z[:, :, k]) .* 1e3; digits = 4) # [Ω/km]

# Series reactance matrix:
round.(imag(result.Z[:, :, k]) .* 1e3; digits = 4) # [Ω/km]

# Shunt conductance matrix:
round.(real(result.Y[:, :, k]) .* 1e9; digits = 4) # [µS/km]

# Shunt susceptance matrix:
round.(imag(result.Y[:, :, k]) .* 1e9; digits = 1) # [µS/km]

#=
## Sheaths grounded at both ends

Sheaths grounded at both ends of the line are represented by setting their voltage to zero.
This condition is declared with the index 0, and the Kron reduction eliminates the sheaths
from the matrices. The resulting system retains the three cores, with the effect of the
sheath currents included in the reduced matrices.
=#

# Cable system with grounded sheaths:
grounded = build(CableSystem, fill(design, 3), placements;
    connections = [(core = i, sheath = 0) for i in 1:3],
    id = "cable_380kv_grounded_sheaths", line_length = 20_000.0) # [m]
grounded_problem = LineParametersProblem(grounded, earth; temperature = 70.0, frequencies)

# Compute the line parameters and inspect the row labels:
grounded_result = compute(grounded_problem, CoaxialFormulation())
grounded_result.terminals

# Series impedance of the three cores:
grounded_result.Z[:, :, k] # [Ω/m]

#=
## Ideal transposition

Since the central cable of a flat formation occupies a position different from the outer
cables, the matrices are not symmetric under a cyclic permutation of the phases. A
transposed line permutes the positions of the phases along its length. The option
`transposition = true` replaces each matrix by the mean over the cyclic permutations of its
rows and columns, so that the entries of each cyclic diagonal become equal.
=#

# Compute the transposed line parameters:
transposed = compute(grounded_problem, CoaxialFormulation(transposition = true))
transposed.Z[:, :, k] # [Ω/m]

# Size of the impedance array and number of frequencies:
size(result.Z), length(result.frequencies)

#=
## Conclusion

This example has demonstrated how to:

1. Represent a compacted stranded core with rectangular strands bent into annular sectors.
2. Examine the lay parameters, the equivalent resistivity and the solenoid correction of
   the stranded core.
3. Compute the line parameters of a three-phase system in flat formation.
4. Eliminate grounded sheaths with the Kron reduction.
5. Apply the ideal transposition to the reduced matrices.
=#
