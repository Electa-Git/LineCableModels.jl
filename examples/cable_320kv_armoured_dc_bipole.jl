# # 320 kV armored cables of a bipole
#=
This example demonstrates how to model a bipole composed of two 320 kV armored single-core
cables. Each cable has the stranded copper core and the insulation system of the 380 kV
design, followed by a lead sheath, a bedding and an armor of 68 steel wires under a jacket.
The cables are buried at 1 m depth, with a horizontal spacing of 1 m, in a homogeneous earth
of 1000 Ω·m.

This example covers:

1. Declaring an armor of round wires over a bedding.
2. Examining the equivalent tubular conductor of the armor.
3. Computing the line parameters of the six conductors of the bipole.
4. Bonding the sheath and the armor of each cable into one conductor.
5. Separating the internal impedance of the cables from the earth-return contribution by
   means of the trace of the computation.
=#

# Load the package:
using LineCableModels

#=
## Materials

The materials of the core, insulation system, sheath and jackets are the same as in the
380 kV example. The armor wires are made of galvanized steel with a relative permeability
of 36.6.
=#

# Conductor, semicon, insulation and armor materials:
core_metal = Material(:conductor, 2.3853e-8, 1.0, 1.01663, 20.0, 0.00393) # [Ω·m]
semicon_inner = Material(:semicon, 1_000.0, 1_000.0) # [Ω·m]
semicon_outer = Material(:semicon, 500.0, 1_000.0) # [Ω·m]
xlpe = Material(:insulator, 1.0e14, 2.77781, 1.59891) # [Ω·m]
lead = Material(:conductor, 2.14e-7, 1.0, 0.999983, 20.0, 0.004) # [Ω·m], [°C], [1/°C]
pe = Material(:insulator, 1.97e14, 2.51862) # [Ω·m]
pp = Material(:insulator, 1.0e15, 2.8, 1.1263) # [Ω·m]
steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045) # [Ω·m], [°C], [1/°C]

#=
## Stranded core

The core is identical to that of the 380 kV example, composed of a central wire and six
layers of rectangular strands. Each strand is bent into an annular sector with the area of
one wire of 3.66 mm.
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
## Cable design with armor

The [`armor`](@ref) builder lays `n` round wires in one layer over the outline of the part
beneath it, in this case the bedding, with the lay ratio `lay`. The space between the wires
is completed with an insulating filler. The thickness of the bedding is chosen so that the
68 wires fit side by side in their layer.
=#

# Build the cable design and inspect its outer radius:
design = @cable "cable_320kv_armoured" begin
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
    bedding(pp; t = 0.009376520171274612, tag = :sheath_bedding) # bedding_thickness [m]
    @terminal :armor begin
        armor(steel; shape = Circle(5.82e-3 / 2), n = 68, # armor_wire_diameter / 2 [m]
            lay = LayRatio(15.0), tag = :armor_wires) # armor_lay_ratio
    end
    jacket(pp; t = 10.0e-3, tag = :armor_jacket) # jacket_thickness [m]
end
r_ex(design.outer) # [m]

#=
The last row in the table of lays corresponds to the armor, containing the radius of the
wire centers, the pitch of the helix and the overlength of the wires.
=#

# Lay radius, pitch and overlength of the armor wires:
armor_lay = last(design.lays)
(radius = armor_lay.radii[1], pitch = armor_lay.pitches[1],
    overlength = armor_lay.overlengths[1]) # [m], [m], [dimensionless]

#=
## Equivalent coaxial representation

The function `LineCableModels.flatten` reduces the armor to a tubular conductor that spans
the layer of wires. The resistance of the tube equals that of the 68 wires along their
helices. Its relative permeability, which reproduces the geometric mean radius of the wires,
differs from the permeability of steel.
=#

# Equivalent conductors, with radii in m and resistivity in Ω·m:
blueprint = LineCableModels.flatten(design)
[(terminal = conductor.terminal, r_in = conductor.r_in, r_ex = conductor.r_ex,
     rho = conductor.material.rho, mu_r = conductor.material.mu_r)
 for conductor in blueprint.conductors] # [m], [m], [Ω·m], [dimensionless]

#=
## Bipole system

The cables are placed at `x = ±0.5` m. The core, sheath and armor of cable `i` are
connected to conductors `3i − 2`, `3i − 1` and `3i`, respectively.
=#

# Homogeneous earth model and frequency range:
earth = homogeneous(rho = 1_000.0, eps_r = 10.0) # [Ω·m]
frequencies = 10.0 .^ range(-1, 7; length = 101) # [Hz]

# Cable positions and cable system for a line 150 km long:
placements = [Placement(-0.5, -1.0), Placement(0.5, -1.0)] # [m]
system = build(CableSystem, [design, design], placements;
    connections = [(core = 3i - 2, sheath = 3i - 1, armor = 3i) for i in 1:2],
    id = "cable_320kv_armoured_dc_bipole", line_length = 150_000.0) # [m]

# Line parameters problem at 70 °C:
problem = LineParametersProblem(system, earth; temperature = 70.0, frequencies) # [°C]

# Compute the line parameters and inspect the row labels:
result = compute(problem, CoaxialFormulation())
result.terminals

#=
## Line parameters at power frequency

Since no two terminals share an index and no terminal has the index zero, the bundle and
Kron reductions leave the six conductors unchanged. The matrices at the frequency closest to
50 Hz are presented below, with `Z` in Ω/km and `Y` in µS/km, rounded.
=#

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
## Bonded sheath and armor

Terminals that share a conductor index are bonded: they have the same voltage and their
currents add. With the default option `bundle = true`, the computation merges the
bonded terminals into one conductor, labeled with the names of both terminals.
=#

# Cable system with the sheath and armor of each cable bonded:
bundled = build(CableSystem, [design, design], placements;
    connections = [(core = 2i - 1, sheath = 2i, armor = 2i) for i in 1:2],
    id = "cable_320kv_bundled_screens", line_length = 150_000.0) # [m]
bundled_problem = LineParametersProblem(bundled, earth; temperature = 70.0, frequencies)

# Compute the line parameters and inspect the row labels:
bundled_result = compute(bundled_problem, CoaxialFormulation())
bundled_result.terminals

# Series resistance matrix of the bonded system:
round.(real(bundled_result.Z[:, :, k]) .* 1e3; digits = 4) # [Ω/km]

#=
## Internal and earth-return matrices

The computation of the [`LCMBackend`](@ref LineCableModels.LCMBackend) accepts the option
`trace`. With `trace = true`, the result retains the four matrices from which `Z` and `Y`
are assembled at each frequency:

- `Zin`, the series impedance in Ω/m inside each cable, from its conductors to its outer
  surface, per conductor.
- `Pin`, the matrix of potential coefficients in m/F inside each cable, per conductor.
- `Zg`, the earth-return impedance in Ω/m from the outer surfaces of the cables, per cable.
- `Pg`, the matrix of earth-return potential coefficients in m/F, per cable.
=#

# Compute the line parameters with the trace and inspect the size of each matrix:
traced = compute(problem, CoaxialFormulation(); trace = true)
map(size, traced.trace)

# Earth-return impedance of the two cables:
traced.trace.Zg[:, :, k] # [Ω/m]

#=
The entry of `Zg` for a pair of cables is added to every pair of their conductors.
In the absence of a reduction, `Z` is then the sum of `Zin` and this expansion of `Zg`.
The largest absolute difference between the two, in Ω/m, is computed below.
=#

# Expand `Zg` to the conductors and compare the sum with `Z`:
expanded = [traced.trace.Zg[cld(i, 3), cld(j, 3), k] for i in 1:6, j in 1:6] # [Ω/m]
maximum(abs, traced.trace.Zin[:, :, k] + expanded - traced.Z[:, :, k]) # [Ω/m]

# Size of the impedance array and number of frequencies:
size(result.Z), length(result.frequencies)

#=
## Conclusion

This example has demonstrated how to:

1. Declare an armor of round wires over a bedding.
2. Examine the lay parameters and the equivalent tubular conductor of the armor.
3. Compute the line parameters of a bipole with core, sheath and armor in each cable.
4. Bond the sheath and the armor of each cable into one conductor.
5. Retrieve the internal and earth-return matrices of the computation and verify that they
   add up to the series impedance.
=#
