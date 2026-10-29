# # 380 kV cables in flat formation
#
# Each of these three cables has a stranded copper core, semicon screens, XLPE insulation, a
# water-blocking layer, a lead sheath and two polymer jackets, with no armor. The cables lie
# 1 m deep and 1 m apart in a flat formation, in earth of 1000 Ω·m. The tutorial computes
# the line parameters of the six conductors, then of the three cores with the sheaths
# grounded, and the ideal transposition of that system.

using LineCableModels

# ## Materials
#
# The semicon screens are materials of the class `:semicon`, with a resistivity of 1000 Ω·m
# inside the insulation and 500 Ω·m outside it, and a relative permittivity of 1000. Each
# insulator also takes the relative permeability that the design gives.

core_metal = Material(:conductor, 2.3853e-8, 1.0, 1.01663, 20.0, 0.00393) # [Ω·m]
semicon_inner = Material(:semicon, 1_000.0, 1_000.0) # [Ω·m]
semicon_outer = Material(:semicon, 500.0, 1_000.0) # [Ω·m]
xlpe = Material(:insulator, 1.0e14, 2.77781, 1.59891) # [Ω·m]
lead = Material(:conductor, 2.14e-7, 1.0, 0.999983, 20.0, 0.004) # [Ω·m], [°C], [1/°C]
pe = Material(:insulator, 1.97e14, 2.51862) # [Ω·m]
pp = Material(:insulator, 1.0e15, 2.8, 1.1263) # [Ω·m]

# ## The stranded core
#
# The core holds 127 copper wires of 3.66 mm: a center wire and six courses of 6, 12, 18,
# 24, 30 and 36 wires, each course laid at lay ratio 15. In the cable the wires are
# compacted until they fill the circle of 127 wire areas. This release represents the
# compacted core with rectangular strands. Course `k` is a ring between two radii, and each
# of its `6k` rectangular strands has one wire's area, which the [`CircularPattern`](@ref)
# bends into an annular sector. The center wire fills the disk of one 3.66 mm wire, and the
# sectors of one wire's area tile the six rings of the conductor, the rectangular-strand
# construction of this release.

strand_radius = 3.66e-3 / 2 # strand_diameter / 2 [m]
strand_area = π * strand_radius^2 # [m²]
rings = [strand_radius * sqrt(1 + 3k * (k + 1)) for k in 0:6] # [m]
courses = [Lay(Helix(LayRatio(15.0)), # core_lay_ratio
               Repeat(CircularPattern(6k; r = Offset(0)),
                   Region(Symbol(:core_strands_, k),
                       Rectangle(strand_area / (rings[k + 1] - rings[k]),
                           rings[k + 1] - rings[k]),
                       core_metal)))
           for k in 1:6]
rings[end] # [m]

# ## The cable
#
# The center wire is a [`core`](@ref) of the wire's radius. The courses follow it in the
# same [`@terminal`](@ref) block. [`screen`](@ref) declares the two semicon screens and the
# water-blocking layer, [`sheath`](@ref) the lead sheath, and [`jacket`](@ref) the inner
# polyethylene and the outer jacket.

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

# Each course has one row in the table of lays, with the lay radius of its rectangular
# strands, the pitch of their helix and their overlength. The overlength is the length of a
# strand per unit length of cable, and at one lay ratio every course has the same
# overlength.

[(radius = row.radii[1], pitch = row.pitches[1], overlength = row.overlengths[1])
 for row in design.lays] # [m], [m], [dimensionless]

# ## The equivalent tubes
#
# `LineCableModels.flatten` reduces each terminal to one tube. The center wire and the six
# courses conduct in parallel. The center wire is straight, and a strand of a course is
# longer than the cable by its overlength, so it has a higher resistance per unit length of
# cable. The tube of the core takes a resistivity above that of copper, by the ratio in the
# second cell. Its relative permeability reproduces the geometric mean radius of the center
# wire and the rectangular strands.

blueprint = LineCableModels.flatten(design)
[(terminal = conductor.terminal, r_in = conductor.r_in, r_ex = conductor.r_ex,
     rho = conductor.material.rho, mu_r = conductor.material.mu_r)
 for conductor in blueprint.conductors] # [m], [m], [Ω·m], [dimensionless]

#-

blueprint.conductors[1].material.rho / core_metal.rho

# The dielectric layers between the core and the sheath form the interval of the core. The
# semicon screens enter it as dielectrics of relative permittivity 1000, in the series chain
# of potential coefficients, so they add little to the potential coefficient of the
# interval.

[(conductor = layer.conductor, r_in = layer.r_in, r_ex = layer.r_ex,
     eps_r = layer.material.eps_r, mu_r = layer.material.mu_r)
 for layer in blueprint.dielectrics] # [m], [m], [dimensionless], [dimensionless]

# The current of a stranded core follows the helices of its rectangular strands, so it flows
# around the axis like the current of a solenoid. That adds magnetic energy in the interval.
# The relative permeability of each interval, `mu_r_ins`, includes this solenoid factor on
# top of the mean permeability of its layers. The interval of the sheath has no lay, so its
# value is the mean of its layers.

blueprint.mu_r_ins

# ## The system
#
# The cables lie at `x = −1, 0, 1` m. Each terminal takes its own conductor index: the core
# and the sheath of cable `i` are conductors `2i − 1` and `2i`.

earth = homogeneous(rho = 1_000.0, eps_r = 10.0) # [Ω·m]
frequencies = 10.0 .^ range(-1, 7; length = 101) # [Hz]
placements = [Placement(x, -1.0) for x in (-1.0, 0.0, 1.0)] # [m]
system = build(CableSystem, fill(design, 3), placements;
    connections = [(core = 2i - 1, sheath = 2i) for i in 1:3],
    id = "cable_380kv_no_armour_ac_flat", line_length = 20_000.0) # [m]
problem = LineParametersProblem(system, earth; temperature = 70.0, frequencies) # [°C]
result = compute(problem, CoaxialFormulation())
result.terminals

# The bundle and Kron reductions are on by default. Here no two terminals share an index and
# none has the index zero, and the result keeps the six conductors. At the frequency nearest
# 50 Hz, the cells show the real and imaginary parts of `Z` in Ω/km and of `Y` in µS/km,
# rounded.

k = argmin(abs.(result.frequencies .- 50.0))
result.frequencies[k] # [Hz]

#-

round.(real(result.Z[:, :, k]) .* 1e3; digits = 4) # [Ω/km]

#-

round.(imag(result.Z[:, :, k]) .* 1e3; digits = 4) # [Ω/km]

#-

round.(real(result.Y[:, :, k]) .* 1e9; digits = 4) # [µS/km]

#-

round.(imag(result.Y[:, :, k]) .* 1e9; digits = 1) # [µS/km]

# ## Grounded sheaths
#
# A sheath bonded to the earth at both ends of the line has zero voltage. The index `0`
# declares it. The Kron reduction then removes the sheaths from the matrices, and the result
# keeps the three cores, with the effect of the sheath currents included.

grounded = build(CableSystem, fill(design, 3), placements;
    connections = [(core = i, sheath = 0) for i in 1:3],
    id = "cable_380kv_grounded_sheaths", line_length = 20_000.0) # [m]
grounded_problem = LineParametersProblem(grounded, earth; temperature = 70.0, frequencies)
grounded_result = compute(grounded_problem, CoaxialFormulation())
grounded_result.terminals

#-

grounded_result.Z[:, :, k] # [Ω/m]

# ## Ideal transposition
#
# In a flat formation the middle cable differs from the outer ones, so the matrix is not
# symmetric under a cyclic exchange of the phases. A transposed line exchanges the positions
# of the phases along its length. The option `transposition = true` replaces each matrix by
# the mean over the cyclic permutations of its rows and columns, so the entries of each
# cyclic diagonal are equal.

transposed = compute(grounded_problem, CoaxialFormulation(transposition = true))
transposed.Z[:, :, k] # [Ω/m]

#-

size(result.Z), length(result.frequencies)
