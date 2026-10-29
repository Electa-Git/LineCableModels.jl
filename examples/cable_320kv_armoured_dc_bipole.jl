# # 320 kV armored cables of a bipole
#
# Each of the two cables of this bipole has the stranded copper core and the insulation
# system of the 380 kV design, a lead sheath, a bedding and an armor of 68 steel wires under
# a jacket. The cables lie 1 m deep and 1 m apart, in earth of 1000 Ω·m. The tutorial
# computes the line parameters of the six conductors, then with the sheath and the armor of
# each cable bundled, and reads the matrices from which the computation builds them.

using LineCableModels

# ## Materials
#
# The armor wires are galvanized steel with a relative permeability of 36.6.

core_metal = Material(:conductor, 2.3853e-8, 1.0, 1.01663, 20.0, 0.00393) # [Ω·m]
semicon_inner = Material(:semicon, 1_000.0, 1_000.0) # [Ω·m]
semicon_outer = Material(:semicon, 500.0, 1_000.0) # [Ω·m]
xlpe = Material(:insulator, 1.0e14, 2.77781, 1.59891) # [Ω·m]
lead = Material(:conductor, 2.14e-7, 1.0, 0.999983, 20.0, 0.004) # [Ω·m], [°C], [1/°C]
pe = Material(:insulator, 1.97e14, 2.51862) # [Ω·m]
pp = Material(:insulator, 1.0e15, 2.8, 1.1263) # [Ω·m]
steel = Material(:conductor, 1.7475e-7, 1.0, 36.6326, 20.0, 0.0045) # [Ω·m], [°C], [1/°C]

# ## The cable
#
# The core repeats that of the 380 kV tutorial, a center wire and six courses of rectangular
# strands, each bent into an annular sector of one 3.66 mm wire's area.

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

# The builder [`armor`](@ref) lays `n` round wires in one course on the outline of the part
# inside it, here the bedding, at the lay ratio `lay`. A fill of an insulator completes the
# layer between the wires. The bedding is as thick as the 68 wires need to fit side by side
# in their course.

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

# The last row of the table of lays is the armor's course: the radius of the wires' centers,
# the pitch of their helix and their overlength.

armor_lay = last(design.lays)
(radius = armor_lay.radii[1], pitch = armor_lay.pitches[1],
    overlength = armor_lay.overlengths[1]) # [m], [m], [dimensionless]

# ## The equivalent tubes
#
# `LineCableModels.flatten` reduces the armor to a tube that spans its course. The tube
# takes the resistance of the 68 wires along their helices. Its relative permeability is the
# value that gives the tube the geometric mean radius of the wires, so it differs from that
# of steel.

blueprint = LineCableModels.flatten(design)
[(terminal = conductor.terminal, r_in = conductor.r_in, r_ex = conductor.r_ex,
     rho = conductor.material.rho, mu_r = conductor.material.mu_r)
 for conductor in blueprint.conductors] # [m], [m], [Ω·m], [dimensionless]

# ## The system
#
# The cables lie at `x = ±0.5` m. The core, the sheath and the armor of cable `i` are
# conductors `3i − 2`, `3i − 1` and `3i`.

earth = homogeneous(rho = 1_000.0, eps_r = 10.0) # [Ω·m]
frequencies = 10.0 .^ range(-1, 7; length = 101) # [Hz]
placements = [Placement(-0.5, -1.0), Placement(0.5, -1.0)] # [m]
system = build(CableSystem, [design, design], placements;
    connections = [(core = 3i - 2, sheath = 3i - 1, armor = 3i) for i in 1:2],
    id = "cable_320kv_armoured_dc_bipole", line_length = 150_000.0) # [m]
problem = LineParametersProblem(system, earth; temperature = 70.0, frequencies) # [°C]
result = compute(problem, CoaxialFormulation())
result.terminals

# No two terminals share an index and none has the index zero. The bundle and Kron
# reductions, on by default, then leave the six conductors as computed. At the frequency
# nearest 50 Hz, the cells show the real and imaginary parts of `Z` in Ω/km and of `Y` in
# µS/km, rounded.

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

# ## Sheath and armor in one bundle
#
# Terminals that share a conductor index are bonded: they have one voltage, and their
# currents add. With the default `bundle = true`, the computation merges them into one
# conductor, labeled with both terminals.

bundled = build(CableSystem, [design, design], placements;
    connections = [(core = 2i - 1, sheath = 2i, armor = 2i) for i in 1:2],
    id = "cable_320kv_bundled_screens", line_length = 150_000.0) # [m]
bundled_problem = LineParametersProblem(bundled, earth; temperature = 70.0, frequencies)
bundled_result = compute(bundled_problem, CoaxialFormulation())
bundled_result.terminals

#-

round.(real(bundled_result.Z[:, :, k]) .* 1e3; digits = 4) # [Ω/km]

# ## The trace
#
# The computation of [`LCMBackend`](@ref LineCableModels.LCMBackend) takes the option
# `trace`. With `trace = true` the result keeps the four matrices from which it builds `Z`
# and `Y` at each frequency:
#
# - `Zin` is the series impedance in Ω/m inside each cable from its conductors to its outer
#   surface, by conductor.
# - `Pin` is the matrix of potential coefficients in m/F inside each cable, by conductor.
# - `Zg` is the earth-return impedance in Ω/m from the outer surfaces of the cables, by
#   cable.
# - `Pg` is the matrix of earth-return potential coefficients in m/F, by cable.

traced = compute(problem, CoaxialFormulation(); trace = true)
map(size, traced.trace)

#-

traced.trace.Zg[:, :, k] # [Ω/m]

# The entry of `Zg` for two cables adds to every pair of their conductors. Without a
# reduction that acts, `Z` is `Zin` plus that expansion of `Zg`. The cell gives the largest
# difference between the two in Ω/m.

expanded = [traced.trace.Zg[cld(i, 3), cld(j, 3), k] for i in 1:6, j in 1:6] # [Ω/m]
maximum(abs, traced.trace.Zin[:, :, k] + expanded - traced.Z[:, :, k]) # [Ω/m]

#-

size(result.Z), length(result.frequencies)
