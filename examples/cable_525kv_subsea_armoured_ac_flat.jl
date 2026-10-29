# # 525 kV subsea armored cables in flat formation
#
# Each of these three cables has a solid copper core, semicon screens, XLPE insulation, a
# water-blocking layer, a lead sheath, a polyethylene jacket, a bedding and an armor of 68
# steel wires under an outer jacket. The cables lie 1 m deep and 1 m apart in a flat
# formation, in earth of 1000 Ω·m. The tutorial computes the line parameters of the nine
# conductors, then of the three cores with the sheaths and armors grounded, and the effect
# of a prescribed longitudinal wavenumber on the nine conductors.

using LineCableModels

# ## A core sized by its resistance
#
# The design gives the direct-current resistance of the core at 20 °C, 0.0072 Ω/km, and its
# diameter, 60 mm. A solid core of that diameter has this resistance when its resistivity is
# the resistance per unit length times its area.

core_r20 = 0.0072 / 1_000 # core_r20 [Ω/km], in [Ω/m]
core_rho = core_r20 * π * (0.06 / 2)^2 # from core_diameter [m], in [Ω·m]

# ## Materials
#
# The water-blocking layer takes the material of the outer semicon screen. The steel of the
# armor has a relative permeability of 10.

copper = Material(:conductor, core_rho, 1.0, 1.0, 20.0, 0.00393) # [Ω·m], [°C], [1/°C]
semicon_inner = Material(:semicon, 1_000.0, 1_000.0) # [Ω·m]
semicon_outer = Material(:semicon, 500.0, 1_000.0) # [Ω·m]
xlpe = Material(:insulator, 1.0e14, 2.4) # [Ω·m]
lead = Material(:conductor, 2.14e-7, 1.0, 0.999983, 20.0, 0.004) # [Ω·m], [°C], [1/°C]
pe = Material(:insulator, 1.97e14, 2.5) # [Ω·m]
pp = Material(:insulator, 1.0e15, 2.8) # [Ω·m]
steel = Material(:conductor, 1.38e-7, 1.0, 10.0, 20.0, 0.0045) # [Ω·m], [°C], [1/°C]

# ## The cable
#
# The terminals are the core, the lead sheath and the armor.

design = @cable "cable_525kv_subsea_armoured" begin
    @terminal :core begin
        core(copper; r = 0.06 / 2, tag = :core_conductor) # core_diameter / 2 [m]
    end
    screen(semicon_inner; t = 2.0e-3, tag = :core_semicon_inner) # inner_semicon_thickness [m]
    insulation(xlpe; t = 26.0e-3, tag = :core_insulation) # insulation_thickness [m]
    screen(semicon_outer; t = 1.8e-3, tag = :core_semicon_outer) # outer_semicon_thickness [m]
    screen(semicon_outer; t = 0.6e-3, tag = :core_water_blocking) # water_blocking_thickness [m]
    @terminal :sheath begin
        sheath(lead; t = 3.2e-3, tag = :sheath_lead) # lead_sheath_thickness [m]
    end
    jacket(pe; t = 3.3e-3, tag = :sheath_inner_pe) # inner_pe_thickness [m]
    bedding(pp; t = 1.0e-3, tag = :sheath_bedding) # bedding_thickness [m]
    @terminal :armor begin
        armor(steel; shape = Circle(6.0e-3 / 2), n = 68, # armor_wire_diameter / 2 [m]
            lay = LayRatio(10.019920548139554), tag = :armor_wires) # armor_lay_ratio
    end
    jacket(pp; t = 6.6e-3, tag = :armor_jacket) # jacket_thickness [m]
end
r_ex(design.outer) # [m]

# `LineCableModels.flatten` gives one tube per terminal. The tube of the solid core keeps
# the resistivity of the core, and the tube of the armor takes the resistance of its wires
# along their helices.

blueprint = LineCableModels.flatten(design)
[(terminal = conductor.terminal, r_in = conductor.r_in, r_ex = conductor.r_ex,
     rho = conductor.material.rho, mu_r = conductor.material.mu_r)
 for conductor in blueprint.conductors] # [m], [m], [Ω·m], [dimensionless]

# ## The system
#
# The core, the sheath and the armor of cable `i` are conductors `3i − 2`, `3i − 1` and
# `3i`.

earth = homogeneous(rho = 1_000.0, eps_r = 10.0) # [Ω·m]
frequencies = 10.0 .^ range(-1, 7; length = 101) # [Hz]
placements = [Placement(x, -1.0) for x in (-1.0, 0.0, 1.0)] # [m]
system = build(CableSystem, fill(design, 3), placements;
    connections = [(core = 3i - 2, sheath = 3i - 1, armor = 3i) for i in 1:3],
    id = "cable_525kv_subsea_armoured_ac_flat", line_length = 20_000.0) # [m]
problem = LineParametersProblem(system, earth; temperature = 70.0, frequencies) # [°C]
result = compute(problem, CoaxialFormulation())
result.terminals

# No two terminals share an index and none has the index zero, so the bundle and Kron
# reductions leave the nine conductors as computed. At the frequency nearest 50 Hz, the
# cells show the real and imaginary parts of `Z` in Ω/km and of `Y` in µS/km, rounded.

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

# ## Sheaths and armors grounded
#
# With the sheath and the armor of each cable at index `0`, the Kron reduction removes them,
# and the three cores remain.

grounded = build(CableSystem, fill(design, 3), placements;
    connections = [(core = i, sheath = 0, armor = 0) for i in 1:3],
    id = "cable_525kv_grounded_screens", line_length = 20_000.0) # [m]
grounded_problem = LineParametersProblem(grounded, earth; temperature = 70.0, frequencies)
grounded_result = compute(grounded_problem, CoaxialFormulation())
grounded_result.terminals

#-

grounded_result.Z[:, :, k] # [Ω/m]

#-

grounded_result.Y[:, :, k] # [S/m]

# ## A prescribed longitudinal wavenumber
#
# The unified earth-return formula takes the longitudinal wavenumber Γ of the fields along
# the line as an input, zero by default. A nonzero Γ describes fields that propagate along
# the line: here, slightly slower than light in vacuum, one value per frequency. The cell
# computes the nine-conductor system again with it.

Γ = im .* 1.01 .* 2π .* frequencies ./ 299792458 # [1/m]
gamma_result = compute(problem, CoaxialFormulation(Γ = Γ))
gamma_result.Z[1:3, 1:3, k] # [Ω/m]

# At each frequency, `change` is the largest change of `Z`, relative to the largest entry
# computed with Γ = 0. The cell gives its value at 50 Hz, its largest value over the scan
# and the frequency of that value.

change = [maximum(abs, gamma_result.Z[:, :, j] - result.Z[:, :, j]) /
          maximum(abs, result.Z[:, :, j]) for j in eachindex(frequencies)]
(at_50Hz = change[k], largest = maximum(change),
    frequency = frequencies[argmax(change)]) # [Hz]

#-

size(result.Z), length(result.frequencies)
