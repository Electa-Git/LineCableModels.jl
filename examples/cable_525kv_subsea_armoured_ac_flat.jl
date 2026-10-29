# # 525 kV subsea armored cables in flat formation
#=
This example demonstrates how to model a three-phase system composed of 525 kV subsea
armored cables. Each cable consists of a solid copper core, inner and outer semicon screens,
XLPE insulation, a water-blocking layer, a lead sheath, a polyethylene jacket, a bedding and
an armor of 68 steel wires under an outer jacket. The cables are buried at 1 m depth in a
flat formation, with a horizontal spacing of 1 m, in a homogeneous earth of 1000 Ω·m.

This example covers:

1. Obtaining the resistivity of a solid core from its specified DC resistance.
2. Computing the line parameters of the nine conductors of the system.
3. Eliminating the sheaths and armors of a system with both grounded at both ends.
4. Evaluating the effect of a prescribed longitudinal wavenumber on the series impedance.
=#

# Load the package:
using LineCableModels

#=
## Core resistivity from the DC resistance

The cable data specify the DC resistance of the core at 20 °C, 0.0072 Ω/km, and its
diameter, 60 mm. The resistivity of a solid core with this resistance is the product of the
resistance per unit length and the cross-section of the core.
=#

# DC resistance per unit length and equivalent resistivity of the core:
core_r20 = 0.0072 / 1_000 # core_r20 [Ω/km], in [Ω/m]
core_rho = core_r20 * π * (0.06 / 2)^2 # from core_diameter [m], in [Ω·m]

#=
## Materials

The water-blocking layer is made of the same material as the outer semicon screen. The steel
of the armor has a relative permeability of 10.
=#

# Conductor, semicon, insulation and armor materials:
copper = Material(:conductor, core_rho, 1.0, 1.0, 20.0, 0.00393) # [Ω·m], [°C], [1/°C]
semicon_inner = Material(:semicon, 1_000.0, 1_000.0) # [Ω·m]
semicon_outer = Material(:semicon, 500.0, 1_000.0) # [Ω·m]
xlpe = Material(:insulator, 1.0e14, 2.4) # [Ω·m]
lead = Material(:conductor, 2.14e-7, 1.0, 0.999983, 20.0, 0.004) # [Ω·m], [°C], [1/°C]
pe = Material(:insulator, 1.97e14, 2.5) # [Ω·m]
pp = Material(:insulator, 1.0e15, 2.8) # [Ω·m]
steel = Material(:conductor, 1.38e-7, 1.0, 10.0, 20.0, 0.0045) # [Ω·m], [°C], [1/°C]

#=
## Cable design

The cable has three terminals, namely the core, the lead sheath and the armor.
=#

# Build the cable design and inspect its outer radius:
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

#=
## Equivalent coaxial representation

The function `LineCableModels.flatten` returns one tubular conductor per terminal. The tube
of the solid core retains the resistivity of the core, whereas the tube of the armor takes
the resistance of the wires along their helices.
=#

# Equivalent conductors, with radii in m and resistivity in Ω·m:
blueprint = LineCableModels.flatten(design)
[(terminal = conductor.terminal, r_in = conductor.r_in, r_ex = conductor.r_ex,
     rho = conductor.material.rho, mu_r = conductor.material.mu_r)
 for conductor in blueprint.conductors] # [m], [m], [Ω·m], [dimensionless]

#=
## Three-phase system in flat formation

The core, sheath and armor of cable `i` are connected to conductors `3i − 2`, `3i − 1` and
`3i`, respectively.
=#

# Homogeneous earth model and frequency range:
earth = homogeneous(rho = 1_000.0, eps_r = 10.0) # [Ω·m]
frequencies = 10.0 .^ range(-1, 7; length = 101) # [Hz]

# Cable positions and cable system for a line 20 km long:
placements = [Placement(x, -1.0) for x in (-1.0, 0.0, 1.0)] # [m]
system = build(CableSystem, fill(design, 3), placements;
    connections = [(core = 3i - 2, sheath = 3i - 1, armor = 3i) for i in 1:3],
    id = "cable_525kv_subsea_armoured_ac_flat", line_length = 20_000.0) # [m]

# Line parameters problem at 70 °C:
problem = LineParametersProblem(system, earth; temperature = 70.0, frequencies) # [°C]

# Compute the line parameters and inspect the row labels:
result = compute(problem, CoaxialFormulation())
result.terminals

#=
## Line parameters at power frequency

Since no two terminals share an index and no terminal has the index zero, the bundle and
Kron reductions leave the nine conductors unchanged. The matrices at the frequency closest
to 50 Hz are presented below, with `Z` in Ω/km and `Y` in µS/km, rounded.
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
## Sheaths and armors grounded at both ends

The sheath and the armor of each cable are connected to the index 0, and the Kron reduction
eliminates them from the matrices. The resulting system retains the three cores.
=#

# Cable system with grounded sheaths and armors:
grounded = build(CableSystem, fill(design, 3), placements;
    connections = [(core = i, sheath = 0, armor = 0) for i in 1:3],
    id = "cable_525kv_grounded_screens", line_length = 20_000.0) # [m]
grounded_problem = LineParametersProblem(grounded, earth; temperature = 70.0, frequencies)

# Compute the line parameters and inspect the row labels:
grounded_result = compute(grounded_problem, CoaxialFormulation())
grounded_result.terminals

# Series impedance of the three cores:
grounded_result.Z[:, :, k] # [Ω/m]

# Shunt admittance of the three cores:
grounded_result.Y[:, :, k] # [S/m]

#=
## Prescribed longitudinal wavenumber

The unified earth-return formula takes the longitudinal wavenumber Γ of the fields along the
line as an input, equal to zero by default. A nonzero Γ represents fields that propagate
along the line, in this case slightly slower than light in vacuum, with one value per
frequency. The system of nine conductors is computed again with this value of Γ.
=#

# Longitudinal wavenumber and the resulting series impedance of the cores:
Γ = im .* 1.01 .* 2π .* frequencies ./ 299792458 # [1/m]
gamma_result = compute(problem, CoaxialFormulation(Γ = Γ))
gamma_result.Z[1:3, 1:3, k] # [Ω/m]

#=
At each frequency, `change` is the largest absolute change of `Z`, relative to the largest
entry computed with Γ = 0. The value at 50 Hz, the largest value over the frequency range and
the frequency at which it occurs are computed below.
=#

# Relative change of the series impedance due to Γ:
change = [maximum(abs, gamma_result.Z[:, :, j] - result.Z[:, :, j]) /
          maximum(abs, result.Z[:, :, j]) for j in eachindex(frequencies)]
(at_50Hz = change[k], largest = maximum(change),
    frequency = frequencies[argmax(change)]) # [Hz]

# Size of the impedance array and number of frequencies:
size(result.Z), length(result.frequencies)

#=
## Conclusion

This example has demonstrated how to:

1. Obtain the resistivity of a solid core from its specified DC resistance.
2. Compute the line parameters of a three-phase system of armored cables.
3. Eliminate grounded sheaths and armors with the Kron reduction.
4. Evaluate the effect of a prescribed longitudinal wavenumber on the series impedance.
=#
