# # Solid 1000 mm² single-core cable
#
# This cable has one solid aluminum conductor of 1000 mm² under 8.3 mm of XLPE insulation.
# It lies 1 m deep in earth of 100 Ω·m. The tutorial builds the cable, places it in a line
# of 1000 m and computes its series impedance and shunt admittance per unit length from 0.1
# Hz to 10 MHz.

using LineCableModels

# ## Materials
#
# A [`Material`](@ref) takes its class first: `:conductor`, `:insulator` or `:semicon`. The
# positional values that follow are the resistivity, the relative permittivity, the relative
# permeability, the temperature at which the resistivity applies and the temperature
# coefficient of the resistivity. Values left out take their defaults.

aluminum = Material(:conductor, 2.8264e-8, 1.0, 1.000022, 20.0, 0.00429) # [Ω·m]
xlpe = Material(:insulator, 1.97e14, 2.5) # [Ω·m]

# ## The cable
#
# The [`@cable`](@ref) block lists the parts of a cable from the inside out. The
# [`@terminal`](@ref) block names the conductive regions that form one terminal of the
# cable, here its core. The builder [`core`](@ref) declares a solid circular conductor by
# its radius. The builder [`insulation`](@ref) declares a layer by its thickness, laid on
# the outline of the part inside it. The last value of the cell is the outer radius of the
# cable in meters.

core_radius = sqrt(1000.0e-6 / π) # from core_cross_section [m²], in [m]
design = @cable "solid_1000mm2" begin
    @terminal :core begin
        core(aluminum; r = core_radius, tag = :core_metal)
    end
    insulation(xlpe; t = 8.3e-3, tag = :core_insulation) # insulation_thickness [m]
end
r_ex(design.outer) # [m]

# A coaxial formulation reads each terminal as a tube. `LineCableModels.flatten` returns
# that view of the design: here one solid tube of the core's radius and one dielectric layer
# around it. Radii are in meters and the resistivity in Ω·m.

blueprint = LineCableModels.flatten(design)
[(terminal = conductor.terminal, r_in = conductor.r_in, r_ex = conductor.r_ex,
     rho = conductor.material.rho) for conductor in blueprint.conductors]

#-

[(conductor = layer.conductor, r_in = layer.r_in, r_ex = layer.r_ex,
     eps_r = layer.material.eps_r) for layer in blueprint.dielectrics]

# ## The system and the problem
#
# A [`CableSystem`](@ref) places cables in one cross-section, with `y` negative below the
# surface of the earth. Each terminal takes a conductor index, here 1 for the core. The
# [`LineParametersProblem`](@ref) adds the earth, the temperature of the conductors and the
# frequencies of the scan.

earth = homogeneous(rho = 100.0, eps_r = 10.0) # [Ω·m]
system = build(CableSystem, [design], [Placement(0.0, -1.0)]; # [m]
    connections = [(core = 1,)], id = "solid_1000mm2_single", line_length = 1000.0) # [m]
frequencies = 10.0 .^ range(-1, 7; length = 101) # [Hz]
problem = LineParametersProblem(system, earth; temperature = 20.0, frequencies) # [°C]
length(problem.frequencies)

# ## The computation
#
# A [`CoaxialFormulation`](@ref) selects one formula for each physical contribution, by
# default the formulas of this release. Its result, from [`compute`](@ref), is a
# [`LineParameters`](@ref) record that stores `Z` in Ω/m and `Y` in S/m by row, column and
# frequency, with one label per row.

result = compute(problem, CoaxialFormulation())
result.terminals

# The frequencies of the scan do not include 50 Hz exactly. The cells below read the
# frequency nearest to it, then `Z` and `Y` at that frequency.

k = argmin(abs.(result.frequencies .- 50.0))
result.frequencies[k] # [Hz]

#-

result.Z[:, :, k] # [Ω/m]

#-

result.Y[:, :, k] # [S/m]

# The real part of `Z` is the resistance per unit length `R` and its imaginary part is `ωL`,
# with `ω = 2πf`. The imaginary part of `Y` is `ωC`. `Z` is the impedance of the loop that
# the core forms with the earth, so `R` and `L` include the earth return. `C` comes from the
# potential coefficients of the insulation and of the earth around the cable.

ω = 2π * result.frequencies[k] # [rad/s]
(R = real(result.Z[1, 1, k]), L = imag(result.Z[1, 1, k]) / ω,
    C = imag(result.Y[1, 1, k]) / ω) # [Ω/m], [H/m], [F/m]

# ## Resistance and frequency
#
# The resistance of the core to direct current at 20 °C is its resistivity over its area. At
# 0.1 Hz the loop resistance is close to that value. As the frequency rises, the current in
# the core moves toward the outer surface of the core (the skin effect), and the return
# current in the earth flows closer to the cable. Both raise the loop resistance, at 50 Hz
# and up to 10 MHz.

(R_dc = aluminum.rho / (π * core_radius^2), R_first = real(result.Z[1, 1, 1]),
    R_50Hz = real(result.Z[1, 1, k]), R_last = real(result.Z[1, 1, end])) # [Ω/m]

# The tutorial of the 320 kV bipole separates the matrices of the cable from those of the
# earth return.

size(result.Z), length(result.frequencies)
