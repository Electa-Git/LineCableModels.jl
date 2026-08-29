"""
    LineCableModels.DataModel.BaseParams

Calculate cable geometry and base-state RLCG properties. Each function
promotes its inputs with Base numeric rules.
"""
module BaseParams

export equivalent_alpha, parallel, helix, strip_resistance, temperature_factor
export tubular_resistance, tubular_inductance, wire_coordinates
export trefoil_inductance, strand_gmr, tubular_gmr, equivalent_gmr, equivalent_mu
export shunt_capacitance, shunt_conductance, series_shunt_admittance
export solenoid_factor, equivalent_rho, equivalent_eps, loss_tangent
export equivalent_conductivity

using DocStringExtensions: TYPEDSIGNATURES

include("geometry.jl")
include("resistance.jl")
include("inductance.jl")
include("dielectrics.jl")

end
