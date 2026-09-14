"""
    LineCableModels.DataModel.BaseParams

Calculate cable-geometry equivalences, reference-state resistance, and
export-oriented dielectric reductions. Each function promotes its inputs with
Base numeric rules.
"""
module BaseParams

export equivalent_alpha, parallel, helix, strip_resistance
export tubular_resistance, wire_coordinates
export strand_gmr, tubular_gmr, equivalent_gmr, equivalent_mu
export shunt_capacitance, shunt_conductance, series_shunt_admittance
export solenoid_factor, equivalent_rho, equivalent_eps
export equivalent_conductivity

using DocStringExtensions: TYPEDSIGNATURES

include("geometry.jl")
include("resistance.jl")
include("inductance.jl")
include("dielectrics.jl")

end
