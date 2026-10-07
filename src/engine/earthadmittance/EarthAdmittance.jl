"""
    LineCableModels.Engine.EarthAdmittance

Define earth-return admittance recipes, numerical primitives, and
formula-owned frequency functors.

# Dependencies

$(IMPORTS)

"""
module EarthAdmittance
import ...Commons: FormulationOptions, formulas

# Export public API
export Formula, formula_id, earth_potential_coefficient, formulas

# Module-specific dependencies
#! explicit-imports: off
# These abbreviations are expanded in this module docstring and included files.
using DocStringExtensions: IMPORTS, TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
#! explicit-imports: on
import ...Earth: EquivalentHomogeneous
import ..Engine: EarthAdmittanceFormulation, formula_id
#! explicit-imports: off
# Explicitly included equations share these physical and numerical operations.
import ...LineCableModels: validate
import ...Commons: Functor
import ..Engine: EarthPair
import ..Engine: EarthPlan, earth!, same_physical_state, layer_index
using ...Earth: EarthModel
import ...Commons: initialize_buffers
import ..Engine: computation_type, EarthImpedanceFormulation, special_besselix,
                 SpectralIntegral, integrate
using LinearAlgebra: lu!, ldiv!
import ..EarthImpedance
import ...LineCableModels: FormulaDefinition, Expression, nominal
import ..Engine: description, conductivity
import ..Engine: formulation_options
import ..Engine: AirVoltageSpectrum, earth_spectral_term, earth_spectral_value,
                 earth_spectral_points!, earth_contour_angle, earth_direct,
                 outgoing_root, bessel_i0m1, bessel_current_ratio
using ...Commons: vacuum_permittivity
#! explicit-imports: on

public source_potential_coefficient, earth!

include("interface.jl")

#! explicit-imports: off
const FORMULAS = (
    include("formulas/unified.jl"),
    include("formulas/default.jl"),
    include("formulas/ideal.jl"),
    include("formulas/pollaczek1926.jl"),
    include("formulas/wise1948.jl"),
    include("formulas/xue2018.jl")
)
#! explicit-imports: on

"""
Return registered earth-admittance identities, including unimplemented equations.
"""
formulas(::Type{<:Formula}) = FORMULAS

end # module EarthAdmittance
