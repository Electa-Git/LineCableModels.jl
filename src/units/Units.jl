"""
    LineCableModels.Units

Define the physical-quantity identities and units used by calculation results,
plots, and reports.

# Public actions

- `quantity` maps a scientific selector to its physical meaning.
- `native_unit` and `display_unit` select calculation and presentation units.
- `scale_factor` converts compatible unit expressions.
- `label` and `symbol` format unit and quantity metadata.
"""
module Units

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: R, L, C

export Unit, UnitExpr, Quantity, units
export quantity, native_unit, display_unit, scale_factor, label, symbol

include("units.jl")
include("quantities.jl")
include("accessors.jl")
include("scaling.jl")
include("definitions.jl")

end # module Units
