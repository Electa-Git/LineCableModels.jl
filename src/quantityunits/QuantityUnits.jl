"""
    LineCableModels.QuantityUnits

Define physical-quantity tags, display-unit descriptions, labels, and metric
prefix conversion factors.

# Public actions

- `quantity` resolves semantic quantity identity.
- `default_unit` and `display_unit` select storage and presentation units.
- `scale_factor` converts values between compatible unit descriptions.
"""
module QuantityUnits

using Base: @kwdef
using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES
import ..LineCableModels: R, L, C

export Unit, Units, units, get_label, get_symbol, get_exp
export QuantityTag, quantity, default_unit, display_unit, scale_factor

include("units.jl")
include("quantities.jl")
include("accessors.jl")
include("scaling.jl")
include("definitions.jl")

end # module QuantityUnits
