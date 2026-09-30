"""
    LineCableModels.ParametricBuilder.WirePatterns

Estimate deterministic wire patterns and retain ranked patterns with their
geometric packing limits.
"""
module WirePatterns

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

import ...LineCableModels: nominal

"""
Return the maximum wire count admitted by one estimate geometry.
"""
function maxfill end

export WireEstimate, estimate_stranding, estimate_screen
public HexaPattern, ScreenPattern

include("types.jl")
include("gauges.jl")
include("stranded.jl")
include("screened.jl")

end # module WirePatterns
