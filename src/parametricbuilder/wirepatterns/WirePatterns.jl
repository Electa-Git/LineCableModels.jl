"""
    LineCableModels.ParametricBuilder.WirePatterns

Estimate deterministic wire patterns and retain ranked candidates with their
geometric packing limits.
"""
module WirePatterns

using DocStringExtensions: TYPEDEF, TYPEDFIELDS, TYPEDSIGNATURES

import ...LineCableModels: nominal

"Return the maximum wire count admitted by one estimate geometry."
function maxfill end

export WireEstimate, make_stranded, make_screened

include("types.jl")
include("gauges.jl")
include("stranded.jl")
include("screened.jl")

end # module WirePatterns
