"""
$(TYPEDEF)

Represent the scientific meaning selected by `Q` without storing numerical
values or presentation metadata.
"""
struct Quantity{Q} end

"""
    quantity(selector::Function)
    quantity(selector::Function, transform::Function)

Return the typed [`Quantity`](@ref) identity registered for a scientific
selector or selector/transform pair.

Unsupported selectors have no fallback method.
"""
function quantity end
