"""
$(TYPEDEF)

Identify one scientific quantity without storing its value.
"""
struct Quantity{Q} end

"Return the physical quantity selected by a scientific function object."
function quantity end
