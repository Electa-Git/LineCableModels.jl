"""
$(TYPEDEF)

Identify one scientific quantity without storing its value.
"""
struct QuantityTag{Q} end

QuantityTag(::Val{Q}) where {Q} = QuantityTag{Q}()
QuantityTag(::Type{QuantityTag{Q}}) where {Q} = QuantityTag{Q}()
_quantity_name(::QuantityTag{Q}) where {Q} = Q

"Return the physical quantity selected by a scientific function object."
function quantity end
