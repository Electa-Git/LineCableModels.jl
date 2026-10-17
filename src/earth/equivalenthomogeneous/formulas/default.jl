"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:bottommost` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :bottommost"

Formula{:default}(; kwargs...) = Formula{:bottommost}(; kwargs...)

:default
