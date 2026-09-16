"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:bottommost` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :bottommost"

Formula(::Val{:default}; kwargs...) = Formula(Val(:bottommost); kwargs...)

:default
