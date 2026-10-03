"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:linear` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :linear"

Formula(::Val{:default}; kwargs...) = Formula(Val(:linear); kwargs...)

:default
