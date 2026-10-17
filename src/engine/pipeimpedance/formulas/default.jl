"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:none` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :none"

Formula{:default}(; kwargs...) = Formula{:none}(; kwargs...)

:default
