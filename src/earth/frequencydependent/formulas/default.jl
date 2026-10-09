"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:constant` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :constant"

Formula{:default}(; kwargs...) = Formula{:constant}(; kwargs...)

:default
