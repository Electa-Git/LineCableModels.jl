"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:unified` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? description(Formula{:unified};compact=true) : "Default routing to :unified"

Formula{:default}(; kwargs...) = Formula{:unified}(; kwargs...)

:default
