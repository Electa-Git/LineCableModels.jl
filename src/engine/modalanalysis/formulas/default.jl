"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:chrysochos2014` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :chrysochos2014"

Formula{:default}(; kwargs...) = Formula{:chrysochos2014}(; kwargs...)

:default
