"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:ulm` implementation.
The `:default` identifier selects an equation without defining one.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :ulm"

Formula{:default}(; kwargs...) = Formula{:ulm}(; kwargs...)

:default
