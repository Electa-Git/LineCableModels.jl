"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:nlt` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :nlt"

Formula{:default}(; kwargs...) = Formula{:nlt}(; kwargs...)

:default
