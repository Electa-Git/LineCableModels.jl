"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:lossless` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :lossless"

Formula{:default}(; kwargs...) = Formula{:lossless}(; kwargs...)

:default
