"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:constant` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :constant"

Formula(::Val{:default}; kwargs...) = Formula(Val(:constant); kwargs...)

:default
