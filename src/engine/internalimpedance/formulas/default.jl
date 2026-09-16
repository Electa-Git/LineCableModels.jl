"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:schelkunoff1934` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :schelkunoff1934"

Formula(::Val{:default}; kwargs...) = Formula(Val(:schelkunoff1934); kwargs...)

:default
