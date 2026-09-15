"""
$(TYPEDSIGNATURES)

Route the default selection to the explicit `:ametani1980` implementation.
No numerical equation is owned by `:default`.
"""
description(::Type{<:Formula{:default}}; compact::Bool=false) =
    compact ? "Default" : "Default routing to :ametani1980"

Formula(::Val{:default}; kwargs...) = Formula(Val(:ametani1980); kwargs...)

:default
