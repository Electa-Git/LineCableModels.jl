"""
$(TYPEDSIGNATURES)

Select time-domain simulation of a fitted line model.
"""
description(::Type{<:Formula{:td}}; compact::Bool=false) =
    compact ? "Time domain" : "Time-domain simulation of a fitted line model"

function formulation_options(::Expression{<:Formula{:td}, typeof(respond!)})
    return FormulationOptions()
end

function formulation_options(::Expression{<:Formula{:td}, typeof(respond!)},
        supplied::FormulationOptions)
    isempty(supplied.data) || throw(ArgumentError(":td has no formulation options yet"))
    return supplied
end

:td
