"""
$(TYPEDSIGNATURES)

Select time-domain simulation of a fitted line model.
"""
description(::Type{<:Formula{:td}}; compact::Bool=false) =
    compact ? "Time domain" : "Time-domain simulation of a fitted line model"

function formulation_options(::Expression{<:Formula{:td}, typeof(respond!)})
    return FormulationOptions((interp=:single, enforce_sign=false, initial_solve=false))
end

function formulation_options(expression::Expression{<:Formula{:td}, typeof(respond!)},
        supplied::FormulationOptions)
    defaults = formulation_options(expression).data
    isempty(setdiff(keys(supplied.data), keys(defaults))) ||
        throw(ArgumentError("unknown :td formulation options"))
    controls = merge(defaults, supplied.data)
    controls.interp::Symbol
    controls.enforce_sign::Bool
    controls.initial_solve::Bool
    controls.interp in (:single, :double) ||
        throw(ArgumentError("interp must be :single or :double"))
    return FormulationOptions(controls)
end

validate(problem::TransientAnalysisProblem{<:LineModel}, ::Formula{:td}) = problem

:td
