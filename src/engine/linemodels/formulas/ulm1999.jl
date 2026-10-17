"""
$(TYPEDSIGNATURES)

Select the universal line model (1999). UniversalLineModel fits the model.
The `yc`, `h`, `passivity` and `dc` option sections pass through unchanged.
"""
description(::Type{<:Formula{:ulm}}; compact::Bool=false) =
    compact ? "ulm" : "Universal line model (1999)"

function formulation_options(::Expression{<:Formula{:ulm}, typeof(fit)})
    return FormulationOptions((yc=(;), h=(;), passivity=(;), dc=(;)))
end

function formulation_options(::Expression{<:Formula{:ulm}, typeof(fit)},
        supplied::FormulationOptions)
    isempty(setdiff(keys(supplied.data), (:yc, :h, :passivity, :dc))) ||
        throw(ArgumentError("unknown :ulm formulation options"))
    validate(supplied.data, Formula)
    return supplied
end

:ulm
