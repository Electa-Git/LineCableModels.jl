"""
$(TYPEDSIGNATURES)

Select the frequency-dependent cable model (2016). UniversalLineModel fits the model.
The `yc`, `h`, `passivity` and `dc` option sections pass through unchanged.
"""
description(::Type{<:Formula{:fdcm}}; compact::Bool=false) =
    compact ? "fdcm" : "Frequency-dependent cable model (2016)"

function formulation_options(::Expression{<:Formula{:fdcm}, typeof(fit)})
    return FormulationOptions((yc=(;), h=(;), passivity=(;), dc=(;)))
end

function formulation_options(::Expression{<:Formula{:fdcm}, typeof(fit)},
        supplied::FormulationOptions)
    isempty(setdiff(keys(supplied.data), (:yc, :h, :passivity, :dc))) ||
        throw(ArgumentError("unknown :fdcm formulation options"))
    validate(supplied.data, Formula)
    return supplied
end

:fdcm
