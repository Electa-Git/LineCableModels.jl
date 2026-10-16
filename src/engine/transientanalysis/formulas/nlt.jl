"""
$(TYPEDSIGNATURES)

Select the numerical Laplace transform of the line two-port.
"""
description(::Type{<:Formula{:nlt}}; compact::Bool=false) =
    compact ? "NLT" : "Numerical Laplace transform of the line two-port"

function formulation_options(::Expression{<:Formula{:nlt}, typeof(respond!)})
    return FormulationOptions((line=Formulation(), N=4096, Tfactor=1.25, taper=:hanning))
end

function formulation_options(expression::Expression{<:Formula{:nlt}, typeof(respond!)},
        supplied::FormulationOptions)
    defaults = formulation_options(expression).data
    isempty(setdiff(keys(supplied.data), keys(defaults))) ||
        throw(ArgumentError("unknown :nlt formulation options"))
    controls = merge(defaults, supplied.data)
    controls.line::LineParametersFormulation
    controls.N::Int
    controls.taper in (:hanning, :blackman, :lanczos) ||
        throw(ArgumentError("taper must be :hanning, :blackman or :lanczos"))
    return FormulationOptions(controls)
end

:nlt
