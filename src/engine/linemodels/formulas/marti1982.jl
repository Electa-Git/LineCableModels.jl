"""
$(TYPEDSIGNATURES)

Select the constant-transformation modal line model (1982).
The options specify the frequency [Hz] for the transformation and
the Bode tolerance [dB]. This selection throws an error until an implementation exists.
"""
description(::Type{<:Formula{:marti}}; compact::Bool=false) =
    compact ? "Marti" : "Constant-transformation modal line model (1982)"

function formulation_options(::Expression{<:Formula{:marti}, typeof(fit)})
    return FormulationOptions((transformation=(frequency=nothing,),
        bode=(tolerance_db=3.0,)))
end

function formulation_options(expression::Expression{<:Formula{:marti}, typeof(fit)},
        supplied::FormulationOptions)
    defaults = formulation_options(expression).data
    isempty(setdiff(keys(supplied.data), keys(defaults))) ||
        throw(ArgumentError("unknown :marti formulation options"))
    sections = map(keys(defaults)) do name
        explicit = get(supplied.data, name, (;))::NamedTuple
        section = getproperty(defaults, name)
        isempty(setdiff(keys(explicit), keys(section))) ||
            throw(ArgumentError("unknown :marti $name options"))
        return merge(section, explicit)
    end
    return FormulationOptions(NamedTuple{keys(defaults)}(sections))
end

:marti
