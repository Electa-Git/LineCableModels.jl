"""
$(TYPEDSIGNATURES)

Select the universal line model (1999). UniversalLineModel fits the model.
The `yc`, `h` and `passivity` sections supply fitting controls.
The `dc` section selects DC correction with targets from the physical line.
"""
description(::Type{<:Formula{:ulm}}; compact::Bool=false) =
    compact ? "ulm" : "Universal line model (1999)"

function formulation_options(::Expression{<:Formula{:ulm}, typeof(fit)})
    return FormulationOptions((yc=(;), h=(;), passivity=(;), dc=(;), line=LineParametersFormulation()))
end

function formulation_options(::Expression{<:Formula{:ulm}, typeof(fit)},
        supplied::FormulationOptions)
    isempty(setdiff(keys(supplied.data), (:yc, :h, :passivity, :dc, :line))) ||
        throw(ArgumentError("unknown :ulm formulation options"))
    if haskey(supplied.data, :line)
        supplied.data.line::LineParametersFormulation
    end
    dc = get(supplied.data, :dc, (;))::NamedTuple
    if !isempty(dc)
        isempty(setdiff(keys(dc), (:method, :k))) ||
            throw(ArgumentError("unknown :ulm dc options"))
        get(dc, :method, nothing) in (:add_pole, :functional_form) ||
            throw(ArgumentError("dc.method must be :add_pole or :functional_form"))
        if haskey(dc, :k)
            dc.method == :add_pole || throw(ArgumentError("dc.k requires :add_pole"))
            k = dc.k::Real
            isfinite(k) && 0 < k < 1 || throw(ArgumentError("dc.k must be in (0, 1)"))
        end
    end
    return supplied
end

function validate(line::Nothing, problem::LineModelProblem, selected::Formula{:ulm})
    validate(line, problem)
    isempty(get(selected.options.data, :dc, (;))) || throw(ArgumentError(
        "dc requires LineModelProblem(segment; line=...) and its line formulation option"))
    get(get(selected.options.data, :h, (;)), :td, false) && throw(ArgumentError(
        "td=true requires LineModelProblem(segment; line=...)"))
    return line
end

function validate(line::LineParametersProblem, problem::LineModelProblem,
        selected::Formula{:ulm})
    validate(line, problem)
    recorded = details(problem.segment).data
    inputs = get(recorded, :inputs, nothing)
    selections = get(recorded, :formulations, nothing)
    if inputs === nothing || selections === nothing
        @warn "line consistency cannot be checked: the segment has no input or formulation record"
        return line
    end
    validate(completed_inputs(line), problem, inputs, (:inputs,))
    line_formulation = get(selected.options.data, :line, LineParametersFormulation())
    declaration = NamedTuple(line_formulation)
    validate(declaration.formulas, problem, selections.formulas, (:formulations, :formulas))
    validate(declaration.options, problem, selections.options, (:formulations, :options))
    return line
end

:ulm
