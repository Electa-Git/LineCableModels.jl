"Record the consumed FEM material laws and fixed field assumptions."
function formulation_record(formulation::LineCableModelsFEM)
    # Hook descriptions identify the supplied callable, without claiming that
    # arbitrary Julia closures can be reconstructed from a saved record.
    Selection = NamedTuple{(:identifier, :parameters, :options, :hooks, :replayable),
        Tuple{Symbol, NamedTuple, NamedTuple, NamedTuple, Bool}}
    Selections = NamedTuple{keys(formulation.methods),
        NTuple{length(formulation.methods), Union{Nothing, Selection}}}
    selections = Selections(map(formulation.methods) do selected
        selected === nothing && return nothing
        hooks = map(selected.hooks) do hook
            (type=string(typeof(hook)), representation=repr(hook), replayable=false)
        end
        Selection((identifier=LineCableModels.formula_id(selected),
            parameters=selected.parameters, options=selected.options, hooks,
            replayable=isempty(selected.hooks)))
    end)
    return (
        schema_version = 3,
        selections,
        assumptions = (
            impedance = "Fixed quasi-TEM finite-element field equations",
            earth = "Horizontal air and one semi-infinite soil; soil constitutive properties evaluated at each frequency",
            propagation = "Fixed package quasi-TEM propagation approximation",
            semicon_domain = "Passive material region, without electrical terminal ownership",
            enclosure = "Supported enclosures are represented by their material and terminal domains"
        )
    )
end
