"Record the consumed FEM material laws and selected field assumptions."
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
    return merge((
        schema_version = 4,
        selections,
        assumptions = (
            impedance = "Axial current-driven A_z/u_r finite-element equations",
            admittance = _quasi_full(formulation.options.physics) ?
                "Coupled first-order Maxwell A_z/A_t/phi equations; axial current supplies normalized leakage; vertical path voltage includes A_t/Gamma; Y = inv(P)" :
                "Scalar electrodynamic Helmholtz equation in surrounding media; equipotential terminals with unit transverse-current excitation; Y = inv(P)",
            earth = "Horizontal air and one semi-infinite soil; soil constitutive properties evaluated at each frequency",
            propagation = _quasi_full(formulation.options.physics) ?
                "Gamma -> 0 with A_t/Gamma and phi/Gamma retained; conduction and displacement; one coupled factorization" :
                "Gamma = 0; medium diffusion and displacement retained; independent Z/P blocks in one factorization",
            semicon_domain = "Passive material region, without electrical terminal ownership",
            enclosure = "Supported enclosures are represented by their material and terminal domains"
        )
    ),NamedTuple(formulation))
end
