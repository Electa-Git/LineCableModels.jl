struct FEMOnelabSnapshot
    names::Vector{String}
    parameters::Dict{String, String}
end

const FEM_ONELAB_ROOT = "LineCableModels/FEM"

function _onelab_name(suffix::AbstractString)
    return string(FEM_ONELAB_ROOT, "/", suffix)
end

function _snapshot_onelab()
    names = String[name for name in gmsh.onelab.get_names("^LineCableModels/FEM/")]
    parameters = Dict{String, String}()
    for name in names
        parameters[name] = gmsh.onelab.get(name)
    end
    return FEMOnelabSnapshot(names, parameters)
end

function _restore_onelab(snapshot::FEMOnelabSnapshot)
    for name in gmsh.onelab.get_names("^LineCableModels/FEM/")
        gmsh.onelab.clear(name)
    end
    for name in snapshot.names
        gmsh.onelab.set(snapshot.parameters[name])
    end
    return nothing
end

function _number_parameter(
        name::String,
        value::Real;
        visible::Bool = false,
        read_only::Bool = true,
        label::Union{Nothing, String} = nothing,
        attributes = Dict{String, String}()
)
    parameter = Dict{String, Any}(
        "type" => "number",
        "name" => name,
        "values" => [Float64(value)],
        "visible" => visible,
        "readOnly" => read_only,
        "attributes" => attributes
    )
    label === nothing || (parameter["label"] = label)
    return parameter
end

function _string_parameter(
        name::String,
        value::AbstractString;
        visible::Bool = false,
        read_only::Bool = true,
        label::Union{Nothing, String} = nothing,
        attributes = Dict{String, String}()
)
    parameter = Dict{String, Any}(
        "type" => "string",
        "name" => name,
        "values" => [String(value)],
        "visible" => visible,
        "readOnly" => read_only,
        "attributes" => attributes
    )
    label === nothing || (parameter["label"] = label)
    return parameter
end

function _publish_parameters(parameters)
    gmsh.onelab.set(JSON3.write(parameters))
    return nothing
end

function _publish_transport!(
        run::FEMRun,
        model::FEMResolvedModel,
        model_data_path::String,
        mesh_path::String,
        formulation::LineCableModelsFEM
)
    parameters = Any[
        _string_parameter(_onelab_name("run_directory"), run.path),
        _string_parameter(
            _onelab_name("model_data_path"), model_data_path),
        _string_parameter(
            _onelab_name("mesh_path"), mesh_path),
        _number_parameter(
            _onelab_name("terminal_count"), length(model.terminal_ids)
        ),
        _number_parameter(
            _onelab_name("frequency_count"), length(model.problem.frequencies)
        ),
        _number_parameter(
            _onelab_name("field_maps"), formulation.execution.plot_field_maps
        ),
        _number_parameter(
            _onelab_name("completion_status"), 0),
        _number_parameter(
            _onelab_name("current_source_index"), 0),
        _number_parameter(
            _onelab_name("current_frequency_index"), 0)
    ]
    for (index, frequency) in enumerate(model.problem.frequencies)
        push!(parameters, _number_parameter(
            _onelab_name("frequency/$index"), frequency
        ))
    end
    _publish_parameters(parameters)
    return nothing
end

function _button_parameter(name::String, label::String, payload::String)
    action = _onelab_name("ui/action")
    command = "SetString[\"$action\", \"$payload\"];"
    return _string_parameter(
        name,
        command;
        visible = true,
        read_only = false,
        label,
        attributes = Dict(
            "Macro" => "GmshParseString",
            "Aspect" => "Button",
            "AutoCheck" => "1",
            "ServerAction" => ""
        )
    )
end

function _ui_state_labels(state::Symbol)
    mesh_state = if state in (:mesh_ready, :running, :results_ready)
        "ready"
    elseif state === :mesh_requested
        "generating"
    elseif state === :mesh_required
        "required"
    else
        "not generated"
    end
    solve_state = if state === :results_ready
        "completed"
    elseif state === :running
        "running"
    else
        "not run"
    end
    return (; mesh_state, solve_state)
end

function _publish_ui!(run::FEMRun, model::FEMResolvedModel, state::Symbol)
    frequency_summary = join(model.problem.frequencies, ", ") * " Hz"
    states = _ui_state_labels(state)
    parameters = Any[
        _string_parameter(
            _onelab_name("ui/problem"),
            model.problem.system.system_id;
            visible = true,
            label = "Problem"
        ),
        _string_parameter(
            _onelab_name("ui/frequencies"),
            frequency_summary;
            visible = true,
            label = "Frequencies"
        ),
        _number_parameter(
            _onelab_name("ui/cable_count"),
            length(model.problem.system.designs);
            visible = true,
            label = "Cables"
        ),
        _number_parameter(
            _onelab_name("ui/terminal_count"),
            length(model.terminal_ids);
            visible = true,
            label = "Terminals"
        ),
        _string_parameter(
            _onelab_name("ui/run_directory"),
            run.path;
            visible = true,
            label = "Run directory"
        ),
        _string_parameter(
            _onelab_name("ui/status"),
            replace(string(state), '_' => ' ');
            visible = true,
            label = "Status"
        ),
        _string_parameter(
            _onelab_name("ui/mesh_state"),
            states.mesh_state;
            visible = true,
            label = "Mesh state"
        ),
        _string_parameter(
            _onelab_name("ui/solve_state"),
            states.solve_state;
            visible = true,
            label = "Solve state"
        ),
        _string_parameter(
            _onelab_name("ui/action"), ""),
        _button_parameter(
            _onelab_name("ui/generate_mesh"), "Generate mesh", "generate_mesh"
        ),
        _button_parameter(
            _onelab_name("ui/run_model"), "Run model", "run_model")
    ]
    _publish_parameters(parameters)
    return nothing
end

function _set_ui_status(state::Symbol)
    states = _ui_state_labels(state)
    gmsh.onelab.set_string(
        _onelab_name("ui/status"), [replace(string(state), '_' => ' ')]
    )
    gmsh.onelab.set_string(
        _onelab_name("ui/mesh_state"), [states.mesh_state]
    )
    gmsh.onelab.set_string(
        _onelab_name("ui/solve_state"), [states.solve_state]
    )
    return state
end

function _take_ui_action()
    values = gmsh.onelab.get_string(_onelab_name("ui/action"))
    action = isempty(values) ? "" : only(values)
    isempty(action) || gmsh.onelab.set_string(_onelab_name("ui/action"), [""])
    return Symbol(action)
end

function _ui_transition(state::Symbol, action::Symbol, available::Bool = true)
    available || return state === :geometry_ready ? :closed_before_mesh :
           state === :mesh_ready ? :closed_before_solve : state
    action === Symbol("") && return state
    state === :geometry_ready && action === :generate_mesh && return :mesh_requested
    state === :geometry_ready && action === :run_model && return :mesh_required
    state === :mesh_ready && action === :generate_mesh && return :mesh_requested
    state === :mesh_ready && action === :run_model && return :solve_requested
    return state
end
