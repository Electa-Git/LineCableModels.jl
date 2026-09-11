"""
$(TYPEDSIGNATURES)

Write one fully materialised line-parameter problem as versioned JSON.

# Arguments

- `problem`: Completed scalar line-parameter problem.

# Keywords

- `file_name`: Destination JSON file.

# Returns

- Absolute path of the written file.
"""
function export_data(
        ::Val{:json},
        problem::Engine.LineParametersProblem;
        file_name::AbstractString
)
    path = _json_path(String(file_name))
    open(path, "w") do io
        #! explicit-imports: off
        JSON3.pretty(io, _json_document(problem); allow_inf = true)
        #! explicit-imports: on
    end
    return abspath(path)
end

"""
$(TYPEDSIGNATURES)

Read one fully materialised line-parameter problem from versioned JSON.

# Arguments

- `LineParametersProblem`: Requested result type.

# Keywords

- `file_name`: Source JSON file.

# Returns

- A validated scalar [`Engine.LineParametersProblem`](@ref).
"""
function import_data(
        ::Val{:json},
        ::Type{Engine.LineParametersProblem};
        file_name::AbstractString
)
    isfile(file_name) || throw(ArgumentError(
        "line-parameter problem file not found: '$(_display_path(file_name))'",
    ))
    document = _read_document(file_name, PROBLEM_SCHEMA)
    problem = deserialize_value(_required(document, "root", PROBLEM_SCHEMA))
    problem isa Engine.LineParametersProblem || throw(ArgumentError(
        "line-parameter problem document did not decode as LineParametersProblem",
    ))
    return problem
end
