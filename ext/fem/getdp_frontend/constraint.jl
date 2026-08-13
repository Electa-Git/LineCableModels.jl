# Constraint: defining constraints

using .GetDPFrontend: add_raw_code, comment, make_args

"""
    ConstraintItem

Represents an individual constraint with cases and optional parameters.
"""
mutable struct ConstraintItem
    name::String
    type::String
    cases::Vector{Dict{String, Any}}
    comment::Union{String, Nothing}
    for_loop::Union{Nothing, Tuple{String, String}}  # (index, range) for For loops

    function ConstraintItem(name; type="Assign", comment=nothing)
        new(name, type, Dict{String, Any}[], comment, nothing)
    end
end

"""
    Constraint

Defining constraints.
"""
mutable struct Constraint <: AbstractGetDPObject
    name::String
    constraints::Vector{ConstraintItem}
    comment::Union{String, Nothing}
    indent::String

    function Constraint()
        new("Constraint", ConstraintItem[], nothing, " "^4)
    end
end

"""
    add!(constraint::Constraint, name; type="Assign", comment=nothing)

Add a constraint to the Constraint object.
Returns the ConstraintItem for adding cases.
"""
function add!(constraint::Constraint, name; type="Assign", comment=nothing)
    item = ConstraintItem(name, type=type, comment=comment)
    push!(constraint.constraints, item)
    return item
end

"""
    assign!(constraint::Constraint, name; comment=nothing)

Add an Assign-type constraint to the Constraint object (shortcut for add!).
Returns the ConstraintItem for adding cases.
"""
function assign!(constraint::Constraint, name; comment=nothing)
    add!(constraint, name, type="Assign", comment=comment)
end

"""
    case!(item::ConstraintItem, region; value=nothing, time_function=nothing, comment=nothing)

Add a case to a ConstraintItem with specified region and parameters.
"""
function case!(item::ConstraintItem, region; value=nothing, time_function=nothing, comment=nothing)
    case_dict = Dict{String, Any}("Region" => region)
    if value !== nothing
        case_dict["Value"] = value
    end
    if time_function !== nothing
        case_dict["TimeFunction"] = time_function
    end
    if comment !== nothing
        case_dict["Comment"] = comment
    end
    push!(item.cases, case_dict)
end

"""
    for_loop!(item::ConstraintItem, index, range)

Set a For loop for the ConstraintItem (e.g., For k In {1:3}).
"""
function for_loop!(item::ConstraintItem, index, range)
    item.for_loop = (index, range)
    push!(item.cases, Dict("LoopStart" => true))  # Add marker for loop start
    return item
end

"""
    add_raw_code!(constraint::Constraint, raw_code, newline=true)

Add raw code to the Constraint object.
"""
function add_raw_code!(constraint::Constraint, raw_code, newline=true)
    # Store raw code as a special ConstraintItem with no cases
    item = ConstraintItem("RawCode")
    item.cases = [Dict("Raw" => raw_code)]
    push!(constraint.constraints, item)
end

"""
    add_comment!(constraint::Constraint, comment_text, newline=true)

Add a comment to the Constraint object.
"""
function add_comment!(constraint::Constraint, comment_text, newline=true)
    add_raw_code!(constraint, comment(comment_text, newline=false), newline)
end

"""
    code(constraint::Constraint)

Generate GetDP code for a Constraint object.
"""
function code(constraint::Constraint)
    code_lines = String[]
    push!(code_lines, "\nConstraint{")

    for item in constraint.constraints
        # Handle raw code items
        if item.name == "RawCode" && length(item.cases) == 1 && haskey(item.cases[1], "Raw")
            push!(code_lines, "  " * item.cases[1]["Raw"])
            continue
        end

        # Add item comment if present
        if item.comment !== nothing
            push!(code_lines, "  " * comment(item.comment, newline=false))
        end

        # Start constraint block
        push!(code_lines, "  { Name $(item.name); Type $(item.type);")
        push!(code_lines, "    Case {")

        # Split cases into pre-loop and loop cases
        loop_start_idx = findfirst(c -> haskey(c, "LoopStart"), item.cases)
        pre_loop_cases = loop_start_idx === nothing ? item.cases : item.cases[1:loop_start_idx-1]
        loop_cases = loop_start_idx === nothing ? Dict{String, Any}[] : item.cases[loop_start_idx+1:end]

        # Pre-loop cases
        for case in pre_loop_cases
            if haskey(case, "Comment") && case["Comment"] !== nothing
                push!(code_lines, "      " * comment(case["Comment"], newline=false))
            end
            if !haskey(case, "Region") || case["Region"] === ""
                continue
            end
            line = "      { Region $(case["Region"])"
            if haskey(case, "Value")
                line *= "; Value $(case["Value"])"
            end
            if haskey(case, "TimeFunction")
                line *= "; TimeFunction $(case["TimeFunction"])"
            end
            line *= "; }"
            push!(code_lines, line)
        end

        # Loop cases
        if item.for_loop !== nothing
            index, range = item.for_loop
            push!(code_lines, "      For $(index) In {$(range)}")
            for case in loop_cases
                if haskey(case, "Comment") && case["Comment"] !== nothing
                    push!(code_lines, "        " * comment(case["Comment"], newline=false))
                end
                if !haskey(case, "Region") || case["Region"] === ""
                    continue
                end
                line = "        { Region $(case["Region"])"
                if haskey(case, "Value")
                    line *= "; Value $(case["Value"])"
                end
                if haskey(case, "TimeFunction")
                    line *= "; TimeFunction $(case["TimeFunction"])"
                end
                line *= "; }"
                push!(code_lines, line)
            end
            push!(code_lines, "      EndFor")
        end

        push!(code_lines, "    }")
        push!(code_lines, "  }")
    end

    push!(code_lines, "}")
    if constraint.comment !== nothing
        return comment(constraint.comment) * "\n" * join(code_lines, "\n")
    else
        return join(code_lines, "\n")
    end
end
