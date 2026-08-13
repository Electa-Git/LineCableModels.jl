# PostOperation: defining post-operations

using .GetDPFrontend: add_raw_code, comment, make_args

"""
    OperationItem

An operation item (Print or Echo) in a post-operation.
"""
mutable struct OperationItem
    code::String
    comment::Union{String,Nothing}

    function OperationItem(code; comment=nothing)
        new(code, comment)
    end
end

function code(item::OperationItem)
    if item.comment !== nothing
        return comment(item.comment, newline=false) * "\n    " * item.code
    else
        return "    " * item.code
    end
end

"""
    POBase_

A collection of operations (Print, Echo, etc.) in a post-operation.
"""
mutable struct POBase_
    operations::Vector{OperationItem}
    raw_codes::Vector{String}
    comment::Union{String,Nothing}

    function POBase_(; comment=nothing)
        new(OperationItem[], String[], comment)
    end
end

function add_operation!(pobase::POBase_, code; comment=nothing)
    item = OperationItem(code; comment=comment)
    push!(pobase.operations, item)
    item
end

function add_raw_code!(pobase::POBase_, raw_code)
    push!(pobase.raw_codes, raw_code)
end

function code(pobase::POBase_)
    code_lines = String[]
    push!(code_lines, "Operation {")
    for operation in pobase.operations
        push!(code_lines, code(operation))
    end
    push!(code_lines, "}")
    join(code_lines, "\n")
end

"""
    PostOperationItem

An item in a post-operation section.
"""
mutable struct PostOperationItem
    id::String
    NameOfPostProcessing::String
    operation::POBase_
    raw_codes::Vector{String}
    comment::Union{String,Nothing}
    kwargs::Dict

    function PostOperationItem(id; NameOfPostProcessing, comment=nothing, kwargs...)
        new(id, NameOfPostProcessing, POBase_(), String[], comment, Dict(kwargs))
    end
end

function add_operation!(poitem::PostOperationItem, header="Operation"; comment=nothing)
    poitem.operation = POBase_(; comment=comment)
    poitem.operation
end

function code(poitem::PostOperationItem)
    code_lines = ["{ Name $(poitem.id); NameOfPostProcessing $(poitem.NameOfPostProcessing);"]
    if poitem.comment !== nothing
        push!(code_lines, "  $(comment(poitem.comment, newline=false))")
    end
    operation_code = code(poitem.operation)
    for line in split(operation_code, '\n')
        push!(code_lines, "    $line")
    end
    push!(code_lines, "}")
    for raw_code in poitem.raw_codes
        push!(code_lines, raw_code)
    end
    join(code_lines, "\n")
end

"""
    PostOperation

Defining post-operations.
"""
mutable struct PostOperation <: AbstractGetDPObject
    name::String
    content::String
    items::Vector{PostOperationItem}
    comment::Union{String,Nothing}
    indent::String
    raw_codes::Vector{String}

    function PostOperation()
        new("PostOperation", "", PostOperationItem[], nothing, " "^4, String[])
    end
end

"""
    add!(postoperation::PostOperation, id, NameOfPostProcessing; kwargs...)

Add a post-operation item.
"""
function add!(postoperation::PostOperation, id, NameOfPostProcessing; comment=nothing, kwargs...)
    item = PostOperationItem(id; NameOfPostProcessing=NameOfPostProcessing, comment=comment, kwargs...)
    push!(postoperation.items, item)
    postoperation.content = code(postoperation)
    item
end

"""
    add_raw_code!(postoperation::PostOperation, raw_code, newline=true)

Add raw code to the PostOperation object.
"""
function add_raw_code!(postoperation::PostOperation, raw_code, newline=true)
    push!(postoperation.raw_codes, raw_code * (newline ? "\n" : ""))
end

"""
    add_raw_code!(poitem::PostOperationItem, raw_code)

Add raw code to a PostOperationItem.
"""
function add_raw_code!(poitem::PostOperationItem, raw_code)
    push!(poitem.raw_codes, raw_code)
end

"""
    add_comment!(postoperation::PostOperation, comment_text, newline=true)

Add a comment to the PostOperation object.
"""
function add_comment!(postoperation::PostOperation, comment_text, newline=true)
    add_raw_code!(postoperation, comment(comment_text, newline=false), newline)
end

"""
    code(postoperation::PostOperation)

Generate GetDP code for a PostOperation object.
"""
function code(postoperation::PostOperation)
    code_lines = ["\nPostOperation {"]
    for raw_code in postoperation.raw_codes
        push!(code_lines, "  " * raw_code)
    end
    for item in postoperation.items
        item_code = code(item)
        for line in split(item_code, '\n')
            if !isempty(line)
                push!(code_lines, "  $line")
            end
        end
    end
    push!(code_lines, "}")
    if postoperation.comment !== nothing
        return comment(postoperation.comment) * "\n" * join(code_lines, "\n") * "\n"
    else
        return join(code_lines, "\n") * "\n"
    end
end
