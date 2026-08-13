# Base GetDP objects

using .GetDPFrontend: add_raw_code, comment, make_args

"""
    GetDPObject

Base GetDP object that other objects inherit from.
"""
abstract type AbstractGetDPObject end

mutable struct GetDPObject <: AbstractGetDPObject
    name::String
    content::String
    comment::Union{String,Nothing}
    indent::String

    function GetDPObject(name="Group", content="", comment=nothing)
        indent = " "^4  # Same as Python's len(self.name)
        new(name, content, comment, indent)
    end
end

"""
    code(obj::GetDPObject)

Generate GetDP code for the object.
"""
function code(obj::GetDPObject)
    code_lines = String[]
    push!(code_lines, "$(obj.name){")

    for line in split(obj.content, '\n')
        push!(code_lines, line)
    end

    push!(code_lines, "}")

    code_str = join(code_lines, "\n" * obj.indent) * "\n"

    if obj.comment !== nothing
        code_str = comment(obj.comment) * "\n" * code_str
    end

    return code_str
end




"""
    add_raw_code!(obj::GetDPObject, raw_code, newline=true)

Add raw code to the object's content.
"""
function add_raw_code!(obj::GetDPObject, raw_code, newline=true)
    obj.content = add_raw_code(obj.content, raw_code, newline)
end

"""
    add_comment!(obj::GetDPObject, comment_text, newline=true)

Add a comment to the object's content.
"""
function add_comment!(obj::GetDPObject, comment_text, newline=true)
    add_raw_code!(obj, comment(comment_text; newline=false), newline)
end

"""
    SimpleItem

A simple item in a GetDP object.
"""
mutable struct SimpleItem
    comment::String
    code::String

    function SimpleItem(; comment=nothing, kwargs...)
        code_str = " { "

        for (k, v) in kwargs
            code_str *= " $k " * make_args(v, sep=",") * ";"
        end

        code_str *= " } "

        if comment !== nothing
            code_str *= GetDP.comment(comment)
        end

        new("", code_str)
    end
end

function code(item::SimpleItem)
    return item.code
end

"""
    Base_

Base class for GetDP objects with items.
"""
abstract type AbstractBase_ end

mutable struct Base_ <: AbstractBase_
    comment::String
    code::String
    items::Vector{SimpleItem}
    _code0::String

    function Base_(header; comment=nothing, kwargs...)
        code_str = header * " { "

        if comment !== nothing
            code_str *= GetDP.comment(comment)
        end

        code_str *= " \n        }"

        new("", code_str, SimpleItem[], code_str)
    end
end

"""
    add!(base::Base_, args...; kwargs...)

Add a SimpleItem to a Base_ object.
"""
function add!(base::Base_, args...; kwargs...)
    item = SimpleItem(args...; kwargs...)
    s = base.code
    n = 10
    base.code = s[1:end-n] * "\n       " * item.code * s[end-n+1:end]
    push!(base.items, item)
    return item
end

"""
    ObjectItem

An object item in a GetDP object.
"""
mutable struct ObjectItem
    Name::String
    comment::String
    cases::Vector{SimpleItem}

    function ObjectItem(Name; comment=nothing)
        new(Name, comment !== nothing ? comment : "", SimpleItem[])
    end
end

"""
    code(obj::ObjectItem)

Generate GetDP code for the object item.
"""
function code(item::ObjectItem)
    case_codes = join([code(case) for case in item.cases], "\n")
    _code = "{ Name $(item.Name); Case {\n" * case_codes * "\n} }"
    if !isempty(item.comment)
        _code *= " " * GetDP.comment(item.comment)
    end
    return _code
end

"""
    CaseItem_

A case item in a GetDP object.
"""
mutable struct CaseItem_
    Region::Any
    code::String

    function CaseItem_(Region; comment=nothing, kwargs...)
        code_str = "{" * " Region $(Region);"

        for (k, v) in kwargs
            code_str *= " $k " * make_args(v, sep=",") * ";"
        end

        code_str *= " } "

        if comment !== nothing
            code_str *= GetDP.comment(comment)
        end

        new(Region, code_str)
    end
end

"""
    Case_

A case in a GetDP object.
"""
mutable struct Case_
    Name::Union{String,Nothing}
    comment::String
    code::String
    case_items::Vector{CaseItem_}

    function Case_(; Name=nothing, comment=nothing, kwargs...)
        case_name = Name === nothing ? "" : Name
        code_str = "Case $(case_name) "
        code_str *= "{ "

        if comment !== nothing
            code_str *= GetDP.comment(comment)
        end

        code_str *= " \n     }"

        new(Name, "", code_str, CaseItem_[])
    end
end

"""
    add!(case::Case_, args...; kwargs...)

Add a CaseItem_ to a Case_ object.
"""
function add!(case::Case_, args...; kwargs...)
    case_item = CaseItem_(args...; kwargs...)
    s = case.code
    n = 7
    case.code = s[1:end-n] * "\n       " * case_item.code * s[end-n+1:end]
    push!(case.case_items, case_item)
    return case_item
end

"""
    add!(obj::ObjectItem, args...; kwargs...)

Add a Case_ to an ObjectItem.
"""
function add!(obj::ObjectItem, args...; kwargs...)
    obj._code0 = code(obj)
    case = Case_(args...; kwargs...)
    push!(obj.cases, case)
    return case
end
