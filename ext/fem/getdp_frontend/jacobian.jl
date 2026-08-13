# Jacobian: defining jacobians

using .GetDPFrontend: add_raw_code, comment, make_args

# function code(item::SimpleItem)
#     return item.code
# end

mutable struct CaseItem
    regions::Vector{Dict{String, Any}}
    comment::Union{String, Nothing}
    
    function CaseItem(; comment=nothing)
        new([], comment)
    end
end

function code(case::CaseItem)
    result = "Case {\n"
    for region in case.regions
        result *= "    { Region $(region["Region"]); Jacobian $(region["Jacobian"]);}"
        if haskey(region, "comment") && region["comment"] !== nothing
            result *= " " * comment(region["comment"])
        end
        result *= "\n"
    end
    result *= "}"
    return result
end

function add!(case::CaseItem; Region, Jacobian, comment=nothing)
    push!(case.regions, Dict("Region" => Region, "Jacobian" => Jacobian, "comment" => comment))
    return case
end

function add!(item::ObjectItem; kwargs...)
    case = SimpleItem(; kwargs...)
    push!(item.cases, case)
    return item
end
function update_code!(item::ObjectItem)
    item._code = "{ Name $(item.Name); \n  "
    
    if !isempty(item.cases)
        item._code *= "Case {\n"
        for case in item.cases
            for region in case.regions
                item._code *= "    { Region $(region["Region"]); Jacobian $(region["Jacobian"]);"
                item._code *= "}"
                if haskey(region, "comment") && region["comment"] !== nothing
                    item._code *= " " * comment(region["comment"])
                end
                item._code *= "\n"
            end
        end
        item._code *= "}"
    end
    
    item._code *= "}"
end

# function code(item::ObjectItem)
#     case_codes = join([code(case) for case in item.cases], "\n")
#     _code = "{ Name $(item.Name); Case {\n" * case_codes * "\n} }"
#     if !isempty(item.comment)
#         _code *= " " * GetDP.comment(item.comment)
#     end
#     return _code
# end
# Jacobian struct
mutable struct Jacobian <: AbstractGetDPObject
    name::String
    content::String
    items::Vector{ObjectItem}
    comment::Union{String,Nothing}

    function Jacobian()
        new("Jacobian", "", ObjectItem[], nothing)
    end
end

function add!(jacobian::Jacobian, name::String; kwargs...)
    item = ObjectItem(name; kwargs...)
    push!(jacobian.items, item)
    update_content!(jacobian)
    return item
end

function VolSphShell(; Rint, Rext, center_X=nothing, center_Y=nothing, center_Z=nothing, power=nothing, inv_inf=nothing)
    params = [Rint, Rext]
    for param in [center_X, center_Y, center_Z, power, inv_inf]
        if param !== nothing
            push!(params, param)
        end
    end
    param_str = join(params, ", ")
    return "VolSphShell{$param_str}"
end
function update_content!(jacobian::Jacobian)
    content = ""
    for item in jacobian.items
        content *= code(item) * "\n"
    end
    jacobian.content = rstrip(content)
end

function code(jac::Jacobian)
    code_lines = ["\nJacobian{"]
    if !isempty(jac.items)
        for item in jac.items
            for line in split(code(item), '\n')
                push!(code_lines, "  " * line)
            end
        end
    elseif !isempty(jac.content)
        for line in split(jac.content, '\n')
            push!(code_lines, "  " * line)
        end
    end
    push!(code_lines, "}")
    if jac.comment !== nothing
        return comment(jac.comment) * "\n" * join(code_lines, "\n")
    else
        return join(code_lines, "\n")
    end
end

function add_raw_code!(jacobian::Jacobian, raw_code, newline=true)
    jacobian.content = add_raw_code(jacobian.content, raw_code, newline)
end

function add_comment!(jacobian::Jacobian, comment_text, newline=true)
    jacobian.comment = comment_text
end
