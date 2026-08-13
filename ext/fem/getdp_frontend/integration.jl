# Integration: defining integration methods

using .GetDPFrontend: comment, make_args

# IntegrationItemCase for nested Case blocks (e.g., Case { { Type Gauss; Case { ... } } })
mutable struct IntegrationItemCase
    items::Vector{Union{SimpleItem, Pair{String, IntegrationItemCase}}}
    comment::String

    function IntegrationItemCase(; comment="")
        new([], comment)
    end
end

function add!(case::IntegrationItemCase; kwargs...)
    item = SimpleItem(; kwargs...)
    push!(case.items, item)
    return case
end

function add_nested_case!(case::IntegrationItemCase; type::String, kwargs...)
    nested_case = IntegrationItemCase(; kwargs...)
    push!(case.items, type => nested_case)
    return nested_case
end

function code(case::IntegrationItemCase; indent_level::Int=0)
    indent = "  " ^ indent_level
    case_codes = []
    for item in case.items
        if item isa SimpleItem
            # Indent SimpleItem code (e.g., { GeoElement Point; NumberOfPoints 1; })
            push!(case_codes, indent * "  " * code(item))
        elseif item isa Pair
            type, nested_case = item
            # Start typed block (e.g., { Type Gauss; )
            push!(case_codes, indent * "  " * "{ Type $type;")
            # Recursively generate nested case code with increased indent
            nested_code = code(nested_case; indent_level=indent_level+1)
            push!(case_codes, nested_code)
            # Close typed block
            push!(case_codes, indent * "  " * "}")
        end
    end
    # Join case codes with newlines
    codes = join(case_codes, "\n")
    # Wrap in Case block
    _code = indent * "Case {\n" * codes * "\n" * indent * "}"
    if !isempty(case.comment)
        _code *= " " * GetDP.comment(case.comment)
    end
    return _code
end
# IntegrationItem for each { Name I1; ... } block
mutable struct IntegrationItem
    name::String
    comment::String
    cases::Vector{IntegrationItemCase}

    function IntegrationItem(name; comment=nothing)
        new(name, comment !== nothing ? comment : "", IntegrationItemCase[])
    end
end

function add!(item::IntegrationItem; kwargs...)
    case = IntegrationItemCase(; kwargs...)
    push!(item.cases, case)
    return case
end

function code(item::IntegrationItem; indent_level::Int=0)
    indent = "  " ^ indent_level
    case_codes = join([code(case; indent_level=indent_level+1) for case in item.cases], "\n")
    _code = indent * "{ Name $(item.name);\n" * case_codes * "\n" * indent * "}"
    if !isempty(item.comment)
        _code *= " " * GetDP.comment(item.comment)
    end
    return _code
end
# Integration struct
mutable struct Integration <: AbstractGetDPObject
    name::String
    content::String
    items::Vector{IntegrationItem}
    comment::Union{String,Nothing}

    function Integration()
        new("Integration", "", IntegrationItem[], nothing)
    end
end

function add!(integration::Integration, name::String; comment=nothing, kwargs...)
    item = IntegrationItem(name; comment)
    push!(integration.items, item)
    update_content!(integration)  # Update content when adding an item
    return item
end

function update_content!(integration::Integration)
    content = ""
    for item in integration.items
        content *= code(item) * "\n"
    end
    integration.content = rstrip(content)
end

function code(integration::Integration; indent_level::Int=0)
    indent = "  " ^ indent_level
    code_lines = [indent * "\nIntegration {"]
    if !isempty(integration.items)
        for item in integration.items
            for line in split(code(item; indent_level=indent_level+1), '\n')
                push!(code_lines, line)
            end
        end
    elseif !isempty(integration.content)
        for line in split(integration.content, '\n')
            push!(code_lines, indent * "  " * line)
        end
    end
    push!(code_lines, indent * "}")
    _code = join(code_lines, "\n")
    if integration.comment !== nothing
        return comment(integration.comment) * "\n" * _code
    else
        return _code
    end
end
function add_raw_code!(integration::Integration, raw_code, newline=true)
    integration.content = add_raw_code(integration.content, raw_code, newline)
end

function add_comment!(integration::Integration, comment_text, newline=true)
    integration.comment = comment_text
end
