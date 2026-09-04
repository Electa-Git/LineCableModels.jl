# These AST helpers are core composition machinery for `@gridspace`. They
# operate on exactly one struct and reject ambiguous inputs.
function _strip_escapes(expression)
    if expression isa Expr && expression.head === :escape
        return _strip_escapes(expression.args[1])
    elseif expression isa Expr
        return Expr(expression.head, map(_strip_escapes, expression.args)...)
    end
    return expression
end

function _struct_nodes(expression, nodes = Expr[])
    expression isa Expr || return nodes
    expression.head === :struct && push!(nodes, expression)
    for argument in expression.args
        _struct_nodes(argument, nodes)
    end
    return nodes
end

function _get_struct_node(expression)
    nodes = _struct_nodes(expression)
    length(nodes) == 1 || throw(ArgumentError(
        "macro input must contain exactly one struct definition; found $(length(nodes))",
    ))
    return only(nodes)
end

function _replace_struct(expression, old_struct, new_struct)
    expression === old_struct && return new_struct
    expression isa Expr || return expression
    return Expr(
        expression.head,
        map(
            argument -> _replace_struct(argument, old_struct, new_struct),
            expression.args
        )...
    )
end

function _parse_fields(struct_body)
    fields = NamedTuple[]
    clean_body = Any[]
    for argument in struct_body.args
        if argument isa LineNumberNode || argument isa String
            push!(clean_body, argument)
            continue
        end
        has_default = argument isa Expr && argument.head === :(=)
        field_expression = has_default ? argument.args[1] : argument
        default = has_default ? argument.args[2] : nothing
        if field_expression isa Symbol
            name = field_expression
            declared_type = nothing
        elseif field_expression isa Expr && field_expression.head === :(::)
            name = field_expression.args[1]
            declared_type = field_expression.args[2]
        else
            push!(clean_body, argument)
            continue
        end
        push!(clean_body, field_expression)
        push!(fields, (; name, declared_type, has_default, default))
    end
    isempty(fields) &&
        throw(ArgumentError("macro struct must declare at least one field"))
    return Tuple(fields), clean_body
end

function _extract_struct_name(struct_node)
    signature = struct_node.args[2]
    signature isa Expr && signature.head === :(<:) &&
        (signature = signature.args[1])
    signature isa Symbol && return signature
    signature isa Expr && signature.head === :curly &&
        return signature.args[1]
    throw(ArgumentError("malformed struct signature"))
end

function _rebuild_ast(expression, old_struct, new_struct, generated)
    replaced = _replace_struct(expression, old_struct, new_struct)
    return Expr(:block, replaced, generated...)
end

"""
$(SIGNATURES)

Retain the strict positional struct and add a keyword constructor that lifts
only explicit Grid or Gridspace fields into a `Gridspace`. With `Target`, the
generated space materialises that target instead of the vault struct itself.
"""
macro gridspace(arguments...)
    length(arguments) in (1, 2) ||
        throw(ArgumentError("@gridspace accepts a struct and optional target"))
    target_expression,
    expression = length(arguments) == 1 ? (nothing, arguments[1]) :
                 arguments
    raw = _strip_escapes(expression)
    struct_node = _get_struct_node(raw)
    struct_name = _extract_struct_name(struct_node)
    fields, clean_body = _parse_fields(struct_node.args[3])
    target = target_expression === nothing ? struct_name : target_expression

    clean_struct = Expr(
        :struct,
        struct_node.args[1],
        struct_node.args[2],
        Expr(:block, clean_body...)
    )
    keywords = Any[]
    for field in fields
        push!(
            keywords,
            field.has_default ?
            Expr(:kw, field.name, field.default) : field.name
        )
    end
    push!(keywords, Expr(:kw, :combine, QuoteNode(:product)))
    signature = Expr(:call, struct_name, Expr(:parameters, keywords...))
    values = Expr(:tuple, map(field -> field.name, fields)...)
    constructor = Expr(:function,
        signature,
        quote
            combine in (:product, :zip) ||
                throw(ArgumentError("combine must be :product or :zip"))
            values = $values
            if any(
                value -> value isa Union{
                    $(GlobalRef(@__MODULE__, :AbstractGrid)),
                    $(GlobalRef(@__MODULE__, :Gridspace))
                },
                values
            )
                grids = map(values) do value
                    value isa Union{
                        $(GlobalRef(@__MODULE__, :AbstractGrid)),
                        $(GlobalRef(@__MODULE__, :Gridspace))
                    } ? value : $(GlobalRef(@__MODULE__, :Grid))((value,))
                end
                return $(GlobalRef(@__MODULE__, :Gridspace)){$target}(
                    $target,
                    grids;
                    combine
                )
            end
            return $target(values...)
        end)
    return esc(_rebuild_ast(raw, struct_node, clean_struct, (constructor,)))
end
