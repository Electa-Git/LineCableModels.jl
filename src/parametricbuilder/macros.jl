# These AST helpers are core composition machinery for `@gridspace` and
# `@relax`. They operate on exactly one struct and reject ambiguous inputs.
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
        elseif field_expression isa Expr && field_expression.head === :(::)
            name = field_expression.args[1]
        else
            push!(clean_body, argument)
            continue
        end
        push!(clean_body, field_expression)
        push!(fields, (; name, has_default, default))
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

function _extract_type_parameters(struct_node)
    signature = struct_node.args[2]
    signature isa Expr && signature.head === :(<:) &&
        (signature = signature.args[1])
    signature isa Expr && signature.head === :curly || return ()
    return Tuple(signature.args[2:end])
end

function _extract_parametric_supertype(struct_node)
    signature = struct_node.args[2]
    signature isa Expr && signature.head === :(<:) || return nothing
    supertype = signature.args[2]
    supertype isa Expr && supertype.head === :curly || return nothing
    return supertype.args[1]
end

function _type_parameter_name(parameter)
    parameter isa Symbol && return parameter
    parameter isa Expr && parameter.head in (:<:, :>:) &&
        return parameter.args[1]
    return nothing
end

function _rebuild_ast(expression, old_struct, new_struct, generated)
    replaced = _replace_struct(expression, old_struct, new_struct)
    return Expr(:block, replaced, generated...)
end

recast(::Type{T}, value::Number) where {T <: Real} = convert(T, value)
function recast(::Type{T}, values::AbstractArray) where {T <: Real}
    map(value -> recast(T, value), values)
end
recast(::Type{T}, values::Tuple) where {T <: Real} = map(value -> recast(T, value), values)
function recast(::Type{T}, values::NamedTuple) where {T <: Real}
    map(value -> recast(T, value), values)
end
recast(::Type{<:Real}, value) = value

_relax_eltype(value::Number) = typeof(value)
_relax_eltype(values::AbstractArray{<:Number}) = eltype(values)
_relax_eltype(values::Tuple) = isempty(values) ? nothing : _promoted_numeric_type(values)
_relax_eltype(value) =
    try
        candidate = eltype(typeof(value))
        candidate <: Number ? candidate : nothing
    catch
        nothing
    end

function _promoted_numeric_type(values::Tuple)
    types = filter(!isnothing, map(_relax_eltype, values))
    isempty(types) && throw(ArgumentError(
        "@relax constructor received no numeric fields to promote",
    ))
    return promote_type(types...)
end

"""
$(SIGNATURES)

Retain the strict positional struct and add a keyword constructor that lifts
only explicit Grid/spec fields into a `Gridspace`. With `Target`, the
generated space materializes that target instead of the vault struct itself.
"""
macro gridspace(arguments...)
    length(arguments) in (1, 2) ||
        throw(ArgumentError("@gridspace accepts a struct and optional target"))
    target_expression, expression = length(arguments) == 1 ? (nothing, arguments[1]) :
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

"""
$(SIGNATURES)

Add promoted numeric construction and conversion methods to a parametric
structure.
"""
macro relax(expression)
    raw = _strip_escapes(expression)
    struct_node = _get_struct_node(raw)
    struct_name = _extract_struct_name(struct_node)
    fields, _ = _parse_fields(struct_node.args[3])
    parameters = _extract_type_parameters(struct_node)
    isempty(parameters) &&
        throw(ArgumentError("@relax requires a parametric struct"))
    scalar_parameter = _type_parameter_name(first(parameters))
    scalar_parameter === nothing && throw(ArgumentError(
        "@relax could not identify the leading scalar type parameter",
    ))

    field_names = map(field -> field.name, fields)
    values = Expr(:tuple, field_names...)
    recasts = [:($(GlobalRef(@__MODULE__, :recast))(promoted_type, $(field.name)))
               for field in fields]
    target_recasts = [:($(GlobalRef(@__MODULE__, :recast))(TargetScalar, value.$(field.name)))
                      for field in fields]
    abstract_supertype = _extract_parametric_supertype(struct_node)

    promoted_constructor = quote
        @inline function $struct_name($(field_names...))
            promoted_type = $(GlobalRef(@__MODULE__, :_promoted_numeric_type))($values)
            return $struct_name{promoted_type}($(recasts...))
        end
    end
    concrete_converter = quote
        @inline function Base.convert(
                ::Type{<:$struct_name{TargetScalar}},
                value::$struct_name
        ) where {TargetScalar <: Real}
            return $struct_name{TargetScalar}($(target_recasts...))
        end
        @inline Base.eltype(::Type{<:$struct_name{Scalar}}) where {Scalar} = Scalar
        @inline Base.eltype(::$struct_name{Scalar}) where {Scalar} = Scalar
        @inline $(GlobalRef(@__MODULE__, :recast))(
            ::Type{TargetScalar},
            value::$struct_name
        ) where {TargetScalar <: Real} = $struct_name{TargetScalar}($(target_recasts...))
    end

    abstract_converter = abstract_supertype === nothing ? nothing :
                         quote
        @inline function Base.convert(
                ::Type{<:$abstract_supertype{TargetScalar}},
                value::$struct_name
        ) where {TargetScalar <: Real}
            return $struct_name{TargetScalar}($(target_recasts...))
        end
    end

    generated = abstract_converter === nothing ?
                (promoted_constructor, concrete_converter) :
                (promoted_constructor, concrete_converter, abstract_converter)
    rebuilt = _rebuild_ast(
        raw,
        struct_node,
        struct_node,
        generated
    )
    return esc(rebuilt)
end
