function _macro_block(value, name::Symbol)
    value isa Expr && value.head === :block || throw(ArgumentError(
        "@$name requires a begin/end block"
    ))
    parts = Any[]
    for item in value.args
        item isa LineNumberNode && continue
        if item isa Expr && item.head === :macrocall &&
                first(item.args) in (GlobalRef(Core, Symbol("@doc")), Symbol("@doc"))
            item = last(item.args)
        end
        item isa AbstractString || item === nothing || item === :nothing ||
            push!(parts, item)
    end
    return parts
end

_macro_call(name::Symbol, args...) = Expr(
    :call, GlobalRef(@__MODULE__, name), args...
)

function _macro_keywords(values, name::Symbol)
    keywords = Pair{Symbol, Any}[]
    for value in values
        value isa Expr && value.head === :(=) && value.args[1] isa Symbol ||
            throw(ArgumentError("@$name accepts keyword assignments"))
        push!(keywords, value.args[1] => value.args[2])
    end
    allunique(first.(keywords)) || throw(ArgumentError(
        "@$name does not accept duplicate keyword assignments"
    ))
    return keywords
end

function _call_with_keywords(name::Symbol, arguments, keywords)
    parameters = Expr(:parameters, [Expr(:kw, key, value) for (key, value) in keywords]...)
    return Expr(:call, GlobalRef(@__MODULE__, name), parameters, arguments...)
end

"Build a completed cable from the expressions in one ordered block."
macro cable(identifier, block)
    parts = _macro_block(block, :cable)
    isempty(parts) && throw(ArgumentError("@cable requires at least one part"))
    return esc(_macro_call(:build, GlobalRef(DataModel, :CableDesign), identifier, parts...))
end

"Coalesce the conductive descendants in one ordered block as a terminal."
macro terminal(name, block)
    parts = _macro_block(block, :terminal)
    isempty(parts) && throw(ArgumentError("@terminal requires at least one part"))
    return esc(_macro_call(:terminal, name, parts...))
end

"Assemble the independent members in one block."
macro assembly(block)
    members = _macro_block(block, :assembly)
    isempty(members) && throw(ArgumentError("@assembly requires at least one member"))
    return esc(_macro_call(:assembly, members...))
end

"Contain the members in one block as a duct."
macro duct(arguments...)
    isempty(arguments) && throw(ArgumentError("@duct requires keywords and a block"))
    block = last(arguments)
    members = _macro_block(block, :duct)
    keywords = _macro_keywords(arguments[1:(end - 1)], :duct)
    keys = first.(keywords)
    :shape in keys || throw(ArgumentError("@duct requires shape"))
    :fill in keys || throw(ArgumentError("@duct requires fill"))
    return esc(_call_with_keywords(:duct, members, keywords))
end

function _coordinate_tuple(value, name::Symbol)
    value isa Expr && value.head === :tuple && length(value.args) in (2, 3) ||
        throw(ArgumentError("@$name requires a two- or three-coordinate tuple"))
    return copy(value.args)
end

"Place a physical subject using coordinate-tuple notation."
macro at(subject, coordinates, assignments...)
    values = _coordinate_tuple(coordinates, :at)
    keywords = _macro_keywords(assignments, :at)
    keys = first.(keywords)
    tuple_has_angle = length(values) == 3
    tuple_has_angle && :φ in keys && throw(ArgumentError(
        "@at cannot receive φ both in the coordinate tuple and as a keyword"
    ))
    angle = tuple_has_angle ? pop!(values) : nothing
    angle === nothing || pushfirst!(keywords, :φ => angle)
    return esc(_call_with_keywords(:at, (subject, values...), keywords))
end

function _formation_macro(name::Symbol, design, assignments)
    keywords = _macro_keywords(assignments, name)
    reserved = name === :trefoil ? (:spacing, :center, :φ0, :connections) :
               (:spacing, :center, :connections)
    explicit_connections = findall(pair -> first(pair) === :connections, keywords)
    length(explicit_connections) <= 1 || throw(ArgumentError(
        "@$name accepts one connections assignment"
    ))
    shorthand = filter(pair -> !(first(pair) in reserved), keywords)
    isempty(explicit_connections) || isempty(shorthand) || throw(ArgumentError(
        "@$name cannot mix connections with connection shorthand"
    ))
    if isempty(explicit_connections) && !isempty(shorthand)
        named = Expr(:tuple, Expr(
            :parameters,
            [Expr(:kw, key, value) for (key, value) in shorthand]...
        ))
        keywords = filter(pair -> first(pair) in reserved, keywords)
        push!(keywords, :connections => named)
    end
    for index in eachindex(keywords)
        key, value = keywords[index]
        if key === :center && value isa Expr && value.head === :tuple
            coordinates = _coordinate_tuple(value, name)
            center = length(coordinates) == 2 ?
                     _macro_call(:at, coordinates...) :
                     _call_with_keywords(
                :at,
                coordinates[1:2],
                [:φ => coordinates[3]]
            )
            keywords[index] = key => center
        end
    end
    return _call_with_keywords(name, (design,), keywords)
end

macro trefoil(design, assignments...)
    return esc(_formation_macro(:trefoil, design, assignments))
end

macro hflat(design, assignments...)
    return esc(_formation_macro(:hflat, design, assignments))
end

macro vflat(design, assignments...)
    return esc(_formation_macro(:vflat, design, assignments))
end

"Insert the deferred `n = capacity()` policy into a repeated-member call."
macro distribute(call)
    call isa Expr && call.head === :call || throw(ArgumentError(
        "@distribute requires a constructor call"
    ))
    rewritten = copy(call)
    rewritten.args = copy(call.args)
    callee = first(rewritten.args)
    callee isa Symbol && callee in (:wires, :stranded, :rope, :tape) || throw(
        ArgumentError("@distribute supports wires, stranded, rope, and tape")
    )
    parameters = findfirst(
        argument -> argument isa Expr && argument.head === :parameters,
        rewritten.args
    )
    if parameters === nothing
        insert!(rewritten.args, 2, Expr(:parameters, Expr(
            :kw, :n, _macro_call(:capacity)
        )))
    else
        original_parameters = rewritten.args[parameters]
        rewritten.args[parameters] = copy(original_parameters)
        kwargs = copy(original_parameters.args)
        rewritten.args[parameters].args = kwargs
        any(keyword -> keyword isa Expr && keyword.head === :kw &&
                       keyword.args[1] === :n, kwargs) && throw(ArgumentError(
            "@distribute cannot overwrite an explicit n"
        ))
        push!(kwargs, Expr(:kw, :n, _macro_call(:capacity)))
    end
    return esc(rewritten)
end
