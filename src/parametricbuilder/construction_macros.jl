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

"""
    @cable identifier [nominal_data=value] [combine=:product] begin
        parts...
    end

Build a completed cable from physical declarations ordered from the centre
outward.

# Arguments

- `identifier`: Stable cable identifier.
- `parts`: Ordered physical cable declarations.

# Keywords

- `nominal_data=nothing`: Descriptive catalogue data stored with the design.
- `combine=:product`: Gridspace composition rule.

# Returns

- A completed `CableDesign`, or a `Gridspace{CableDesign}` when a direct input
  varies.
"""
macro cable(identifier, arguments...)
    isempty(arguments) && throw(ArgumentError(
        "@cable requires a begin/end block"
    ))
    block = last(arguments)
    parts = _macro_block(block, :cable)
    isempty(parts) && throw(ArgumentError("@cable requires at least one part"))
    keywords = _macro_keywords(arguments[1:(end - 1)], :cable)
    allowed = (:nominal_data, :combine)
    unsupported = filter(pair -> !(first(pair) in allowed), keywords)
    isempty(unsupported) || throw(ArgumentError(
        "@cable accepts only nominal_data and combine"
    ))
    call = if isempty(keywords)
        _macro_call(
            :build,
            GlobalRef(DataModel, :CableDesign),
            identifier,
            parts...
        )
    else
        _call_with_keywords(
            :build,
            (GlobalRef(DataModel, :CableDesign), identifier, parts...),
            keywords
        )
    end
    return esc(call)
end

"""
    @system identifier [keyword=value ...] begin
        placements...
    end

Build a completed line-cable system from placed cable declarations.

# Arguments

- `identifier`: Stable system identifier.
- `placements`: Ordered expressions that produce placed cable declarations.

# Keywords

- `environment=nothing`: Optional physical environment declaration.
- `line_length=1`: Physical line length \\[m\\].
- `combine=:product`: Gridspace composition rule.

# Returns

- A completed `LineCableSystem`, or a `Gridspace{LineCableSystem}` when a
  placement or keyword value varies.
"""
macro system(identifier, arguments...)
    isempty(arguments) && throw(ArgumentError(
        "@system requires a begin/end block"
    ))
    block = last(arguments)
    placements = _macro_block(block, :system)
    isempty(placements) && throw(ArgumentError(
        "@system requires at least one placed cable"
    ))
    keywords = _macro_keywords(arguments[1:(end - 1)], :system)
    allowed = (:environment, :line_length, :combine)
    unsupported = filter(pair -> !(first(pair) in allowed), keywords)
    isempty(unsupported) || throw(ArgumentError(
        "@system accepts only environment, line_length, and combine"
    ))
    pushfirst!(keywords, :system_id => identifier)
    declarations = Expr(:tuple, placements...)
    return esc(_call_with_keywords(
        :build,
        (GlobalRef(DataModel, :LineCableSystem), declarations),
        keywords
    ))
end

"""
    @earth [keyword=value ...] begin
        layers...
    end

Build a completed earth model from ordered `EarthLayer` declarations. For
horizontal interfaces, the block runs from the surface downward.
Semi-infinite air is implicit.

# Arguments

- `layers`: Ordered expressions that produce completed earth layers.

# Keywords

- `vertical_layers=false`: Whether earth interfaces are vertical.
- `air_layer=nothing`: Optional explicit semi-infinite air layer.
- `combine=:product`: Gridspace composition rule.

# Returns

- An `EarthModel`, or a `Gridspace{EarthModel}` when a layer or keyword value
  varies.
"""
macro earth(arguments...)
    isempty(arguments) && throw(ArgumentError(
        "@earth requires a begin/end block"
    ))
    block = last(arguments)
    layers = _macro_block(block, :earth)
    isempty(layers) && throw(ArgumentError(
        "@earth requires at least one earth layer"
    ))
    keywords = _macro_keywords(arguments[1:(end - 1)], :earth)
    allowed = (:vertical_layers, :air_layer, :combine)
    unsupported = filter(pair -> !(first(pair) in allowed), keywords)
    isempty(unsupported) || throw(ArgumentError(
        "@earth accepts only vertical_layers, air_layer, and combine"
    ))
    declarations = Expr(:tuple, layers...)
    return esc(_call_with_keywords(
        :build,
        (GlobalRef(EarthProps, :EarthModel), declarations),
        keywords
    ))
end

"""
    @terminal name [combine=:product] begin
        parts...
    end

Coalesce the conductive descendants in one ordered block as a terminal.

# Arguments

- `name`: Retained electrical terminal name.
- `parts`: Ordered physical cable declarations.

# Keywords

- `combine=:product`: Gridspace composition rule.

# Returns

- A terminal-owning `Group`, or a `Gridspace{Group}` when a direct input
  varies.
"""
macro terminal(name, arguments...)
    isempty(arguments) && throw(ArgumentError(
        "@terminal requires a begin/end block"
    ))
    block = last(arguments)
    parts = _macro_block(block, :terminal)
    isempty(parts) && throw(ArgumentError("@terminal requires at least one part"))
    keywords = _macro_keywords(arguments[1:(end - 1)], :terminal)
    all(pair -> first(pair) === :combine, keywords) || throw(ArgumentError(
        "@terminal accepts only combine"
    ))
    call = if isempty(keywords)
        _macro_call(:terminal, name, parts...)
    else
        _call_with_keywords(:terminal, (name, parts...), keywords)
    end
    return esc(call)
end

"""
    @assembly [combine=:product] begin
        members...
    end

Assemble independent physical members while preserving their terminal
identities.

# Arguments

- `members`: Independent physical declarations, optionally placed with `@at`.

# Keywords

- `combine=:product`: Gridspace composition rule.

# Returns

- An `Assembly`, or a `Gridspace{Assembly}` when a direct input varies.
"""
macro assembly(arguments...)
    isempty(arguments) && throw(ArgumentError(
        "@assembly requires a begin/end block"
    ))
    block = last(arguments)
    members = _macro_block(block, :assembly)
    isempty(members) && throw(ArgumentError("@assembly requires at least one member"))
    keywords = _macro_keywords(arguments[1:(end - 1)], :assembly)
    all(pair -> first(pair) === :combine, keywords) || throw(ArgumentError(
        "@assembly accepts only combine"
    ))
    call = if isempty(keywords)
        _macro_call(:assembly, members...)
    else
        _call_with_keywords(:assembly, members, keywords)
    end
    return esc(call)
end

function _enclosure_macro(name::Symbol, arguments)
    isempty(arguments) && throw(ArgumentError(
        "@$name requires keywords and a begin/end block"
    ))
    block = last(arguments)
    members = _macro_block(block, name)
    isempty(members) && throw(ArgumentError("@$name requires enclosed content"))
    keywords = _macro_keywords(arguments[1:(end - 1)], name)
    keys = first.(keywords)
    :shape in keys || throw(ArgumentError("@$name requires shape"))
    :fill in keys || throw(ArgumentError("@$name requires fill"))
    allowed = name === :pipe ? (:shape, :fill, :wall, :at, :combine) :
              (:shape, :fill, :wall, :formation, :at, :combine)
    unsupported = filter(pair -> !(first(pair) in allowed), keywords)
    isempty(unsupported) || throw(ArgumentError(
        "@$name accepts only $(join(allowed, ", "))"
    ))
    return _call_with_keywords(name, members, keywords)
end

"""
    @pipe shape=value fill=value [keyword=value ...] begin
        members...
    end

Contain one or more physical members inside a pipe cross-section.

# Arguments

- `members`: Enclosed physical declarations.

# Keywords

- `shape`: Intrinsic containing primitive.
- `fill`: Filling material or explicit filling region.
- `wall=nothing`: Optional outward wall declaration.
- `at=nothing`: Pipe pose relative to its parent frame.
- `combine=:product`: Gridspace composition rule.

# Returns

- An `Enclosure`, or a `Gridspace{Enclosure}` when a direct input varies.
"""
macro pipe(arguments...)
    return esc(_enclosure_macro(:pipe, arguments))
end

"""
    @duct shape=value fill=value [keyword=value ...] begin
        members...
    end

Contain one or more physical members inside a duct cross-section.

# Arguments

- `members`: Enclosed physical declarations.

# Keywords

- `shape`: Intrinsic containing primitive.
- `fill`: Filling material or explicit filling region.
- `wall=nothing`: Optional outward wall declaration.
- `formation=nothing`: Placement pattern for one repeated prototype.
- `at=nothing`: Duct pose relative to its parent frame.
- `combine=:product`: Gridspace composition rule.

# Returns

- An `Enclosure`, or a `Gridspace{Enclosure}` when a direct input varies.
"""
macro duct(arguments...)
    return esc(_enclosure_macro(:duct, arguments))
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
    reserved = (:φ, :connections, :combine)
    explicit_connections = findall(pair -> first(pair) === :connections, keywords)
    shorthand = filter(pair -> !(first(pair) in reserved), keywords)
    isempty(explicit_connections) || isempty(shorthand) || throw(ArgumentError(
        "@at cannot mix connections with connection shorthand"
    ))
    if isempty(explicit_connections) && !isempty(shorthand)
        named = Expr(:tuple, Expr(
            :parameters,
            [Expr(:kw, key, value) for (key, value) in shorthand]...
        ))
        keywords = filter(pair -> first(pair) in reserved, keywords)
        push!(keywords, :connections => named)
    end
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
    reserved = name === :trefoil ?
               (:spacing, :center, :φ0, :connections, :combine) :
               (:spacing, :center, :connections, :combine)
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

"""
    @trefoil design spacing=value [keyword=value ...]

Place three copies of `design` in a trefoil formation. `center`, `φ0`,
`connections`, and `combine` are forwarded to [`trefoil`](@ref); every other
keyword is terminal-connection shorthand.
"""
macro trefoil(design, assignments...)
    return esc(_formation_macro(:trefoil, design, assignments))
end

"""
    @hflat design spacing=value [keyword=value ...]

Place three copies of `design` in a horizontal flat formation. `center`,
`connections`, and `combine` are forwarded to [`hflat`](@ref); every other
keyword is terminal-connection shorthand.
"""
macro hflat(design, assignments...)
    return esc(_formation_macro(:hflat, design, assignments))
end

"""
    @vflat design spacing=value [keyword=value ...]

Place three copies of `design` in a vertical flat formation. `center`,
`connections`, and `combine` are forwarded to [`vflat`](@ref); every other
keyword is terminal-connection shorthand.
"""
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
