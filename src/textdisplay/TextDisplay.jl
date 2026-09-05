# Package-local formatting support. Domain owners select semantic names, fields,
# units, and child relationships; this module supplies bounded writers only.
module TextDisplay

import ..Units

"Return the explicitly declared semantic display name for an owned type."
function name end

function _fallback(value)
    return sprint(show, value; context = :compact => true)
end

"Format one numeric value with at most `sigdigits` significant digits."
function value(number::Real; sigdigits::Integer = 6)
    sigdigits > 0 || throw(ArgumentError("sigdigits must be positive"))
    isnan(number) && return "NaN"
    isinf(number) && return signbit(number) ? "−∞" : "∞"
    iszero(number) && return "0"
    rounded = try
        round(number; sigdigits)
    catch
        return _fallback(number)
    end
    if rounded isa AbstractFloat && isinteger(rounded) &&
            abs(rounded) <= typemax(Int)
        return string(round(Int, rounded))
    end
    return replace(string(rounded), "e+" => "e")
end

function value(number::Complex; sigdigits::Integer = 6)
    real_text = value(real(number); sigdigits)
    imag_text = value(abs(imag(number)); sigdigits)
    sign = signbit(imag(number)) ? " − " : " + "
    return string(real_text, sign, imag_text, "im")
end

value(item; sigdigits::Integer = 6) = _fallback(item)

const _ENGINEERING_PREFIX = (
    -24 => :yocto,
    -21 => :zepto,
    -18 => :atto,
    -15 => :femto,
    -12 => :pico,
    -9 => :nano,
    -6 => :micro,
    -3 => :milli,
    0 => :base,
    3 => :kilo,
    6 => :mega,
    9 => :giga,
    12 => :tera,
    15 => :peta,
    18 => :exa,
    21 => :zetta,
    24 => :yotta
)

function _engineering_exponent(number::Real)
    iszero(number) && return 0
    magnitude = try
        abs(Float64(number))
    catch
        return 0
    end
    isfinite(magnitude) || return 0
    exponent = 3 * floor(Int, log10(magnitude) / 3)
    return clamp(exponent, first(first(_ENGINEERING_PREFIX)), first(last(_ENGINEERING_PREFIX)))
end

_prefix(exponent::Integer) = only(last(pair) for pair in _ENGINEERING_PREFIX
    if first(pair) == exponent)

function _physical_unit(::Val{:meter}, prefix::Symbol)
    return Units.label(Units.Unit(:meter, prefix))
end
function _physical_unit(::Val{:ohm_meter}, prefix::Symbol)
    return string(Units.label(Units.Unit(:ohm, prefix)), "·m")
end
function _physical_unit(::Val{:hertz}, prefix::Symbol)
    return Units.label(Units.Unit(:hertz, prefix))
end
_physical_unit(::Val{:celsius}, ::Symbol) = "°C"
_physical_unit(::Val{:kelvin_inverse}, ::Symbol) = "K⁻¹"
_physical_unit(::Val{:degree}, ::Symbol) = "°"
_physical_unit(::Val{:dimensionless}, ::Symbol) = ""

_adaptive_prefix(::Val{:meter}) = true
_adaptive_prefix(::Val{:ohm_meter}) = true
_adaptive_prefix(::Val{:hertz}) = true
_adaptive_prefix(::Val{:celsius}) = false
_adaptive_prefix(::Val{:kelvin_inverse}) = false
_adaptive_prefix(::Val{:degree}) = false
_adaptive_prefix(::Val{:dimensionless}) = false

function _append_unit(number::AbstractString, unit::AbstractString)
    isempty(unit) && return String(number)
    unit == "°" && return string(number, unit)
    return string(number, " ", unit)
end

"Format a physical value using an adaptive engineering prefix."
function engineering(number::Real, unit::Symbol; sigdigits::Integer = 6)
    if !isfinite(number)
        return _append_unit(value(number; sigdigits), _physical_unit(Val(unit), :base))
    end
    exponent = _adaptive_prefix(Val(unit)) ? _engineering_exponent(number) : 0
    prefix = _prefix(exponent)
    scaled = number / (10.0^exponent)
    return _append_unit(
        value(scaled; sigdigits),
        _physical_unit(Val(unit), prefix)
    )
end

"Format an angle stored in radians."
function angle(radians::Real; sigdigits::Integer = 6)
    if isfinite(radians) && !iszero(radians)
        turns = radians / π
        nearest = round(turns)
        if abs(turns) >= 1 && isapprox(turns, nearest; rtol = 0, atol = 16eps(float(turns)))
            nearest == 1 && return "π"
            nearest == -1 && return "−π"
            return string(value(nearest; sigdigits), "π")
        end
    end
    return _append_unit(value(rad2deg(radians); sigdigits), "°")
end

"Format a result quantity using the units already owned by `Units`."
function quantity(
        observed,
        scientific_quantity::Units.Quantity,
        basis::Symbol = :pul;
        sigdigits::Integer = 6
)
    native = Units.native_unit(scientific_quantity, basis)
    displayed = Units.display_unit(scientific_quantity, basis)
    converted = observed * Units.scale_factor(native, displayed)
    unit = replace(Units.label(displayed), "." => "·")
    return _append_unit(value(converted; sigdigits), unit)
end

function _fit(text::AbstractString, width::Integer)
    width <= 0 && return ""
    textwidth(text) <= width && return String(text)
    width == 1 && return "…"
    characters = collect(text)
    while !isempty(characters) && textwidth(join(characters)) > width - 1
        pop!(characters)
    end
    return string(join(characters), "…")
end

function _display_width(io::IO)
    _, columns = displaysize(io)
    return max(columns, 20)
end

"Write an explicitly named flat field record in compact or aligned form."
function fields(
        io::IO,
        semantic_name::AbstractString,
        displayed::NamedTuple;
        multiline::Bool = false
)
    retained = Pair{Symbol, Any}[
        key => item for (key, item) in pairs(displayed) if item !== nothing
    ]
    if !multiline || get(io, :compact, false)
        print(io, semantic_name, "(")
        for (index, pair) in enumerate(retained)
            index > 1 && print(io, ", ")
            print(io, first(pair), "=", last(pair))
        end
        print(io, ")")
        return nothing
    end
    print(io, semantic_name)
    isempty(retained) && return nothing
    key_width = maximum(textwidth(String(first(pair))) for pair in retained)
    available = max(_display_width(io) - key_width - 4, 1)
    limited = get(io, :limit, true)
    rows, _ = displaysize(io)
    visible_count = limited ? min(length(retained), max(rows - 2, 0)) : length(retained)
    for pair in Iterators.take(retained, visible_count)
        key = String(first(pair))
        print(io, '\n', "  ", rpad(key, key_width), "  ", _fit(string(last(pair)), available))
    end
    omitted = length(retained) - visible_count
    omitted > 0 && print(io, '\n', "  ⋮ $omitted more fields")
    return nothing
end

function _tree_state(io::IO)
    limited = get(io, :limit, true)
    rows, _ = displaysize(io)
    return (
        # Reserve one row for an explicit truncation summary. This keeps the
        # complete output within `displaysize(io)` even when the row budget is
        # exhausted in a nested branch.
        remaining = Ref(limited ? max(rows - 2, 0) : typemax(Int)),
        depth = limited ? clamp(rows ÷ 4, 1, 4) : typemax(Int),
        children = limited ? 8 : typemax(Int),
        truncation = Ref(""),
    )
end

_tree_label(node::NamedTuple) = String(node.label)
_tree_label(node) = string(node)
_tree_children(node::NamedTuple) = get(node, :children, ())
_tree_children(node) = ()
_tree_noun(node::NamedTuple) = String(get(node, :noun, "items"))
_tree_noun(node) = "items"

function _write_tree_line(io::IO, prefix, connector, label)
    available = max(_display_width(io) - textwidth(prefix) - 3, 1)
    print(io, '\n', prefix, connector, _fit(label, available))
    return nothing
end

function _write_tree_children!(io::IO, nodes, prefix, level, state)
    entries = collect(nodes)
    isempty(entries) && return nothing

    visible_count = min(length(entries), state.children)
    omitted_by_count = length(entries) - visible_count
    for index in 1:visible_count
        if state.remaining[] <= 0
            isempty(state.truncation[]) &&
                (state.truncation[] =
                    "$(length(entries) - index + 1) more $(_tree_noun(first(entries)))")
            return nothing
        end
        node = entries[index]
        children = collect(_tree_children(node))
        has_following = index < visible_count || omitted_by_count > 0
        connector = has_following ? "├─ " : "└─ "
        _write_tree_line(io, prefix, connector, _tree_label(node))
        state.remaining[] -= 1

        isempty(children) && continue
        continuation = string(prefix, has_following ? "│  " : "   ")
        if level >= state.depth
            if state.remaining[] <= 0
                isempty(state.truncation[]) &&
                    (state.truncation[] = "$(length(children)) nested $(_tree_noun(node))")
                return nothing
            end
            noun = _tree_noun(node)
            _write_tree_line(
                io,
                continuation,
                "└─ ",
                "⋮ $(length(children)) nested $noun"
            )
            state.remaining[] -= 1
        else
            _write_tree_children!(io, children, continuation, level + 1, state)
        end
    end

    if omitted_by_count > 0
        message = "$omitted_by_count more $(_tree_noun(first(entries)))"
        if state.remaining[] > 0
            _write_tree_line(io, prefix, "└─ ", "⋮ $message")
            state.remaining[] -= 1
        else
            isempty(state.truncation[]) && (state.truncation[] = message)
        end
    end
    return nothing
end

"""
Write a deterministic, depth-, row-, and width-bounded tree.

Each child is either a string or a named tuple with `label`, optional
`children`, and optional `noun` fields. The tuples are detached display input;
the formatter does not inspect domain objects.
"""
function tree(
        io::IO,
        header::AbstractString,
        children;
        noun::AbstractString = "items"
)
    print(io, _fit(header, _display_width(io)))
    get(io, :compact, false) && return nothing
    entries = collect(children)
    isempty(entries) && return nothing
    normalized = map(entries) do entry
        entry isa NamedTuple ? entry : (label = string(entry), noun = noun)
    end
    state = _tree_state(io)
    _write_tree_children!(io, normalized, "", 1, state)
    isempty(state.truncation[]) || _write_tree_line(
        io,
        "",
        "└─ ",
        "⋮ $(state.truncation[])"
    )
    return nothing
end

"Generate the three Base display methods for one explicitly described flat record."
macro showfields(type_expression, semantic_name, mapping)
    mapping isa Expr && mapping.head === :-> || throw(ArgumentError(
        "@showfields expects `value -> (field = formatted, ...)`"
    ))
    binding = mapping.args[1]
    binding isa Symbol || throw(ArgumentError(
        "@showfields requires one untyped value binding"
    ))
    body = mapping.args[2]
    name_ref = GlobalRef(@__MODULE__, :name)
    fields_ref = GlobalRef(@__MODULE__, :fields)
    return quote
        $name_ref(::Type{<:$(esc(type_expression))}) = $(esc(semantic_name))
        function Base.summary(io::IO, value::$(esc(type_expression)))
            print(io, $name_ref(typeof(value)))
        end
        function Base.show(io::IO, value::$(esc(type_expression)))
            displayed = let $(esc(binding)) = value
                $(esc(body))
            end
            $fields_ref(io, $name_ref(typeof(value)), displayed; multiline = false)
        end
        function Base.show(
                io::IO,
                ::MIME"text/plain",
                value::$(esc(type_expression))
        )
            displayed = let $(esc(binding)) = value
                $(esc(body))
            end
            $fields_ref(io, $name_ref(typeof(value)), displayed; multiline = true)
        end
    end
end

public name, value, engineering, angle, quantity, fields, tree

Base.summary(io::IO, unit::Units.Unit) = print(io, "Physical unit ", Units.label(unit))
Base.show(io::IO, unit::Units.Unit) = print(io, Units.label(unit))
Base.show(io::IO, ::MIME"text/plain", unit::Units.Unit) = show(io, unit)

Base.summary(io::IO, unit::Units.UnitExpr) =
    print(io, "Unit expression ", Units.label(unit))
Base.show(io::IO, unit::Units.UnitExpr) = print(io, Units.label(unit))
Base.show(io::IO, ::MIME"text/plain", unit::Units.UnitExpr) = show(io, unit)

function Base.summary(io::IO, quantity::Units.Quantity)
    print(io, Units.label(quantity))
end
function Base.show(io::IO, quantity::Units.Quantity)
    scientific_symbol = Units.symbol(quantity)
    isempty(scientific_symbol) ? print(io, Units.label(quantity)) :
        print(io, scientific_symbol, " · ", Units.label(quantity))
end
Base.show(io::IO, ::MIME"text/plain", quantity::Units.Quantity) = show(io, quantity)

end
