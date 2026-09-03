const REPEATER_STYLES_PATH = normpath(joinpath(@__DIR__, "..", "assets", "repeater.css"))
include_dependency(REPEATER_STYLES_PATH)
const REPEATER_STYLES = read(REPEATER_STYLES_PATH, String)

"""
    RepeatEntry(fields, content; id=uuid4(), extract=_snapshot_value)

Represent one independently mounted row in a [`Repeater`](@ref).

`fields` is an opaque value bundle, normally a named tuple containing the
`.value` observables of the controls rendered by `content`. The component does
not prescribe field names, field count, or control types. `extract(fields)`
returns the row value used by [`snapshot`](@ref).

# Arguments

- `fields`: Arbitrary row state.
- `content`: Bonito-renderable row content.

# Keywords

- `id`: Stable row identity used while inserting, moving, or removing rows.
- `extract`: Function that converts `fields` to a plain value snapshot.

# Returns

- A keyed repeat entry.
"""
struct RepeatEntry{F,C,X}
    "Stable identity used by Bonito's keyed renderer."
    id::UUID

    "Arbitrary state bundle owned by the row."
    fields::F

    "Bonito-renderable controls or content for the row."
    content::C

    "Function that extracts a plain value from the state bundle."
    extract::X
end

function RepeatEntry(fields, content; id::UUID=uuid4(), extract=_snapshot_value)
    applicable(extract, fields) || throw(ArgumentError(
        "repeat-entry extract function does not accept the supplied fields"
    ))
    return RepeatEntry(id, fields, content, extract)
end

_snapshot_value(value::Observable) = _snapshot_value(value[])
_snapshot_value(value::NamedTuple{names}) where {names} =
    NamedTuple{names}(map(_snapshot_value, Tuple(value)))
_snapshot_value(value::Tuple) = map(_snapshot_value, value)
_snapshot_value(value::AbstractArray) = map(_snapshot_value, value)
_snapshot_value(value::AbstractDict) =
    Dict(key => _snapshot_value(item) for (key, item) in pairs(value))
_snapshot_value(value) = value

"""
    snapshot(entry::RepeatEntry)

Extract a plain value from `entry.fields` with the entry's extraction
function.

# Arguments

- `entry`: Repeat entry to inspect.

# Returns

- The extracted row value.
"""
snapshot(entry::RepeatEntry) = entry.extract(entry.fields)

"""
    Repeater(factory; initial=1, duplicate=nothing, min_entries=0,
             max_entries=nothing, label="Entries", add_label="Add entry",
             reorderable=true)

Manage an ordered collection of arbitrary [`RepeatEntry`](@ref) rows.

The zero-argument `factory` creates a fresh entry for the Add action. `initial`
may be a nonnegative count or a vector of preconstructed entries. When supplied,
`duplicate(value)` receives the plain [`snapshot`](@ref) of the source row and
must return a fresh entry.

# Arguments

- `factory`: Zero-argument function that constructs a fresh `RepeatEntry`.

# Keywords

- `initial`: Initial entry count or collection; defaults to `1`.
- `duplicate`: Optional function for cloning a row from its extracted value.
- `min_entries`: Smallest permitted row count; defaults to `0`.
- `max_entries`: Largest permitted row count, or `nothing` for no limit.
- `label`: Accessible collection label.
- `add_label`: Label for the Add action.
- `reorderable`: Whether move actions are rendered.

# Returns

- A reusable repeater component.

# Errors

- Throws `ArgumentError` for inconsistent limits, duplicate row identities,
  incompatible callbacks, or an invalid initial collection.
"""
struct Repeater{F,D}
    "Ordered repeat entries observed by the keyed renderer."
    entries::Observable{Vector{RepeatEntry}}

    "Zero-argument function that constructs a fresh entry."
    factory::F

    "Optional function that reconstructs an entry from an extracted value."
    duplicate_entry::D

    "Smallest permitted entry count."
    min_entries::Int

    "Largest permitted entry count, or `nothing` when unlimited."
    max_entries::Union{Nothing,Int}

    "Accessible label for the repeated collection."
    label::String

    "Visible label for the Add action."
    add_label::String

    "Whether the rendered component exposes row movement actions."
    reorderable::Bool
end

function _new_repeat_entry(factory)
    entry = factory()
    entry isa RepeatEntry || throw(ArgumentError(
        "repeater factory must return a RepeatEntry, got $(typeof(entry))"
    ))
    return entry
end

function _initial_repeat_entries(initial, factory)
    if initial isa Integer
        initial >= 0 || throw(ArgumentError("initial entry count must be nonnegative"))
        return RepeatEntry[_new_repeat_entry(factory) for _ in 1:initial]
    end
    initial isa AbstractVector || throw(ArgumentError(
        "initial must be a nonnegative entry count or a vector of RepeatEntry values"
    ))
    all(entry -> entry isa RepeatEntry, initial) || throw(ArgumentError(
        "every initial repeater item must be a RepeatEntry"
    ))
    return RepeatEntry[initial...]
end

function _validate_repeat_entries(entries)
    ids = getfield.(entries, :id)
    allunique(ids) || throw(ArgumentError("repeat-entry identities must be unique"))
    return entries
end

function Repeater(
        factory;
        initial=1,
        duplicate=nothing,
        min_entries::Integer=0,
        max_entries::Union{Nothing,Integer}=nothing,
        label::AbstractString="Entries",
        add_label::AbstractString="Add entry",
        reorderable::Bool=true
    )
    applicable(factory) || throw(ArgumentError(
        "repeater factory must be callable without arguments"
    ))
    min_entries >= 0 || throw(ArgumentError("min_entries must be nonnegative"))
    if !isnothing(max_entries)
        max_entries >= min_entries || throw(ArgumentError(
            "max_entries must be greater than or equal to min_entries"
        ))
    end
    entries = _validate_repeat_entries(_initial_repeat_entries(initial, factory))
    length(entries) >= min_entries || throw(ArgumentError(
        "initial collection contains fewer than min_entries entries"
    ))
    if !isnothing(max_entries) && length(entries) > max_entries
        throw(ArgumentError("initial collection contains more than max_entries entries"))
    end
    if !isnothing(duplicate) && !isempty(entries)
        applicable(duplicate, snapshot(first(entries))) || throw(ArgumentError(
            "repeater duplicate function must accept an extracted row value"
        ))
    end
    return Repeater(
        Observable(entries),
        factory,
        duplicate,
        Int(min_entries),
        isnothing(max_entries) ? nothing : Int(max_entries),
        string(label),
        string(add_label),
        reorderable
    )
end

Base.length(repeater::Repeater) = length(repeater.entries[])
Base.isempty(repeater::Repeater) = isempty(repeater.entries[])
Base.getindex(repeater::Repeater, index::Integer) = repeater.entries[][index]
Base.iterate(repeater::Repeater, state...) = iterate(repeater.entries[], state...)

function _entry_index(repeater::Repeater, id::UUID)
    index = findfirst(entry -> entry.id == id, repeater.entries[])
    isnothing(index) && throw(KeyError(id))
    return index
end

function _check_repeat_capacity(repeater::Repeater)
    if !isnothing(repeater.max_entries) && length(repeater) >= repeater.max_entries
        throw(ArgumentError("repeater already contains max_entries rows"))
    end
    return nothing
end

function _replace_entries!(repeater::Repeater, entries)
    _validate_repeat_entries(entries)
    repeater.entries[] = RepeatEntry[entries...]
    return repeater
end

"""
    add!(repeater::Repeater)

Append a freshly constructed row.

# Arguments

- `repeater`: Collection to extend.

# Returns

- The new `RepeatEntry`.

# Errors

- Throws `ArgumentError` when `max_entries` has been reached or the factory
  does not return a fresh valid entry.
"""
function add!(repeater::Repeater)
    _check_repeat_capacity(repeater)
    entry = _new_repeat_entry(repeater.factory)
    any(existing -> existing.id == entry.id, repeater.entries[]) && throw(ArgumentError(
        "repeater factory returned an entry identity already in use"
    ))
    _replace_entries!(repeater, [repeater.entries[]; entry])
    return entry
end

"""
    remove!(repeater::Repeater, id::UUID)

Remove the row identified by `id`.

# Arguments

- `repeater`: Collection to modify.
- `id`: Stable identity of the row to remove.

# Returns

- The removed `RepeatEntry`.

# Errors

- Throws `KeyError` when `id` is absent and `ArgumentError` when removing the
  row would violate `min_entries`.
"""
function remove!(repeater::Repeater, id::UUID)
    length(repeater) > repeater.min_entries || throw(ArgumentError(
        "repeater must retain at least $(repeater.min_entries) rows"
    ))
    index = _entry_index(repeater, id)
    entries = copy(repeater.entries[])
    removed = splice!(entries, index)
    _replace_entries!(repeater, entries)
    return removed
end

"""
    duplicate!(repeater::Repeater, id::UUID)

Insert a fresh row after `id`, initialized from the source row snapshot.

# Arguments

- `repeater`: Collection to modify.
- `id`: Stable identity of the source row.

# Returns

- The new `RepeatEntry`.

# Errors

- Throws `ArgumentError` when duplication is unsupported, `max_entries` has
  been reached, or the duplicate callback returns an invalid entry. Throws
  `KeyError` when `id` is absent.
"""
function duplicate!(repeater::Repeater, id::UUID)
    isnothing(repeater.duplicate_entry) && throw(ArgumentError(
        "repeater does not define a duplicate callback"
    ))
    _check_repeat_capacity(repeater)
    index = _entry_index(repeater, id)
    entry = repeater.duplicate_entry(snapshot(repeater.entries[][index]))
    entry isa RepeatEntry || throw(ArgumentError(
        "repeater duplicate callback must return a RepeatEntry, got $(typeof(entry))"
    ))
    any(existing -> existing.id == entry.id, repeater.entries[]) && throw(ArgumentError(
        "repeater duplicate callback returned an entry identity already in use"
    ))
    entries = copy(repeater.entries[])
    insert!(entries, index + 1, entry)
    _replace_entries!(repeater, entries)
    return entry
end

"""
    move!(repeater::Repeater, id::UUID, destination::Integer)

Move the row identified by `id` to the one-based `destination` index.

# Arguments

- `repeater`: Collection to reorder.
- `id`: Stable identity of the row to move.
- `destination`: One-based destination index.

# Returns

- The moved `RepeatEntry`.

# Errors

- Throws `KeyError` when `id` is absent and `BoundsError` when `destination`
  lies outside the collection.
"""
function move!(repeater::Repeater, id::UUID, destination::Integer)
    checkbounds(repeater.entries[], destination)
    source = _entry_index(repeater, id)
    entries = copy(repeater.entries[])
    entry = splice!(entries, source)
    insert!(entries, destination, entry)
    _replace_entries!(repeater, entries)
    return entry
end

"""
    snapshot(repeater::Repeater)

Extract the ordered plain values of all rows.

# Arguments

- `repeater`: Collection to inspect.

# Returns

- A vector containing each row snapshot in display order.
"""
snapshot(repeater::Repeater) = map(snapshot, repeater.entries[])

struct RenderedRepeatEntry
    repeater::Repeater
    entry::RepeatEntry
end

function repeater_icon(name::Symbol)
    children = if name == :add
        (SVG.path(; d="M12 5v14M5 12h14"),)
    elseif name == :up
        (SVG.path(; d="m6 15 6-6 6 6"),)
    elseif name == :down
        (SVG.path(; d="m6 9 6 6 6-6"),)
    elseif name == :duplicate
        (
            SVG.rect(; x="8", y="8", width="11", height="11", rx="1"),
            SVG.path(; d="M16 8V5a1 1 0 0 0-1-1H5a1 1 0 0 0-1 1v10a1 1 0 0 0 1 1h3"),
        )
    elseif name == :remove
        (SVG.path(; d="M5 12h14"),)
    else
        throw(ArgumentError("unknown repeater icon: $name"))
    end
    return SVG.svg(
        children...;
        viewBox="0 0 24 24",
        fill="none",
        stroke="currentColor",
        var"stroke-width"="1.8",
        var"stroke-linecap"="round",
        var"stroke-linejoin"="round",
        var"aria-hidden"="true",
        focusable="false"
    )
end

function repeater_action(callback, session, label, icon, enabled)
    trigger = Observable(false)
    disabled = map(session, enabled) do allowed
        !allowed
    end
    on(session, trigger) do _
        enabled[] && callback()
        return nothing
    end
    return DOM.button(
        repeater_icon(icon);
        type="button",
        class="lc-repeater-action",
        title=label,
        var"aria-label"=label,
        disabled,
        onclick=js"event => $(trigger).notify(true)"
    )
end

function Bonito.jsrender(session::Session, rendered::RenderedRepeatEntry)
    repeater = rendered.repeater
    entry = rendered.entry
    position = map(session, repeater.entries) do entries
        something(findfirst(candidate -> candidate.id == entry.id, entries), 0)
    end
    can_move_up = map(session, position) do index
        repeater.reorderable && index > 1
    end
    can_move_down = map(session, position, repeater.entries) do index, entries
        repeater.reorderable && 0 < index < length(entries)
    end
    can_duplicate = map(session, repeater.entries) do entries
        !isnothing(repeater.duplicate_entry) &&
            (isnothing(repeater.max_entries) || length(entries) < repeater.max_entries)
    end
    can_remove = map(session, repeater.entries) do entries
        length(entries) > repeater.min_entries
    end

    actions = Any[]
    if repeater.reorderable
        push!(actions, repeater_action(session, "Move row up", :up, can_move_up) do
            move!(repeater, entry.id, position[] - 1)
        end)
        push!(actions, repeater_action(session, "Move row down", :down, can_move_down) do
            move!(repeater, entry.id, position[] + 1)
        end)
    end
    if !isnothing(repeater.duplicate_entry)
        push!(actions, repeater_action(session, "Duplicate row", :duplicate, can_duplicate) do
            duplicate!(repeater, entry.id)
        end)
    end
    push!(actions, repeater_action(session, "Remove row", :remove, can_remove) do
        remove!(repeater, entry.id)
    end)

    node = DOM.div(
        DOM.div(
            DOM.span("ENTRY"; class="lc-repeater-entry-label"),
            DOM.span(; class="lc-repeater-entry-number"),
            DOM.div(actions...; class="lc-repeater-entry-actions");
            class="lc-repeater-entry-header"
        ),
        DOM.div(entry.content; class="lc-repeater-entry-content");
        class="lc-repeater-entry",
        var"data-repeat-entry"=string(entry.id)
    )
    return Bonito.jsrender(session, node)
end

function Bonito.jsrender(session::Session, repeater::Repeater)
    rendered = Observable(RenderedRepeatEntry[
        RenderedRepeatEntry(repeater, entry) for entry in repeater.entries[]
    ])
    on(session, repeater.entries) do entries
        rendered[] = RenderedRepeatEntry[
            RenderedRepeatEntry(repeater, entry) for entry in entries
        ]
        return nothing
    end
    keyed_entries = KeyedList(rendered; key=item -> item.entry.id)
    count = map(session, repeater.entries) do entries
        count = length(entries)
        return count == 1 ? "1 row" : "$count rows"
    end
    can_add = map(session, repeater.entries) do entries
        isnothing(repeater.max_entries) || length(entries) < repeater.max_entries
    end
    add_disabled = map(session, can_add) do allowed
        !allowed
    end
    add_trigger = Observable(false)
    on(session, add_trigger) do _
        can_add[] && add!(repeater)
        return nothing
    end

    node = DOM.div(
        DOM.style(REPEATER_STYLES),
        DOM.div(
            DOM.span(repeater.label; class="lc-repeater-label"),
            DOM.output(count; class="lc-repeater-count");
            class="lc-repeater-header"
        ),
        DOM.div(keyed_entries; class="lc-repeater-list"),
        DOM.div(
            DOM.button(
                repeater_icon(:add),
                DOM.span(repeater.add_label);
                type="button",
                class="lc-repeater-add",
                disabled=add_disabled,
                onclick=js"event => $(add_trigger).notify(true)"
            );
            class="lc-repeater-footer"
        );
        class="lc-repeater",
        role="group",
        var"aria-label"=repeater.label
    )
    return Bonito.jsrender(session, node)
end
