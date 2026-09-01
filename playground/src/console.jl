"Channels accepted by `ConsoleEntry` and their visual roles."
const CONSOLE_CHANNELS = (:input, :stdout, :info, :warning, :error, :broker)

"""
    ConsoleEntry

Represent one message displayed by a [`ConsoleView`](@ref).

# Fields

- `channel`: Visual message channel. Supported values are `:input`, `:stdout`,
  `:info`, `:warning`, `:error`, and `:broker`.
- `text`: Message text. Newlines are split into separately aligned entries by
  [`append_console!`](@ref).
"""
struct ConsoleEntry
    "Visual message channel."
    channel::Symbol

    "Message text."
    text::String

    function ConsoleEntry(channel::Symbol, text::AbstractString)
        channel in CONSOLE_CHANNELS || throw(ArgumentError(
            "unsupported console channel: $channel"
        ))
        return new(channel, string(text))
    end
end

ConsoleEntry(text::AbstractString; channel::Symbol=:stdout) =
    ConsoleEntry(channel, text)

"""
    ConsoleView(; entries=ConsoleEntry[], title="Julia output",
                status="ready", max_entries=500)

Create a bounded, reactive, read-only console viewport. The component displays
messages supplied by Julia, a worker, or a broker; it does not evaluate code.

# Keywords

- `entries`: Initial console entries.
- `title`: Header shown above the message viewport.
- `status`: Initial status shown in the console header.
- `max_entries`: Maximum retained entry count. Older entries are discarded.

# Returns

- A `ConsoleView` whose `entries` and `status` observables can be updated
  without replacing the surrounding DOM.
"""
struct ConsoleView
    "Messages currently retained by the viewport."
    entries::Observable{Vector{ConsoleEntry}}

    "Header title."
    title::String

    "Reactive header status."
    status::Observable{String}

    "Monotonic token used to trigger browser-side viewport updates."
    revision::Observable{Int}

    "Maximum number of retained messages."
    max_entries::Int

    "Lock protecting task-driven message updates."
    update_lock::ReentrantLock
end

function ConsoleView(;
        entries=ConsoleEntry[],
        title::AbstractString="Julia output",
        status::AbstractString="ready",
        max_entries::Integer=500
    )
    max_entries > 0 || throw(ArgumentError("max_entries must be positive"))
    initial = ConsoleEntry[entry for entry in entries]
    length(initial) > max_entries &&
        (initial = initial[(end - Int(max_entries) + 1):end])
    return ConsoleView(
        Observable(initial),
        string(title),
        Observable(string(status)),
        Observable(0),
        Int(max_entries),
        ReentrantLock()
    )
end

function console_message_lines(message::AbstractString)
    normalized = replace(string(message), "\r\n" => "\n", '\r' => '\n')
    lines = split(normalized, '\n'; keepempty=true)
    endswith(normalized, '\n') && length(lines) > 1 && pop!(lines)
    return isempty(lines) ? [""] : lines
end

"""
    append_console!(console, message; channel=:stdout)

Append a message to `console`, splitting multiline messages into aligned rows.
Updates are serialized so worker or broker tasks can safely share one view.

# Arguments

- `console`: Destination console viewport.
- `message`: Text to append.

# Keywords

- `channel`: Visual channel assigned to each appended row.

# Returns

- The updated `ConsoleView`.
"""
function append_console!(
        console::ConsoleView,
        message::AbstractString;
        channel::Symbol=:stdout
    )
    channel in CONSOLE_CHANNELS || throw(ArgumentError(
        "unsupported console channel: $channel"
    ))
    additions = ConsoleEntry[
        ConsoleEntry(channel, line) for line in console_message_lines(message)
    ]
    lock(console.update_lock) do
        updated = vcat(console.entries[], additions)
        length(updated) > console.max_entries &&
            (updated = updated[(end - console.max_entries + 1):end])
        console.entries[] = updated
        console.revision[] += 1
    end
    return console
end

"""
    clear_console!(console)

Remove all messages from a console viewport and return the updated view.
"""
function clear_console!(console::ConsoleView)
    lock(console.update_lock) do
        console.entries[] = ConsoleEntry[]
        console.revision[] += 1
    end
    return console
end

"""
    set_console_status!(console, status)

Set the console header status and return the updated view.
"""
function set_console_status!(console::ConsoleView, status::AbstractString)
    console.status[] = string(status)
    return console
end

const CONSOLE_PREFIXES = Dict(
    :input => "julia>",
    :stdout => "",
    :info => "[ Info]",
    :warning => "[ Warn]",
    :error => "[Error]",
    :broker => "[ NATS]",
)

function console_entry_dom(entry::ConsoleEntry)
    return DOM.div(
        DOM.span(CONSOLE_PREFIXES[entry.channel]; class="lc-console-prefix"),
        DOM.span(entry.text; class="lc-console-message");
        class="lc-console-line lc-console-$(entry.channel)"
    )
end

function Bonito.jsrender(session::Session, console::ConsoleView)
    rendered_entries = map(session, console.entries) do entries
        if isempty(entries)
            return DOM.div(
                "No messages";
                class="lc-console-empty"
            )
        end
        return DOM.div(
            (console_entry_dom(entry) for entry in entries)...;
            class="lc-console-lines"
        )
    end
    attributes = Dict{Symbol,Any}(
        Symbol("aria-live") => "polite",
        Symbol("aria-atomic") => "false",
        Symbol("aria-label") => console.title,
    )
    viewport = DOM.div(
        rendered_entries;
        attributes...,
        class="lc-console-viewport",
        role="log",
        tabindex="0"
    )
    follow_output = js"""
    (() => {
        const viewport = $(viewport);
        const revision = $(console.revision);
        const follow = () => requestAnimationFrame(() => {
            viewport.scrollTop = viewport.scrollHeight;
        });
        revision.on(follow);
        follow();
    })();
    """
    return Bonito.jsrender(
        session,
        DOM.section(
            DOM.header(
                DOM.span(console.title; class="lc-console-title"),
                DOM.span(console.status; class="lc-console-status");
                class="lc-console-header"
            ),
            viewport,
            DOM.script(follow_output);
            class="lc-console"
        )
    )
end
