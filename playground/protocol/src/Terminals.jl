"Maximum private terminal data chunk, in bytes."
const MAX_TERMINAL_CHUNK_BYTES=8192
"Maximum encoded private terminal frame, in bytes."
const MAX_TERMINAL_FRAME_BYTES=65536

"""
    TerminalCommand

Carry one transient, assignment-fenced terminal action. Revision orders all
actions within the lease; input_sequence separately orders admitted PTY input.
Unused fields must be empty/zero. The record contains neither executable launch
configuration nor a scientific job. Raw bytes must never enter durable logs.
"""
struct TerminalCommand <: RuntimeRecord
    "Assigned runtime protocol version."
    protocol_version::String
    "Explicit request UUID, preserved for an uncertain retry."
    request_id::String
    "Complete run, owner, worker incarnation and lease authority."
    fence::AssignmentFence
    "Strictly increasing command revision for this assignment."
    revision::Int
    "Open, status, read, input, resize, keepalive, disconnect, stop or restart."
    action::String
    "Exact stream UUID; absent only on initial open."
    session_id::Union{Nothing,String}
    "Sole writer UUID; absent on read-only status and output requests."
    writer_id::Union{Nothing,String}
    "Next input chunk sequence; zero for other actions."
    input_sequence::Int
    "Requested relative output cursor; zero except for read."
    after::Int
    "Terminal columns; zero unless opening, resizing or restarting."
    columns::Int
    "Terminal rows; zero unless opening, resizing or restarting."
    rows::Int
    "Bounded raw input bytes, including escape sequences and Ctrl-C."
    bytes::Vector{UInt8}
end

"""
    TerminalReport

Return one transient terminal acknowledgement or bounded output chunk. An input
acknowledgement confirms queue admission, never Julia evaluation. Startup and
ended states are explicit; neither cached reports nor output renew a live lease.
"""
struct TerminalReport <: RuntimeRecord
    "Assigned runtime protocol version."
    protocol_version::String
    "Exact command UUID."
    request_id::String
    "Exact command assignment, including worker incarnation."
    fence::AssignmentFence
    "Exact command revision."
    revision::Int
    "Whether the requested action was admitted."
    accepted::Bool
    "Fixed non-sensitive result token."
    reason::String
    "Current stream UUID, absent when open was not admitted."
    session_id::Union{Nothing,String}
    "Unknown, starting, ready, closing, exited or failed."
    phase::String
    "Whether the sole writer remains connected."
    writer_connected::Bool
    "Latest admitted input sequence, not an evaluation count."
    input_sequence::Int
    "Output cursor after the returned bytes."
    cursor::Int
    "Latest produced output cursor."
    output_sequence::Int
    "Whether requested earlier output has been evicted."
    gap::Bool
    "Bounded raw output; never ordinary diagnostics."
    bytes::Vector{UInt8}
    "Fixed session failure token, if present."
    failure::Union{Nothing,String}
    "Whether exact physical retirement remains unresolved."
    cleanup_pending::Bool
end

terminal_cursor(value::Integer)=0<=value<=9_007_199_254_740_991 ||
    throw(ArgumentError("invalid terminal cursor"))

function validate(command::TerminalCommand)
    command.protocol_version==RUNTIME_PROTOCOL_VERSION || throw(ArgumentError("unsupported terminal protocol"))
    runtime_uuid(command.request_id);validate(command.fence);runtime_sequence(command.revision)
    command.action in ("open","status","read","input","resize","keepalive","disconnect","stop","restart") ||
        throw(ArgumentError("unsupported terminal action"))
    if command.action=="open"
        command.session_id === nothing || throw(ArgumentError("initial open cannot choose a stream identity"))
    else
        command.session_id !== nothing || throw(ArgumentError("terminal action requires a stream identity"))
        runtime_uuid(command.session_id)
    end
    if command.action in ("read","status")
        command.writer_id === nothing || throw(ArgumentError("read-only terminal action has no writer"))
    else
        command.writer_id !== nothing || throw(ArgumentError("terminal action requires a writer"))
        runtime_uuid(command.writer_id)
    end
    if command.action in ("open","resize","restart")
        1<=command.columns<=1000 && 1<=command.rows<=1000 || throw(ArgumentError("invalid terminal cell dimensions"))
    else
        command.columns==0 && command.rows==0 || throw(ArgumentError("unused terminal dimensions must be zero"))
    end
    if command.action=="input"
        runtime_sequence(command.input_sequence)
        0<length(command.bytes)<=MAX_TERMINAL_CHUNK_BYTES || throw(ArgumentError("terminal input exceeds its bound"))
    else
        command.input_sequence==0 && isempty(command.bytes) || throw(ArgumentError("terminal action cannot carry input"))
    end
    command.action=="read" ? terminal_cursor(command.after) :
        command.after==0 || throw(ArgumentError("unused terminal cursor must be zero"))
    return command
end

function validate(report::TerminalReport)
    report.protocol_version==RUNTIME_PROTOCOL_VERSION || throw(ArgumentError("unsupported terminal protocol"))
    runtime_uuid(report.request_id);validate(report.fence);runtime_sequence(report.revision);runtime_token(report.reason)
    report.session_id === nothing || runtime_uuid(report.session_id)
    report.phase in ("unknown","starting","ready","closing","exited","failed") || throw(ArgumentError("invalid terminal phase"))
    foreach(terminal_cursor,(report.input_sequence,report.cursor,report.output_sequence))
    length(report.bytes)<=MAX_TERMINAL_CHUNK_BYTES && length(report.bytes)<=report.cursor<=report.output_sequence ||
        throw(ArgumentError("invalid terminal output bounds"))
    report.failure === nothing || runtime_token(report.failure)
    report.accepted || (isempty(report.bytes) && !report.gap) || throw(ArgumentError("rejected terminal report cannot carry output"))
    report.phase=="unknown" || report.session_id !== nothing || throw(ArgumentError("terminal state lacks stream identity"))
    report.session_id !== nothing || (!report.writer_connected && report.input_sequence==0 &&
        report.output_sequence==0 && !report.cleanup_pending) || throw(ArgumentError("missing stream cannot report live state"))
    return report
end

"""Decode a strict terminal record within its smaller private-frame bound."""
function decode_terminal_message(::Type{T},payload) where {T<:Union{TerminalCommand,TerminalReport}}
    bytes=payload isa AbstractString ? ncodeunits(payload) : length(payload)
    bytes<=MAX_TERMINAL_FRAME_BYTES || throw(ArgumentError("terminal frame exceeds 64 KiB"))
    return decode_runtime_message(T,payload)
end

"""Return an exact ephemeral channel; no terminal subject is a durable job subject."""
function terminal_subject(fence::AssignmentFence,direction::Symbol)
    validate(fence)
    direction in (:command,:report) || throw(ArgumentError("invalid terminal channel direction"))
    return "lcm.terminal.v2.$(fence.worker_id).$(fence.worker_boot).$(fence.lease_id).$(fence.generation).$direction"
end

Base.show(io::IO,record::Union{TerminalCommand,TerminalReport})=
    print(io,nameof(typeof(record)),"(",record.request_id,", <private>)")
