"Number of complete private terminal frames retained by one browser socket."
const TERMINAL_SOCKET_QUEUE=4

"""
    bound_terminal_socket!(socket)

Replace HTTP 2.6.6's unbounded receive channel before its sticky reader starts.
This is a deliberately isolated compatibility adapter: a different HTTP version
or an already-started reader fails closed. Per-frame and fragmentation limits
remain enforced by HTTP's codec. No global HTTP method or package is modified.
"""
function bound_terminal_socket!(socket::HTTP.WebSockets.WebSocket)
    Base.pkgversion(HTTP)==v"2.6.6" || throw(ArgumentError("terminal WebSocket adapter requires its verified HTTP version"))
    socket.readtask !== nothing && !istaskstarted(socket.readtask) ||
        throw(ArgumentError("terminal socket must be bounded before its reader starts"))
    socket.maxframesize<=Protocol.MAX_TERMINAL_FRAME_BYTES && socket.maxfragmentation==1 ||
        throw(ArgumentError("terminal socket limits were not installed"))
    original=socket.readchannel
    pending=Union{String,Vector{UInt8}}[]
    while isready(original)
        length(pending)<TERMINAL_SOCKET_QUEUE || throw(ArgumentError("too many early terminal frames"))
        push!(pending,take!(original))
    end
    replacement=Channel{Union{String,Vector{UInt8}}}(TERMINAL_SOCKET_QUEUE)
    socket.readchannel=replacement
    foreach(message->put!(replacement,message),pending)
    close(original)
    return nothing
end

"""Release a bounded socket even when a peer stops reading or floods its queue."""
function abort_terminal_socket!(socket::HTTP.WebSockets.WebSocket)
    # HTTP's reader can hold its send lock while queued put! waits. Release that
    # wait and the original transport before asking HTTP to finish its close.
    isopen(socket.readchannel) && close(socket.readchannel)
    try socket.close_transport!() catch end
    try close(socket) catch end
    return nothing
end

"""Own one bounded browser attachment, without owning the remote REPL process."""
mutable struct TerminalAttachment
    "Exact private assignment authority."
    fence::AssignmentFence
    "Fresh connection identity; distinct from the resumable writer UUID."
    id::UUID
    "Original HTTP handler task, joined by coordinator shutdown."
    task::Task
    "Original socket, absent until the authorized upgrade completes."
    socket::Union{Nothing,HTTP.WebSockets.WebSocket}
    "Writer UUID received in the first private frame only."
    writer::Union{Nothing,String}
    "Latest acknowledged stream UUID for best-effort disconnect."
    session::Union{Nothing,String}
    "Latest valid browser frame time in local monotonic seconds."
    seen_at::Float64
    "Current finite socket-send deadline, or infinity between sends."
    write_deadline::Float64
    "Permanent attachment closure; does not release the remote assignment."
    closed::Bool
end
Base.show(io::IO,::TerminalAttachment)=print(io,"TerminalAttachment(<private>)")

function abort_terminal_attachment!(attachment::TerminalAttachment)
    attachment.closed=true
    attachment.socket === nothing || abort_terminal_socket!(attachment.socket)
    return nothing
end
