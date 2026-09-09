terminal_gateway_clock()=time_ns()/1e9

function claim_terminal_attachment!(service::TerminalCoordinator,principal::Principal,id::UUID)
    deadline=terminal_gateway_clock()+5
    while true
        # Recheck ownership/lease after each wait. A closing predecessor is not
        # an active writer, but its slot remains reserved until exact teardown.
        fence=remote_terminal_fence(service,principal,id)
        attachment=lock(service.lock) do
            service.closed && throw(AccessDenied(503,"Terminal channel is stopped"))
            prior=get(service.attachments,fence.lease_id,nothing)
            if prior!==nothing
                prior.closed || throw(AccessDenied(409,"Terminal already has a browser attachment"))
                return nothing
            end
            length(service.attachments)<min(256,service.coordinator.assignments.limits.total) || throw(CapacityUnavailable())
            value=TerminalAttachment(fence,uuid4(),current_task(),nothing,nothing,nothing,
                terminal_gateway_clock(),Inf,false)
            service.attachments[fence.lease_id]=value
            return value
        end
        attachment===nothing || return attachment
        terminal_gateway_clock()<deadline || throw(AccessDenied(409,"Previous terminal attachment cleanup is pending"))
        sleep(0.05)
    end
end

function authorize_terminal_attachment(service,attachment)
    !service.closed && !attachment.closed || throw(AccessDenied(409,"Terminal attachment is closed"))
    remote_terminal_fence(service,Principal(attachment.fence.owner),UUID(attachment.fence.lease_id))==attachment.fence ||
        throw(AccessDenied(409,"Terminal assignment changed"))
    lock(service.lock) do
        get(service.attachments,attachment.fence.lease_id,nothing)===attachment ||
            throw(AccessDenied(409,"Terminal attachment is no longer current"))
    end
    return nothing
end

function read_terminal_browser_frame(message)
    message isa String && ncodeunits(message)<=Protocol.MAX_TERMINAL_FRAME_BYTES ||
        throw(AccessDenied(400,"Expected a bounded terminal JSON frame"))
    try
        object=JSON3.read(message)
        object isa JSON3.Object || throw(ArgumentError("object required"))
        validate_gateway_json(object)
        return JSON3.read(message,Dict{String,Any})
    catch
        throw(AccessDenied(400,"Invalid terminal frame"))
    end
end

function terminal_browser_request(service,attachment,data)
    requested_fields(data,("action","request_id","session_id","input_sequence","after","columns","rows","bytes","retry"))
    data["action"] isa String && data["retry"] isa Bool || throw(AccessDenied(400,"Invalid terminal action"))
    action=data["action"]
    writer=action in ("read","status") ? nothing : attachment.writer
    writer === nothing && !(action in ("read","status")) && throw(AccessDenied(400,"Terminal writer is not attached"))
    # Decode through the same strict protocol grammar as the remote agent. The
    # browser cannot supply an owner, worker, fence, profile, revision or command.
    wire=merge(data,Dict("protocol_version"=>"2.0","fence"=>attachment.fence,"revision"=>1,"writer_id"=>writer))
    delete!(wire,"retry")
    command=Protocol.decode_terminal_message(Protocol.TerminalCommand,JSON3.write(wire))
    return request_terminal!(service,Principal(attachment.fence.owner),UUID(attachment.fence.lease_id),action;
        session_id=command.session_id,writer_id=command.writer_id,input_sequence=command.input_sequence,
        after=command.after,columns=command.columns,rows=command.rows,bytes=command.bytes,
        request_id=UUID(command.request_id),retry=data["retry"])
end

function send_terminal_browser!(service,attachment,payload)
    authorize_terminal_attachment(service,attachment)
    message=JSON3.write(payload)
    ncodeunits(message)<=Protocol.MAX_TERMINAL_FRAME_BYTES || throw(ArgumentError("terminal browser report exceeds its bound"))
    attachment.write_deadline=terminal_gateway_clock()+3
    try
        HTTP.WebSockets.send(attachment.socket,message)
    finally
        attachment.write_deadline=Inf
    end
    return nothing
end

function terminal_browser_report(report::Protocol.TerminalReport)
    (;kind="report",request_id=report.request_id,accepted=report.accepted,reason=report.reason,
        session_id=report.session_id,phase=report.phase,writer_connected=report.writer_connected,
        input_sequence=report.input_sequence,cursor=report.cursor,output_sequence=report.output_sequence,
        gap=report.gap,bytes=report.bytes,failure=report.failure,cleanup_pending=report.cleanup_pending)
end

function disconnect_terminal_attachment!(service,attachment)
    (attachment.writer === nothing || attachment.session === nothing || service.closed) && return nothing
    # Preserve uncertain action identity: cleanup does not supersede its pending
    # flight. Independent agent-side presence expiry covers loss of this message.
    try
        flight=request_terminal!(service,Principal(attachment.fence.owner),UUID(attachment.fence.lease_id),"disconnect";
            session_id=attachment.session,writer_id=attachment.writer)
        deadline=terminal_gateway_clock()+5
        while flight.report === nothing && terminal_gateway_clock()<deadline && !service.closed &&
                terminal_remote_clock(service)<flight.deadline
            sleep(0.02)
        end
    catch
        # Never expose raw input, broker records or remote exception text.
    end
    return nothing
end

function relay_terminal_socket!(service,attachment)
    socket=attachment.socket
    inbox=Channel{Dict{String,Any}}(1)
    receiver=@async try
        rate_start=terminal_gateway_clock();received=0
        while !attachment.closed && isopen(socket.readchannel)
            if !isready(socket);sleep(0.01);continue;end
            authorize_terminal_attachment(service,attachment)
            data=read_terminal_browser_frame(HTTP.WebSockets.receive(socket))
            now=terminal_gateway_clock()
            if now-rate_start>=1;rate_start=now;received=0;end
            received+=1;received<=128 || throw(AccessDenied(429,"Terminal request rate exceeded"))
            # One action may be in flight and one may wait. Additional pipelining
            # is a protocol error, not an unbounded keyboard/output queue.
            isready(inbox) && throw(AccessDenied(429,"Terminal requests are not serialized"))
            attachment.seen_at=now
            put!(inbox,data)
        end
    catch error
        @debug "Private terminal receiver rejected a frame" error_type=typeof(error) reason=(error isa AccessDenied ? error.reason : "internal or transport failure")
    finally
        @debug "Private terminal receiver ended" channel_open=isopen(socket.readchannel)
        isopen(inbox) && close(inbox)
        abort_terminal_attachment!(attachment)
    end
    watchdog=@async try
        while !attachment.closed
            authorize_terminal_attachment(service,attachment)
            now=terminal_gateway_clock()
            if now >= attachment.write_deadline
                @debug "Private terminal socket write deadline elapsed"
                break
            end
            if now-attachment.seen_at >= (attachment.writer === nothing ? 5 : 45)
                @debug "Private terminal socket idle deadline elapsed" attached=(attachment.writer !== nothing)
                break
            end
            sleep(0.05)
        end
    catch error
        @debug "Private terminal watchdog revoked attachment" error_type=typeof(error) reason=(error isa AccessDenied ? error.reason : "internal or transport failure")
    finally
        abort_terminal_attachment!(attachment)
    end
    try
        send_terminal_browser!(service,attachment,(kind="hello",schema_version=1,
            chunk_bytes=Protocol.MAX_TERMINAL_CHUNK_BYTES,serialized=true,keepalive_seconds=5))
        initial=take!(inbox)
        requested_fields(initial,("action","writer_id"))
        initial["action"]=="attach" && initial["writer_id"] isa String || throw(AccessDenied(400,"Expected private writer attachment"))
        Protocol.runtime_uuid(initial["writer_id"])
        attachment.writer=initial["writer_id"]
        send_terminal_browser!(service,attachment,(kind="attached",connection_id=string(attachment.id)))
        while !attachment.closed
            data=take!(inbox)
            authorize_terminal_attachment(service,attachment)
            flight=terminal_browser_request(service,attachment,data)
            deadline=terminal_gateway_clock()+(flight.action=="restart" ? 30 : 5)
            while flight.report === nothing && terminal_gateway_clock()<deadline &&
                    terminal_remote_clock(service)<flight.deadline && !attachment.closed
                sleep(0.02)
            end
            authorize_terminal_attachment(service,attachment)
            report=flight.report
            if report === nothing
                send_terminal_browser!(service,attachment,(kind="uncertain",request_id=flight.request_id,
                    action=flight.action,reason="reply_unavailable",automatic_retry=false))
            else
                report.accepted && report.session_id !== nothing && (attachment.session=report.session_id)
                send_terminal_browser!(service,attachment,terminal_browser_report(report))
                if flight.action=="disconnect" && report.accepted && !report.writer_connected
                    # The worker has confirmed release. Do not issue a second
                    # asynchronous disconnect from finally before freeing the slot.
                    attachment.session=nothing
                    break
                end
            end
        end
    catch error
        # Socket failures remain private and never fall back to an HTTP error
        # written onto an already-upgraded connection.
        @debug "Private terminal relay closed" error_type=typeof(error) reason=(error isa AccessDenied ? error.reason : "internal or transport failure") frames=stacktrace(catch_backtrace())[1:min(end,6)]
    finally
        abort_terminal_attachment!(attachment)
        isopen(inbox) && close(inbox)
        wait(receiver);wait(watchdog)
    end
    return nothing
end

"""Authorize an exact private terminal socket before upgrade and on each frame."""
function terminal_gateway_request(control,principal::Principal,stream,path::String)
    route=match(r"^/runtime/api/assignments/([a-f0-9-]{36})/terminal$",path)
    route === nothing && return false
    control === nothing && throw(AccessDenied(503,"Worker control is not configured"))
    service=control.terminals
    id=requested_uuid(route[1])
    remote_terminal_fence(service,principal,id)
    stream.message.method=="GET" && HTTP.WebSockets.isupgrade(stream.message) ||
        throw(AccessDenied(405,"Terminal requires an authorized WebSocket"))
    attachment=claim_terminal_attachment!(service,principal,id)
    try
        HTTP.WebSockets.upgrade(stream;check_origin=(_...)->true,
                maxframesize=Protocol.MAX_TERMINAL_FRAME_BYTES,maxfragmentation=1,compress=false) do socket
            attachment.socket=socket
            try
                bound_terminal_socket!(socket)
                relay_terminal_socket!(service,attachment)
            catch
                abort_terminal_attachment!(attachment)
            end
        end
    finally
        abort_terminal_attachment!(attachment)
        disconnect_terminal_attachment!(service,attachment)
        lock(service.lock) do
            get(service.attachments,attachment.fence.lease_id,nothing)===attachment &&
                delete!(service.attachments,attachment.fence.lease_id)
        end
    end
    return true
end
