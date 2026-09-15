# Current bounded producer/watcher protocol; no numerical payloads.
include(joinpath(@__DIR__,"../../../gauntlet/Gauntlet.jl"))
using .Gauntlet,LineCableModels
length(ARGS)>=2 || error("expected mode and output directory")
root=ARGS[2]
if ARGS[1]=="watch"
    length(ARGS)==3 || error("expected session identity")
    Base.exit_on_sigint(false)
    write(joinpath(root,"watch.ready"),"current watcher ready")
    Gauntlet.watch_campaign(root;io=stdout,session=ARGS[3])
elseif ARGS[1]=="produce"
    function handshake(name)
        timedwait(()->isfile(joinpath(root,name)),60.0;pollint=.02)===:ok ||
            error("progress handshake timed out: $name")
    end
    # A bounded non-yielding interval checks observation by the separate watcher.
    function occupied(seconds)
        stop=time_ns()+round(UInt64,seconds*1e9)
        accumulator=0.0
        while time_ns()<stop
            accumulator=cos(accumulator+.25)
        end
        accumulator
    end
    Gauntlet._with_campaign_progress(root,[:current],"current-session";io=devnull) do tracker
        receiver=LineCableModels.progress_receiver()
        emit(event)=LineCableModels.report_progress(receiver,merge((benchmark=:current,attempt="fresh"),event))
        emit((kind=:benchmark,state=:running))
        tracker.rows[1]["case"]="expanded current label "*repeat("wide label ",8)
        emit((kind=:operand,role=:reference,state=:running,backend="ProcessControl",mode="ordinary"))
        handshake("compact")
        tracker.rows[1]["case"]="compact"
        Gauntlet._publish_progress!(tracker;force=true)
        handshake("busy")
        emit((role=:reference,stage=:solving))
        Gauntlet._publish_progress!(tracker;force=true)
        occupied(4.0)
        handshake("sample")
        Gauntlet._performance_span(;sample=1,samples=1,role=:reference) do
            occupied(4.0)
        end
        handshake("finish")
        emit((kind=:operand,role=:reference,state=:complete,seconds=8.0))
        emit((kind=:benchmark,state=:complete))
    end
else
    error("unknown current process mode")
end
