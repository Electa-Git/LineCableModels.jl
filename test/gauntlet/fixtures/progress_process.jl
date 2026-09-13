# Synthetic process fixture: no native solver or numerical campaign artifacts.
include(joinpath(@__DIR__, "..", "..", "..", "gauntlet", "Gauntlet.jl"))
using .Gauntlet, LineCableModels
const root=ARGS[2]
if ARGS[1]=="watch"
    Base.exit_on_sigint(false)
    write(joinpath(root,"watch.ready"),"ready")
    Gauntlet.watch_campaign(root;io=stdout,session=ARGS[3])
elseif ARGS[1]=="produce"
    function await_command(command)
        while !isfile(joinpath(root,command))
            sleep(0.05)
        end
    end
    function opaque(seconds)
        deadline=time_ns()+round(UInt64,seconds*1e9)
        value=0.0
        while time_ns()<deadline
            value=sin(value+0.1)
        end
        return value
    end
    Gauntlet._with_campaign_progress(root,[:a],"smoke";io=devnull) do tracker
        receiver=LineCableModels.progress_receiver()
        LineCableModels.report_progress(receiver,(kind=:benchmark,benchmark=:a,attempt="one",state=:running))
        tracker.rows[1]["case"]="long case "*repeat("decreasing text ",6)
        LineCableModels.report_progress(receiver,(kind=:operand,benchmark=:a,attempt="one",role=:reference,
            state=:running,backend="Synthetic",mode="ordinary"))
        await_command("short")
        tracker.rows[1]["case"]="short"
        Gauntlet._publish_progress!(tracker;force=true)
        await_command("opaque")
        LineCableModels.report_progress(receiver,(benchmark=:a,role=:reference,stage=:solving))
        Gauntlet._publish_progress!(tracker;force=true)
        opaque(4.0) # No Julia yield or producer publication during this interval.
        await_command("sample")
        Gauntlet._performance_span(;sample=1,samples=1,role=:reference) do
            opaque(4.0)
        end
        await_command("finish")
        LineCableModels.report_progress(receiver,(kind=:operand,benchmark=:a,attempt="one",role=:reference,state=:complete,seconds=8.0))
        LineCableModels.report_progress(receiver,(kind=:benchmark,benchmark=:a,attempt="one",state=:complete,verdict_failed=true))
    end
else
    error("unknown fixture mode")
end
