# Runtime observation never belongs to a numerical declaration or checkpoint.
const _EXECUTION_OBSERVATION = Base.ScopedValues.ScopedValue{Any}(
    (receiver = nothing, scope = (;), quiet = false))

"""
$(TYPEDSIGNATURES)

Return the current runtime receiver and scope, or `nothing` when observation is
disabled. Owners may retain this value at an outer loop boundary.
"""
function progress_receiver()
    context = _EXECUTION_OBSERVATION[]
    return context.quiet || context.receiver === nothing ? nothing :
           (context.receiver, context.scope)
end

"""
$(TYPEDSIGNATURES)

Publish absolute execution counters to a runtime receiver. `nothing` is a no-op.
Reports contain no numerical results and do not change the computation.
"""
report_progress(::Nothing, event::NamedTuple) = nothing
function report_progress(receiver::Tuple, event::NamedTuple)
    first(receiver)(merge(last(receiver), event))
end

"""
$(TYPEDSIGNATURES)

Evaluate `f()` with a task-scoped execution receiver and identifying `scope`.
Restore the previous context on return or exception. The receiver is not retained
in problems, formulations or results.
"""
function with_progress(f, receiver; scope::NamedTuple = (;))
    previous = _EXECUTION_OBSERVATION[]
    context = (receiver, scope = merge(previous.scope, scope), quiet = previous.quiet)
    return Base.ScopedValues.with(f, _EXECUTION_OBSERVATION => context)
end

"""
$(TYPEDSIGNATURES)

Identify a nested execution without replacing its receiver. With observation
disabled, evaluate `f()` directly.
"""
function with_progress_scope(f; scope...)
    previous = _EXECUTION_OBSERVATION[]
    previous.receiver === nothing && return f()
    return with_progress(f, previous.receiver; scope = (; scope...))
end

"""Return whether the current task is inside a quiet performance sample."""
performance_sample_active() = _EXECUTION_OBSERVATION[].quiet

"""
$(TYPEDSIGNATURES)

Evaluate `f()` without progress publication or ordinary Julia diagnostics.
Owned console loggers and native launchers also honor this task-scoped flag.
The caller must pause its renderer before entering and time only the desired
computation inside `f`. Exceptions and interrupts propagate normally.
"""
function with_performance_sample(f)
    context = (receiver = nothing, scope = (;), quiet = true)
    return Base.ScopedValues.with(_EXECUTION_OBSERVATION => context) do
        Logging.with_logger(f, Logging.NullLogger())
    end
end

"""
$(TYPEDSIGNATURES)

Print diagnostics without overwriting a live progress display. Receivers may
suspend rendering between the paired output notifications. No numerical work or
message text is sent to the receiver.
"""
function with_progress_output(f)
    receiver=progress_receiver()
    receiver === nothing && return f()
    report_progress(receiver, (kind = :output_begin,))
    try
        return f()
    finally
        report_progress(receiver, (kind = :output_end,))
    end
end

public progress_receiver, report_progress, with_progress, with_progress_scope,
       performance_sample_active, with_performance_sample, with_progress_output
