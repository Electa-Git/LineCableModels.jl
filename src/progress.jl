# Runtime observation never belongs to a numerical declaration or checkpoint.
const _EXECUTION_OBSERVATION = Base.ScopedValues.ScopedValue{Any}(
    (receiver = nothing, scope = (;), quiet = false, sequence = nothing))

"""
$(TYPEDSIGNATURES)

Return the current runtime receiver and scope, or `nothing` when observation is
disabled. Owners may retain this value at an outer loop boundary.
"""
function progress_receiver()
    context = _EXECUTION_OBSERVATION[]
    return context.quiet || context.receiver === nothing ? nothing :
           (context.receiver, context.scope, context.sequence)
end

"""
$(TYPEDSIGNATURES)

Publish absolute execution counters to a runtime receiver. `nothing` is a no-op.
Reports contain no numerical results and do not change the computation.
"""
report_progress(::Nothing, event::NamedTuple) = nothing
function report_progress(receiver::Tuple, event::NamedTuple)
    sequence = Threads.atomic_add!(receiver[3], 1) + 1
    first(receiver)(merge(receiver[2], event, (; sequence)))
end

"""
$(TYPEDSIGNATURES)

Evaluate `f()` with a task-scoped execution receiver and identifying `scope`.
Restore the previous context on return or exception. The receiver is not retained
in problems, formulations or results.
"""
function with_progress(f, receiver; scope::NamedTuple = (;))
    previous = _EXECUTION_OBSERVATION[]
    sequence = previous.sequence === nothing ? Threads.Atomic{Int}(0) : previous.sequence
    context = (receiver, scope = merge(previous.scope, scope), quiet = previous.quiet, sequence)
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
The execution owner publishes suspended observation before entering and times
only the desired computation inside `f`. Exceptions and interrupts propagate.
"""
function with_performance_sample(f)
    context = (receiver = nothing, scope = (;), quiet = true, sequence = nothing)
    return Base.ScopedValues.with(_EXECUTION_OBSERVATION => context) do
        Logging.with_logger(f, Logging.NullLogger())
    end
end

"""
$(TYPEDSIGNATURES)

Evaluate `f(receiver)` in a complete-scan scope. `total` is the owner-known scan
target, or `nothing`; `batch` identifies the number of outputs sharing a call.
Owners report absolute accepted counts after validating returned results. Nested
scopes retain their parent identity and cannot replace its primary counter.
With observation disabled, pass `nothing` directly without constructing reports.
"""
function with_scan_progress(f; total=nothing, batch=1)
    receiver = progress_receiver()
    receiver === nothing && return f(nothing)
    scope = Threads.atomic_add!(receiver[3], 1) + 1
    parent = get(receiver[2], :scan_scope, 0)
    return with_progress_scope(; scan_scope=scope, scan_parent=parent) do
        scoped = progress_receiver()
        report_progress(scoped, (kind=:scan_start, total, batch))
        state = :complete
        try
            return f(scoped)
        catch error
            state = error isa InterruptException ? :interrupted : :failed
            rethrow()
        finally
            report_progress(scoped, (kind=:scan_end, state))
        end
    end
end

public progress_receiver, report_progress, with_progress, with_progress_scope,
       performance_sample_active, with_performance_sample, with_scan_progress
