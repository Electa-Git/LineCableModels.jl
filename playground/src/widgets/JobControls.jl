"""
    JobPanel(client; operation, parameters, input_observables=(),
             session_id=string(uuid4()), title="Calculation", priority="normal",
             result_handler=nothing)

Create reusable asynchronous calculation controls. `parameters` is a function
or observable returning JSON-compatible values. Constructing the panel does not
submit work; only the Run or Retry callbacks do so.

Pass every selector observable through `input_observables`; changing any of
them marks the current result dirty without starting a calculation.

When supplied, `result_handler(result)` runs for each successful result. It is
intended to replace data in rendering observables captured by the callback;
figures, axes, plot objects, canvases, and surrounding DOM remain persistent.
An exception in the rendering callback is logged and appended to the panel
console without changing the durable calculation result.
"""
struct JobPanel
    dom::Any
    handle::JobHandle
end

function parameter_snapshot(parameters)
    value = parameters isa Function ? parameters() :
            parameters isa Observable ? parameters[] : parameters
    return normalize_wire(value)
end

function JobPanel(
        client::BrokerClient;
        operation::AbstractString,
        parameters,
        input_observables=(),
        session_id::AbstractString=string(uuid4()),
        title::AbstractString="Calculation",
        priority::AbstractString="normal",
        result_handler::Union{Nothing,Function}=nothing
    )
    operation_subject_token(operation)
    priority in PRIORITIES || throw(ArgumentError("unsupported job priority: $priority"))
    handle = JobHandle()
    observed_inputs = Any[input_observables...]
    parameters isa Observable && push!(observed_inputs, parameters)
    for source in observed_inputs
        source isa Observable || throw(ArgumentError(
            "input_observables must contain Bonito observables"
        ))
        on(source) do _
            mark_dirty!(handle)
            return nothing
        end
    end
    console = ConsoleView(; title="Job activity", status="idle", max_entries=200)
    on(handle.messages) do entries
        console.entries[] = copy(entries)
        console.revision[] += 1
        return nothing
    end
    on(handle.state) do state
        set_console_status!(console, string(state))
        return nothing
    end
    if !isnothing(result_handler)
        on(handle.last_successful_result) do result
            isnothing(result) && return nothing
            try
                result_handler(result)
            catch error
                diagnostic_id = string(uuid4())
                @error "Result renderer failed" job_id=result.job_id diagnostic_id exception=(
                    error,
                    catch_backtrace()
                )
                entry = ConsoleEntry(
                    :error,
                    "Result rendering failed · diagnostic $diagnostic_id"
                )
                handle.messages[] = vcat(handle.messages[], [entry])
            end
            return nothing
        end
    end

    run_disabled = map(handle.state) do state
        state in (:submitting, :queued, :unavailable, :running, :cancelling)
    end
    cancel_disabled = map(handle.state) do state
        !(state in (:queued, :unavailable, :running))
    end
    retry_disabled = map(handle.state) do state
        !(state in (:failed, :canceled, :expired))
    end
    run_button = Button("Run"; disabled=run_disabled)
    cancel_button = Button("Cancel"; disabled=cancel_disabled)
    retry_button = Button("Retry"; disabled=retry_disabled)

    function submit_current!()
        handle.state[] in (:submitting, :queued, :unavailable, :running, :cancelling) &&
            return nothing
        request = new_job_request(
            operation,
            parameter_snapshot(parameters);
            session_id,
            priority
        )
        submit_async!(client, handle, request)
        return nothing
    end
    on(run_button.value) do _
        submit_current!()
    end
    on(retry_button.value) do _
        handle.state[] in (:failed, :canceled, :expired) || return nothing
        submit_current!()
    end
    on(cancel_button.value) do _
        handle.state[] in (:queued, :unavailable, :running) || return nothing
        @async try
            cancel!(client, handle)
        catch error
            submission_failed!(handle, error)
        end
    end

    progress_value = map(handle.progress) do value
        string(round(Int, 100 * clamp(value, 0, 1)))
    end
    progress_style = map(handle.progress) do value
        "--lc-job-progress: $(100 * clamp(value, 0, 1))%"
    end
    state_label = map(handle.state) do state
        uppercasefirst(replace(string(state), '_' => ' '))
    end
    result_status = map(
        handle.state,
        handle.last_successful_result,
        handle.pending_input_hash,
        handle.dirty_inputs
    ) do state, result, pending, dirty
        if isnothing(result)
            return state in (:queued, :running, :unavailable) ?
                "No previous result · request pending" : "No result"
        elseif !isnothing(pending)
            return "Displaying previous result · new calculation pending"
        elseif dirty
            return "Displaying previous result · selectors changed"
        end
        return "Latest result · $(result.cache_status == "hit" ? "cache hit" : "computed")"
    end
    result_preview = map(handle.last_successful_result) do result
        isnothing(result) && return "—"
        if !isnothing(result.inline_result)
            return JSON3.write(result.inline_result)
        end
        artifact = result.artifact
        return DOM.a(
            "Download $(artifact.media_type) artifact ($(artifact.size) bytes)";
            href=artifact.retrieval_reference,
            target="_blank",
            rel="noopener"
        )
    end

    dom = DOM.section(
        DOM.header(
            DOM.div(
                DOM.span("CALCULATION"; class="lc-widget-kicker"),
                DOM.h2(title);
                class="lc-job-heading"
            ),
            WorkerStatus(client, string(operation));
            class="lc-job-header"
        ),
        DOM.div(
            DOM.div(
                DOM.span(state_label; class="lc-job-state"),
                DOM.span(handle.stage; class="lc-job-stage");
                class="lc-job-status-line"
            ),
            DOM.div(
                DOM.span(; class="lc-job-progress-fill"),
                DOM.span(progress_value, "%"; class="lc-job-progress-value");
                Dict{Symbol,Any}(
                    :class => "lc-job-progress",
                    :style => progress_style,
                    :role => "progressbar",
                    Symbol("aria-valuemin") => "0",
                    Symbol("aria-valuemax") => "100",
                    Symbol("aria-valuenow") => progress_value,
                )...
            ),
            DOM.div(
                run_button,
                cancel_button,
                retry_button;
                class="lc-widget-actions lc-job-actions"
            ),
            DOM.div(
                DOM.span(result_status; class="lc-job-result-status"),
                DOM.div(result_preview; class="lc-job-result-preview");
                class="lc-job-result"
            ),
            console;
            class="lc-job-body"
        );
        class="lc-job-panel"
    )
    return JobPanel(dom, handle)
end

Bonito.jsrender(session::Session, panel::JobPanel) = Bonito.jsrender(session, panel.dom)
