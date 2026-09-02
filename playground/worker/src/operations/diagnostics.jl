function validate_echo(parameters)
    haskey(parameters, "message") || throw(PermanentOperationError(
        "invalid_input",
        "system.echo requires a message"
    ))
    return parameters
end

function execute_echo(context, parameters)
    progress!(context, 0.5, "echo"; message="Echoing validated data")
    progress!(context, 1.0, "complete")
    return Dict("echo" => parameters)
end

function validate_delay(parameters)
    seconds = bounded_real(parameters, "seconds", 0.0, 30.0)
    return Dict{String,Any}("seconds" => seconds)
end

function execute_delay(context, parameters)
    seconds = parameters["seconds"]
    steps = max(1, ceil(Int, seconds / 0.1))
    for step in 1:steps
        check_canceled(context)
        seconds > 0 && sleep(seconds / steps)
        progress!(context, step / steps, "delay";
            message="Delay step $step of $steps")
    end
    return Dict("elapsed_seconds" => seconds)
end

function validate_progress(parameters)
    steps = bounded_integer(parameters, "steps", 1, 100)
    interval = get(parameters, "interval_seconds", 0.05)
    interval isa Real || throw(PermanentOperationError(
        "invalid_input",
        "interval_seconds must be numeric"
    ))
    0 <= interval <= 1 || throw(PermanentOperationError(
        "invalid_input",
        "interval_seconds must be between zero and one"
    ))
    steps * interval <= 55 || throw(PermanentOperationError(
        "invalid_input",
        "diagnostic progress duration cannot exceed 55 seconds"
    ))
    return Dict{String,Any}(
        "steps" => steps,
        "interval_seconds" => Float64(interval),
    )
end

function execute_progress(context, parameters)
    steps = parameters["steps"]
    interval = parameters["interval_seconds"]
    for step in 1:steps
        check_canceled(context)
        interval > 0 && sleep(interval)
        log!(context, "completed diagnostic step $step")
        progress!(context, step / steps, "diagnostic";
            message="Step $step of $steps")
    end
    return Dict("steps_completed" => steps)
end

function validate_fail(parameters)
    category = string(get(parameters, "category", "requested_failure"))
    message = string(get(parameters, "message", "Diagnostic failure requested"))
    return Dict{String,Any}("category" => category, "message" => message)
end

function execute_fail(_, parameters)
    throw(PermanentOperationError(parameters["category"], parameters["message"]))
end

function execute_warning(_, parameters)
    @warn parameters["message"]
    return Dict{String,Any}("message" => parameters["message"])
end

function register_diagnostic_operations!(registry)
    register!(registry, OperationSpec(
        "system.echo",
        validate_echo,
        execute_echo;
        timeout_seconds=5,
        cache_policy=:content,
        execution_mode=:daemon
    ))
    register!(registry, OperationSpec(
        "system.delay",
        validate_delay,
        execute_delay;
        timeout_seconds=35,
        cache_policy=:none,
        execution_mode=:daemon
    ))
    register!(registry, OperationSpec(
        "system.executor_delay",
        validate_delay,
        execute_delay;
        timeout_seconds=5,
        cache_policy=:none,
        execution_mode=:supervised
    ))
    register!(registry, OperationSpec(
        "system.executor_warning",
        validate_echo,
        execute_warning;
        timeout_seconds=5,
        cache_policy=:none,
        execution_mode=:supervised
    ))
    register!(registry, OperationSpec(
        "system.progress",
        validate_progress,
        execute_progress;
        timeout_seconds=60,
        cache_policy=:none,
        execution_mode=:daemon
    ))
    register!(registry, OperationSpec(
        "system.fail",
        validate_fail,
        execute_fail;
        timeout_seconds=5,
        cache_policy=:none,
        execution_mode=:daemon
    ))
    return registry
end
