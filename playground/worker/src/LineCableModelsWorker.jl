"""
    LineCableModelsWorker

Consume typed playground jobs and execute only registered scientific
operations. The module owns engine dependencies; the publisher does not.
"""
module LineCableModelsWorker

using AWS
using AWSS3
using Dates
using JSON3
using LineCableModelsPlaygroundProtocol
using Logging
using NATS
using SHA
using Sockets
using TOML
using URIs
using UUIDs

include("OperationRegistry.jl")
include("Cache.jl")
include("Artifacts.jl")
include("Executor.jl")
include("operations/diagnostics.jl")
include("operations/linecablemodels.jl")
include("operations/powerimpedance.jl")
include("Consumer.jl")

export OperationRegistry,
    OperationSpec,
    WorkerConfig,
    default_registry,
    executor_main,
    execute_operation,
    register!,
    run_worker

function default_registry()
    registry = OperationRegistry()
    register_diagnostic_operations!(registry)
    register_linecablemodels_operations!(registry)
    register_powerimpedance_operations!(registry)
    return registry
end

function load_scientific_packages!()
    isdefined(@__MODULE__, :LineCableModels) || Core.eval(
        @__MODULE__,
        :(import LineCableModels)
    )
    isdefined(@__MODULE__, :PowerImpedance) || Core.eval(
        @__MODULE__,
        :(import PowerImpedance)
    )
    return nothing
end

function usage(io::IO=stdout)
    println(io, "LineCableModels worker")
    println(io)
    println(io, "Usage:")
    println(io, "  lcm worker start [--nats URL] [--id NAME] [--capacity N]")
    return nothing
end

function cli_config(arguments)
    nats_url = get(ENV, "NATS_CONNECT_URL", "nats://127.0.0.1:4222")
    worker_id = get(ENV, "LCM_WORKER_ID", "$(Sockets.gethostname())-$(getpid())")
    capacity = parse(Int, get(ENV, "LCM_WORKER_CAPACITY", "1"))
    index = 1
    while index <= length(arguments)
        argument = arguments[index]
        argument in ("-h", "--help") && return nothing
        if argument in ("--nats", "--id", "--capacity")
            index == length(arguments) && throw(ArgumentError("$argument requires a value"))
            value = arguments[index + 1]
            if argument == "--nats"
                nats_url = value
            elseif argument == "--id"
                worker_id = value
            else
                parsed_capacity = tryparse(Int, value)
                isnothing(parsed_capacity) && throw(ArgumentError(
                    "--capacity must be an integer"
                ))
                capacity = parsed_capacity
            end
            index += 2
            continue
        end
        throw(ArgumentError("unknown worker option: $argument"))
    end
    return WorkerConfig(; nats_url, worker_id, capacity)
end

function run_cli(arguments)
    arguments in (
        ["-h"],
        ["--help"],
        ["worker"],
        ["worker", "-h"],
        ["worker", "--help"],
    ) && return usage()
    length(arguments) >= 2 || throw(ArgumentError("expected `worker start`"))
    arguments[1:2] == ["worker", "start"] || throw(ArgumentError(
        "expected `worker start`"
    ))
    config = cli_config(arguments[3:end])
    isnothing(config) && return usage()
    return run_worker(config)
end

function (@main)(arguments)
    try
        run_cli(arguments)
        return 0
    catch error
        error isa InterruptException && return 0
        println(stderr, "worker error: ", sprint(showerror, error))
        return 2
    end
end

end
