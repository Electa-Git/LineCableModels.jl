# Loaded by the fixed Julia command before any numerical profile or interactive
# prompt. Image admission and physical lifetime ownership remain agent duties.
import LineCableModelsExecutionCore

const LCM_CONTAINER_ISOLATION = let
    failure = nothing
    report = nothing
    try
        E = LineCableModelsExecutionCore
        limits = E.ContainerLimits(
            parse(Float64, ENV["LCM_CONTAINER_CPUS"]),
            parse(Int, ENV["LCM_CONTAINER_MEMORY_BYTES"]),
            parse(Int, ENV["LCM_CONTAINER_PIDS"]),
            parse(Int, ENV["LCM_CONTAINER_SCRATCH_BYTES"]))
        report = E.verify_container_isolation(limits)
    catch error
        failure = error isa LineCableModelsExecutionCore.IsolationError ? error.code : :guard_configuration_unverified
    end
    if failure !== nothing
        println(stderr,"LCM container entry denied: ",failure)
        exit(78)
    end
    report
end

# This hook is installed only after kernel isolation verification succeeded.
# Its receipt-bound marker is startup evidence, not writer authorization.
haskey(ENV,"LCM_TERMINAL_READY") && include("terminal-ready.jl")
