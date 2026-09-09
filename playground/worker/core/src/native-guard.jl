# Fixed pre-import guard. The agent supplies these values from the approved
# profile and exact receipt; no browser command or code is accepted here.
import LineCableModelsExecutionCore

const LCM_NATIVE_ISOLATION = let
    failure,report = nothing,nothing
    try
        E = LineCableModelsExecutionCore
        limits = E.ExecutorLimits(parse(Float64,ENV["LCM_NATIVE_CPUS"]),
            parse(Int,ENV["LCM_NATIVE_MEMORY_BYTES"]),parse(Int,ENV["LCM_NATIVE_PIDS"]),
            parse(Int,ENV["LCM_NATIVE_SCRATCH_BYTES"]))
        identity = E.NativeIdentity(parse(Int,ENV["LCM_NATIVE_UID"]),parse(Int,ENV["LCM_NATIVE_GID"]),
            ENV["LCM_NATIVE_CGROUP"])
        report = E.verify_native_isolation(limits,identity)
    catch error
        failure = error isa LineCableModelsExecutionCore.IsolationError ? error.code : :guard_configuration_unverified
    end
    if failure !== nothing
        println(stderr,"LCM native entry denied: ",failure)
        exit(78)
    end
    report
end
