"""
    LineCableModelsLineParameters

Provide the existing cable and line-parameter adapters in a numerical environment
that does not depend on PowerImpedance, Bonito, or NATS.
"""
module LineCableModelsLineParameters

using LineCableModelsExecutionCore
import LineCableModelsExecutionCore: operation_registry, load_profile!, validate_preparation, prepare!, cleanup!
include("../../../src/operations/linecablemodels.jl")

"""Select the isolated line-parameter operation and representative-workload hooks."""
struct LineParameterProfile <: AbstractScientificProfile end

function load_profile!(::LineParameterProfile)
    isdefined(@__MODULE__, :LineCableModels) || Core.eval(@__MODULE__, :(import LineCableModels))
    return nothing
end

function operation_registry(::LineParameterProfile)
    return register_linecablemodels_operations!(OperationRegistry())
end

validate_preparation(::LineParameterProfile, parameters::Dict{String,Any}) = validate_line_parameters(parameters)

function prepare!(::LineParameterProfile, context::ExecutionContext, parameters::Dict{String,Any})
    # Two endpoint frequencies exercise the actual scientific route without
    # claiming that all future specializations or a full scan have been warmed.
    frequencies = unique([first(parameters["frequencies_hz"]), last(parameters["frequencies_hz"])])
    representative = merge(parameters, Dict{String,Any}("frequencies_hz" => frequencies))
    result = compute_line_parameters(context, representative)
    return PreparedWorkload(Dict{String,Any}("operation" => "line.frequency_scan", "frequencies_hz" => frequencies,
        "workload_input_hash" => input_hash("line.frequency_scan", representative),
        "result_hash" => input_hash("line.frequency_scan", result),
        "preparation_kind" => "representative_execution"))
end

cleanup!(::LineParameterProfile, cache::PreparedResourceCache) = clear_prepared!(cache)

"""Run one disposable line-parameter executor on its parent's local input/output channel."""
main() = profile_executor_main(LineParameterProfile())

end
