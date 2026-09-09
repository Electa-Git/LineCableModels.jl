"""
    LineCableModelsPowerFlow

Provide the existing OHL/UGC power-flow and impedance adapters without importing
the line-parameter package, Bonito or NATS. Numerical preparation remains explicit.
"""
module LineCableModelsPowerFlow

using LineCableModelsExecutionCore
import LineCableModelsExecutionCore: operation_registry, load_profile!, validate_preparation, prepare!, cleanup!
include("../../../src/operations/powerimpedance.jl")

"""Select the isolated OHL/UGC power-flow preparation and impedance hooks."""
struct PowerFlowProfile <: AbstractScientificProfile end

function load_profile!(::PowerFlowProfile)
    isdefined(@__MODULE__, :PowerImpedance) || Core.eval(@__MODULE__, :(import PowerImpedance))
    return nothing
end

operation_registry(::PowerFlowProfile) = register_powerimpedance_operations!(OperationRegistry())
validate_preparation(::PowerFlowProfile, parameters::Dict{String,Any}) = validate_powerflow_spec(parameters)

function prepare!(::PowerFlowProfile, context::ExecutionContext, parameters::Dict{String,Any})
    result = execute_powerflow_prepare(context, parameters)
    return PreparedWorkload(Dict{String,Any}("operation" => "powerflow.prepare",
        "workload_input_hash" => input_hash("powerflow.prepare", parameters),
        "prepared_resource_key" => result["prepared_resource_key"],
        "preparation_kind" => "solved_and_linearized_model"); resources=(result["prepared_resource_key"],))
end

cleanup!(::PowerFlowProfile, cache::PreparedResourceCache) = clear_prepared!(cache)

"""Run one disposable power-flow executor on its parent's local input/output channel."""
main() = profile_executor_main(PowerFlowProfile())

end
