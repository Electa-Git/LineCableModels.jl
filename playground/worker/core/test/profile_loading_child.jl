using LineCableModelsExecutionCore
import LineCableModelsExecutionCore: operation_registry, load_profile!, validate_preparation, prepare!, cleanup!

struct DeferredFixtureProfile <: AbstractScientificProfile end
const loads = Ref(0)
function load_profile!(::DeferredFixtureProfile)
    loads[] += 1
    only(ARGS) == "fail" && throw(RetryableOperationError("Fixture environment load failed"))
    # Exercise methods introduced after the command reader's world began.
    isdefined(@__MODULE__, :loaded_value) || Core.eval(@__MODULE__, :(loaded_value() = 7))
    return nothing
end
operation_registry(::DeferredFixtureProfile) = register!(OperationRegistry(),
    OperationSpec("fixture.loaded", identity,
        (_, _) -> Dict{String,Any}("loads"=>loads[], "value"=>loaded_value()); execution_mode=:supervised))
validate_preparation(::DeferredFixtureProfile, parameters::Dict{String,Any}) = parameters
prepare!(::DeferredFixtureProfile, ::ExecutionContext, ::Dict{String,Any}) =
    PreparedWorkload(Dict{String,Any}("loads"=>loads[], "value"=>loaded_value()))
cleanup!(::DeferredFixtureProfile, cache::PreparedResourceCache) = clear_prepared!(cache)
profile_executor_main(DeferredFixtureProfile())
