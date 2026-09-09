# Explicit test profile, not a production resource-policy adapter.
ccall(:getppid, Cint, ()) == parse(Int, ENV["LCM_EXECUTOR_PARENT_PID"]) || exit(70)
using LineCableModelsExecutionCore
import LineCableModelsExecutionCore: operation_registry, validate_preparation, prepare!, cleanup!
struct OwnedFixtureProfile <: AbstractScientificProfile end
function operation_registry(::OwnedFixtureProfile)
    registry = OperationRegistry()
    register!(registry, OperationSpec("fixture.echo", identity, (_, p) -> p; execution_mode=:supervised,schema_version="1.2"))
    register!(registry, OperationSpec("fixture.large", p -> begin
        count=get(p,"count",0)
        count isa Integer && !(count isa Bool) && 1<=count<=100000 || throw(ArgumentError("bounded count required"))
        p
    end, (_, p) -> Dict{String,Any}("values"=>collect(1:p["count"])); execution_mode=:supervised))
    register!(registry, OperationSpec("fixture.delay", identity, (context, p) -> begin
        progress!(context, 0.1, "running")
        log!(context, "fixture private output must not enter status")
        sleep(p["seconds"])
        Dict{String,Any}("done"=>true)
    end; execution_mode=:supervised))
    register!(registry, OperationSpec("fixture.evict", identity, (context, _) -> begin
        delete!(context.prepared_cache.entries, "model")
        Dict{String,Any}("evicted"=>true)
    end; execution_mode=:supervised))
    register!(registry, OperationSpec("fixture.fail", identity, (_, _) ->
        throw(PermanentOperationError("validation", "/private/fixture secret diagnostic"));
        execution_mode=:supervised))
    return registry
end
validate_preparation(::OwnedFixtureProfile, p::Dict{String,Any}) = p
function prepare!(::OwnedFixtureProfile, context::ExecutionContext, p::Dict{String,Any})
    progress!(context, 0.2, "preparing")
    sleep(get(p, "seconds", 0.0))
    model = prepare_resource!(() -> fill(UInt8(1), 32), context.prepared_cache, "model")
    PreparedWorkload(Dict{String,Any}("bytes"=>length(model)); resources=("model",))
end
cleanup!(::OwnedFixtureProfile, cache::PreparedResourceCache) = clear_prepared!(cache)
profile_executor_main(OwnedFixtureProfile(); cache=PreparedResourceCache(ttl_seconds=parse(Float64, only(ARGS))))
