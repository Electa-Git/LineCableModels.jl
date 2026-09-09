using LineCableModelsExecutionCore
import LineCableModelsExecutionCore: operation_registry, validate_preparation, prepare!, cleanup!

struct CacheProfile <: AbstractScientificProfile
    builds::Base.RefValue{Int}
end
operation_registry(::CacheProfile) = register!(OperationRegistry(),
    OperationSpec("fixture.echo", identity, (_, parameters) -> parameters; execution_mode=:supervised))
validate_preparation(::CacheProfile, parameters::Dict{String,Any}) = parameters
function prepare!(profile::CacheProfile, context::ExecutionContext, parameters::Dict{String,Any})
    profile.builds[] += 1
    model = prepare_resource!(context.prepared_cache, "model") do
        fill(UInt8(1), get(parameters,"bytes",32))
    end
    return PreparedWorkload(Dict{String,Any}("bytes"=>length(model)); resources=("model",))
end
cleanup!(::CacheProfile, cache::PreparedResourceCache) = clear_prepared!(cache)
