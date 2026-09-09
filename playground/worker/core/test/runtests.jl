using Test, Dates, LineCableModelsExecutionCore
const E = LineCableModelsExecutionCore
include("profile_fixture.jl")

function context(; cache=PreparedResourceCache())
    ExecutionContext("fixture", CancellationToken(), (_, _, _) -> nothing, _ -> nothing,
        String[], cache, nothing)
end

include("executor_cleanup.jl")
include("profile_loading.jl")
include("container_isolation.jl")
include("container_entry.jl")
include("native_isolation.jl")
include("native_entry.jl")

@testset "core import and framed IO boundary" begin
    @test !any(m -> nameof(m) in (:LineCableModels,:PowerImpedance,:NATS,:Bonito,:AWS), values(Base.loaded_modules))
    @test E.read_executor_line(IOBuffer("abcd\nrest"); maximum_bytes=4) == "abcd"
    @test E.read_executor_line(IOBuffer("last"); maximum_bytes=4) == "last"
    @test_throws ArgumentError E.read_executor_line(IOBuffer("abcde\n"); maximum_bytes=4)
    @test_throws EOFError E.read_executor_line(IOBuffer())
    @test_throws ArgumentError E.read_executor_line(IOBuffer(); maximum_bytes=0)
    spec = OperationSpec("fixture.echo",identity,(_,p)->p;execution_mode=:supervised)
    supervisor = ExecutorSupervisor(;command=`false`)
    @test_throws ArgumentError execute_supervised!(supervisor,spec,context(),
        Dict{String,Any}("oversize"=>repeat("x",E.EXECUTOR_MAX_LINE_BYTES)))
    @test supervisor.process === nothing && supervisor.generation == 0
    @test_throws ArgumentError OperationSpec("fixture.bad",identity,(_,p)->p;timeout_seconds=true)
end

@testset "preparation cache bounds and cleanup" begin
    for kwargs in ((ttl_seconds=Inf,), (ttl_seconds=true,), (max_entries=0,),
            (max_entries=257,), (max_bytes=0,), (max_bytes=true,))
        @test_throws ArgumentError PreparedResourceCache(; kwargs...)
    end
    clock = Ref(0.0)
    cache = PreparedResourceCache(max_entries=1, max_bytes=1024, clock=() -> clock[])
    gate = Channel{Nothing}(1)
    calls = Ref(0)
    builder() = (calls[] += 1; take!(gate); :value)
    first = @async prepare_resource!(builder, cache, "key")
    yield()
    second = @async prepare_resource!(builder, cache, "key")
    yield()
    @test prepared_status(cache,"key") == :warming
    @test_throws ArgumentError prepare_resource!(() -> :other, cache, "other")
    put!(gate,nothing)
    @test fetch(first) == fetch(second) == :value
    @test calls[] == 1
    @test prepare_resource!(() -> error("must reuse"), cache, "key") == :value
    clock[] = 901
    @test prepared_status(cache,"key") == :cold
    @test_throws TaskFailedException prepare_resource!(() -> ones(UInt8,2048),cache,"large")
    @test isempty([e for e in values(cache.entries) if e.state == :hot])
    @test_throws ArgumentError prepare_resource!(() -> :any,cache,repeat("a",257))
    clear_prepared!(cache)
    late = @async prepare_resource!(builder,cache,"late")
    yield();clear_prepared!(cache);put!(gate,nothing)
    @test_throws TaskFailedException fetch(late)
    @test isempty(cache.entries)
    bounded = PreparedResourceCache(max_bytes=250,max_entries=8)
    prepare_resource!(() -> fill(UInt8(1),100),bounded,"one")
    prepare_resource!(() -> fill(UInt8(2),100),bounded,"two")
    @test sum(e.retained_bytes for e in values(bounded.entries)) <= 250
    @test length(bounded.entries) == 1
end

struct IncompleteProfile <: AbstractScientificProfile end
operation_registry(::IncompleteProfile) = OperationRegistry()

@testset "preparation inspection never extends or invents readiness" begin
    clock = Ref(0.0)
    cache = PreparedResourceCache(ttl_seconds=10, clock=()->clock[])
    c = context(; cache)
    profile = CacheProfile(Ref(0))
    prepared = E.prepare_profile!(profile, c, Dict{String,Any}())
    key = prepared["preparation_input_hash"]
    @test E.preparation_snapshot(cache, key)["ready"]
    @test E.preparation_snapshot(cache, key)["remaining_seconds"] == 10
    clock[] = 4
    @test E.preparation_snapshot(cache, key)["remaining_seconds"] == 6
    @test cache.entries["model"].expires_at == 10
    delete!(cache.entries, "model")
    @test !E.preparation_snapshot(cache, key)["ready"]
    @test !haskey(cache.entries, "preparation:" * key)
    @test_throws ArgumentError E.preparation_snapshot(cache, "invalid")
    supervisor = ExecutorSupervisor(; command=`false`)
    @test_throws E.RetryableOperationError inspect_preparation!(supervisor, c, key)
    @test supervisor.process === nothing && supervisor.generation == 0
    spec = OperationSpec("fixture.echo", identity, (_,p)->p;execution_mode=:supervised)
    @test_throws E.RetryableOperationError execute_supervised!(supervisor, spec, c, Dict{String,Any}(); preparation_key=key)
    @test supervisor.process === nothing && supervisor.generation == 0
end

@testset "readiness evidence depends on retained resources" begin
    @test_throws ArgumentError E.checked_profile_registry(IncompleteProfile())
    profile = CacheProfile(Ref(0))
    @test E.checked_profile_registry(profile) isa OperationRegistry
    @test load_profile!(profile) === nothing && profile.builds[] == 0
    c = context()
    first = E.prepare_profile!(profile,c,Dict{String,Any}())
    @test first["cache_status"] == "miss"
    @test first["evidence"]["bytes"] == 32
    @test length(first["preparation_input_hash"]) == 64
    second = E.prepare_profile!(profile,c,Dict{String,Any}())
    @test second["cache_status"] == "hit" && profile.builds[] == 1
    delete!(c.prepared_cache.entries,"model")
    third = E.prepare_profile!(profile,c,Dict{String,Any}())
    @test third["cache_status"] == "miss" && profile.builds[] == 2
    cleanup!(profile,c.prepared_cache)
    @test isempty(c.prepared_cache.entries)
    @test_throws TaskFailedException E.prepare_profile!(profile,
        context(cache=PreparedResourceCache(max_entries=1)),Dict{String,Any}())
    @test_throws ArgumentError PreparedWorkload(Dict{String,Any}(); resources=("a","a"))
    # A completion that evicts its own model must never become a ready label.
    @test_throws ArgumentError E.prepare_profile!(profile,
        context(cache=PreparedResourceCache(max_bytes=1400)),Dict{String,Any}("bytes"=>1200))
end
