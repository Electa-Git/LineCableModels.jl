import Bonito

@testset "UI compilation workload stays disconnected" begin
    server = Bonito.GLOBAL_SERVER[]
    cleanup_servers = Set(keys(Bonito.SERVER_CLEANUP_TASKS))
    uploads = Set(keys(LineCableModelsPlayground.DEFAULT_UPLOAD_REGISTRY.entries))
    @test LineCableModelsPlayground.precompile_ui_workload() === nothing
    @test Bonito.GLOBAL_SERVER[] === server
    @test Set(keys(Bonito.SERVER_CLEANUP_TASKS)) == cleanup_servers
    @test Set(keys(LineCableModelsPlayground.DEFAULT_UPLOAD_REGISTRY.entries)) == uploads
    forbidden = Set(("LineCableModels", "PowerImpedance", "PowerModels", "JuMP", "Ipopt"))
    @test isempty(intersect(forbidden, Set(id.name for id in keys(Base.loaded_modules))))
end
