@testset "runtime dependency and credential boundaries" begin
    forbidden = Set(["Bonito", "Makie", "WGLMakie", "LineCableModels", "PowerImpedance"])
    dependencies = TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))["deps"]
    @test isempty(intersect(Set(keys(dependencies)), forbidden))
    @test isempty(intersect(Set(id.name for id in keys(Base.loaded_modules)), forbidden))
    context = HostContext(uuid4(), "/unused/", "/tmp/owned-fixture/ready.json")
    withenv("NATS_CONNECT_URL"=>"private-broker", "AWS_SECRET_ACCESS_KEY"=>"private-storage",
            "LCM_PROXY_KEY"=>"private-proxy") do
        environment = RT.child_environment(context, "owned-host-credential")
        @test !any(haskey(environment, key) for key in
            ("NATS_CONNECT_URL", "AWS_SECRET_ACCESS_KEY", "LCM_PROXY_KEY"))
        @test environment["TMPDIR"] == dirname(context.ready_file)
        @test environment["LCM_UI_HOST_KEY"] == "owned-host-credential"
    end
    root = normpath(joinpath(@__DIR__, "..", ".."))
    surface = read(joinpath(root, "runtime", "src", "RunSurface.jl"), String)
    @test occursin("forms.css", surface)
    @test occursin("brand.css", surface)
    for file in ("application-catalogue.css", joinpath("..", "runtime", "ui", "run.css"))
        css = read(joinpath(root, "assets", file), String)
        @test !occursin(r"#[a-fA-F0-9]{3,8}\b", css)
    end
end
