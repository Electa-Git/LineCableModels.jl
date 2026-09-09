@testset "default catalogue is passive and launch capability is explicit" begin
    registry = default_applications()
    @test length(registry.definitions) == 6
    @test length(registry.applications) == 6
    @test haskey(registry.definitions, "ichqp-showcase")
    @test haskey(registry.applications, "ichqp-showcase")
    @test haskey(registry.applications, "cable-study")
    @test registry.definitions["starter-deck"].entrypoint == "/presentations/starter.html"
    @test registry.definitions["toolkit-gallery"].entry_surface == :published
    @test_throws ArgumentError ApplicationDefinition("bad", "Bad", :workbench, "/bad"; entry_surface=:guessed)
    context = HostContext(uuid4(), "/unused/", "/unused/ready.json")
    @test all(app -> "--compiled-modules=existing" in ui_command(app, context).exec,
        values(registry.applications))
    @test all(app -> RequiredInterfaces.check_interface_implemented(AbstractApplication, typeof(app)) === true,
        values(registry.applications))
    @test_throws ArgumentError register!(registry, registry.definitions["ichqp-showcase"])
    mktempdir() do dir
        # Unimplemented passive registrations still cannot allocate a host.
        register!(registry, ApplicationDefinition("passive-fixture", "Passive", :presentation, "/passive.html"))
        store = RuntimeStore(joinpath(dir, "runtime.sqlite"))
        supervisor = UIHostSupervisor(store, registry, joinpath(dir, "hosts"))
        try
            @test_throws AccessDenied start_ui!(supervisor, Principal("alice"), "passive-fixture")
            @test isempty(list_runs(store, Principal("alice")))
            @test isempty(supervisor.handles)
        finally
            close(supervisor); close(store)
        end
    end
end

@testset "published site exact routes and confinement" begin
    mktempdir() do dir
        write(joinpath(dir, "index.html"), "public")
        mkpath(joinpath(dir, "presentations"))
        write(joinpath(dir, "presentations", "index.html"), "catalogue")
        site = PublishedSite(dir)
        @test site.files["/"] == joinpath(dir, "index.html")
        @test site.files["/presentations"] == site.files["/presentations/"]
        @test !haskey(site.files, "/../index.html")
        symlink(joinpath(dir, "index.html"), joinpath(dir, "alias.html"))
        @test_throws ArgumentError PublishedSite(dir)
        rm(joinpath(dir, "alias.html"))
        mkpath(joinpath(dir, "runtime"))
        write(joinpath(dir, "runtime", "index.html"), "must not bypass authentication")
        @test_throws ArgumentError PublishedSite(dir)
    end
end

@testset "runtime CLI refuses ambiguous options" begin
    @test_throws ArgumentError runtime_cli(["runtime", "start"])
    @test_throws ArgumentError runtime_cli(["runtime", "unexpected"])
    @test_throws ArgumentError runtime_cli(["runtime", "check", "--xray"])
    @test_throws ArgumentError runtime_cli(["runtime", "check", "--config", "one", "--config", "two"])
end
