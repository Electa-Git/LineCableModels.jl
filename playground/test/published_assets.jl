@testset "committed publication snapshots" begin
    PA = LineCableModelsPlayground.PublishedAssets
    mktempdir() do directory
        write(joinpath(directory, "index.html"), "old document")
        write(joinpath(directory, "old.css"), "old CSS")
        PA.publish!(directory)
        site = PA.PublishedSite(directory)
        text(route) = String(copy(PA.published_asset(site, route)))
        @test text("/") == "old document"
        write(joinpath(directory, "index.html"), "new document")
        write(joinpath(directory, "new.css"), "new CSS")
        # An unfinished or failed render must not leak any partial changes.
        @test text("/") == "old document"
        @test PA.published_asset(site, "/new.css") === nothing
        rm(joinpath(directory, "old.css"))
        @test text("/old.css") == "old CSS"
        PA.publish!(directory)
        @test text("/") == "new document"
        @test text("/new.css") == "new CSS"
        @test text("/old.css") == "old CSS"
        @test PA.published_asset(site, "/../index.html") === nothing
        @test PA.published_asset(site, "/.lcm-publication") === nothing
        symlink(joinpath(directory, "index.html"), joinpath(directory, "alias.html"))
        @test_throws ArgumentError PA.publish!(directory)
        @test text("/") == "new document"
        rm(joinpath(directory, "alias.html"))
        write(joinpath(directory, "index.html"), "third document")
        PA.publish!(directory)
        @test text("/") == "third document"
        @test PA.published_asset(site, "/old.css") === nothing
        mkpath(joinpath(directory, "runtime"))
        write(joinpath(directory, "runtime", "index.html"), "reserved")
        @test_throws ArgumentError PA.publish!(directory)
        @test text("/") == "third document"
    end
end

@testset "publication yields to lazily installed live routes" begin
    LCM = LineCableModelsPlayground
    HTTP = LCM.Bonito.HTTP
    router = LCM.Bonito.HTTPServer
    mktempdir() do directory
        write(joinpath(directory, "index.html"), "home")
        LCM.PublishedAssets.publish!(directory)
        routes = router.Routes()
        LCM.register_static_site_routes!(routes, directory)
        # Bonito installs an asset regex only after rendering its first App.
        router.route!(routes, r"^/assets/[a-f0-9]{40}-" => (_ -> HTTP.Response(200, "live asset")))
        @test first(last(routes.table)) isa LCM.PublishedRoute
        request(path) = router.delegate(routes, nothing, HTTP.Request("GET", path))
        @test String(request("/assets/" * repeat("a", 40) * "-session.js").body) == "live asset"
        write(joinpath(directory, "new.css"), "new CSS")
        LCM.PublishedAssets.publish!(directory)
        @test String(request("/new.css").body) == "new CSS"
        @test request("/unknown.js").status == 404
    end
end
