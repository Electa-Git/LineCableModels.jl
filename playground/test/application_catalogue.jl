@testset "shared application catalogue and developer navigation" begin
    root = normpath(joinpath(@__DIR__, ".."))
    catalogue = LineCableModelsPlayground.ApplicationCatalogue
    entries = catalogue.entries()
    browser = catalogue.public_entries()
    @test allunique(entry.id for entry in entries)
    @test length(entries) == length(browser)
    @test all(entry -> !haskey(entry, :ui), browser)
    @test Set(e.id for e in entries if e.visibility == :public) == Set(("ichqp-showcase", "cable-study"))
    @test only(filter(e -> e.id == "ichqp-showcase", entries)).requirements ==
        only(filter(e -> e.id == "cable-study", entries)).requirements
    @test length(filter(e -> e.kind == :presentation, entries)) == 3
    for entry in entries
        @test startswith(entry.entrypoint, "/") && !startswith(entry.entrypoint, "//")
        if entry.kind == :presentation
            source = joinpath(root, replace(entry.entrypoint[2:end], ".html" => ".qmd"))
            @test isfile(source)
        end
    end
    home = read(joinpath(root, "index.qmd"), String)
    @test occursin("ENGINEERING SHOWCASE", home)
    @test !occursin("NEW FOUNDATION", home)
    @test occursin("NEW FOUNDATION", read(joinpath(root, "dev", "index.qmd"), String))
    @test occursin("Workbench foundation", read(joinpath(root, "dev", "workbench.qmd"), String))
    chooser = read(joinpath(root, "assets", "application-catalogue.js"), String)
    @test occursin("application-catalogue.json", chooser)
    @test occursin("request_id:requestId", chooser)
    @test !occursin("innerHTML", chooser)
    controller = read(joinpath(root, "_extensions", "lcm-deck", "deck-controller.js"), String)
    @test occursin("query.get('lcm-run')", controller)
    @test occursin("'/applications/runs/' + runId + authored", controller)
    @test occursin("query.has('lcm-run') && !runId", controller)
end
