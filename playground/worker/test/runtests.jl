using HTTP: HTTP
using Test: @test, @testset
using quarto_jll: quarto

root = normpath(@__DIR__, "..")
playground = normpath(root, "..")

@testset "PlaygroundWorker" begin
    specimen = joinpath(playground, "presentations", "specimen")
    run(pipeline(`$(quarto()) render $(joinpath(specimen, "index.qmd")) -M engine:markdown`; stdout = devnull))
    lists = [joinpath(@__DIR__, "fixtures", "cells.json"),
        joinpath(playground, "_site", "presentations", "specimen", "cells.json")]
    secret = ["Playground-Secret" => "a secret of the test"]
    logs = [Pipe() for list in lists]
    workers = [open(
                   pipeline(
                       addenv(
                           `$(Base.julia_cmd()) --threads=2,1 --project=$project -m PlaygroundWorker $list --prefix /live`,
                           "JULIA_LOAD_PATH" => "@:$root:@stdlib", "PLAYGROUND_SECRET" => secret[1][2]);
                       stderr = log),
                   "r+") for (project, list, log) in zip([root, specimen], lists, logs)]
    foreach(log -> close(log.in), logs)
    printed = [@async(read(log, String)) for log in logs]
    try
        addresses = [replace(readline(worker), "ready " => "http://127.0.0.1:")
                     for worker in workers]
        fixture, presentation = addresses
        timedwait(120; pollint = 0.5) do
            occursin("cell 3 fails", String(HTTP.get("$fixture/cells/3", secret).body))
        end
        pages = [String(HTTP.get("$fixture/cells/$index", secret).body) for index in 1:3]

        @testset "a value shows" begin
            @test occursin(r"<span class=\"text-plain\"[^>]*>2</span>", pages[1])
        end
        @testset "printed output is captured" begin
            @test occursin("printed by cell 2", pages[2])
        end
        @testset "a definition of cell 1 is read by cell 2" begin
            @test occursin(r"<span class=\"text-plain\"[^>]*>5</span>", pages[2])
        end
        @testset "an error shows its message and the line of the cell" begin
            @test occursin("cell 3 fails", pages[3])
            @test occursin("cell 3:2", pages[3])
            @test !occursin("PlaygroundWorker.jl", pages[3])
        end
        @testset "a second session runs the cell again" begin
            @test occursin(r"<span class=\"text-plain\"[^>]*>7</span>",
                String(HTTP.get("$fixture/cells/2", secret).body))
        end
        @testset "the specimen's plot is a figure of each session's own" begin
            plots = [String(HTTP.get("$presentation/cells/4", secret).body)
                     for session in 1:2]
            @test all(page -> occursin("Earth resistivity ρ", page), plots)
            @test all(page -> occursin("WGLMakie", page), plots)
            @test !any(page -> occursin("--lc-danger", page), plots)
        end
        @testset "each session gets a canvas probe of its own" begin
            probes = [String(HTTP.get("$presentation/cells/5", secret).body)
                      for session in 1:2]
            found = [match(r"data-session-id=\"([^\"]+)\"", page) for page in probes]
            @test all(session -> session !== nothing, found)
            @test found[1].captures[1] != found[2].captures[1]
        end
        @testset "the specimen's first cell is a page of the deck's live frame" begin
            page = String(HTTP.get("$presentation/cells/1", secret).body)
            @test occursin("lc-widget-shell", page)
            @test occursin("lcm:child-ready", page)
            @test occursin(r"href=\"/live/assets/[0-9a-f]+-brand\.css\"", page)
        end
        @testset "a request without the secret is refused" begin
            @test HTTP.get("$fixture/cells/1"; status_exception = false).status == 403
            @test HTTP.get("$fixture/cells/1", ["Playground-Secret" => "another"];
                status_exception = false).status == 403
        end
    finally
        for worker in workers
            close(worker.in)
            wait(worker)
        end
    end
    @testset "the worker returns at the end of its input" begin
        @test all(success, workers)
    end
    @testset "the worker logs the version of Bonito it loaded" begin
        @test all(log -> occursin(r"Bonito \d+\.\d+\.\d+", fetch(log)), printed)
    end
end
