using HTTP: HTTP
using JSON: JSON
using Playground: Playground
using Sockets: Sockets
using Test: @test, @testset
using quarto_jll: quarto

root = pkgdir(Playground)

@testset "Playground" begin
    @testset "julia -m Playground serves the landing page and the presentations" begin
        port, listener = Sockets.listenany(Sockets.localhost, 18094)
        close(listener)
        address = "http://127.0.0.1:$port"
        server = run(`$(Base.julia_cmd()) --project=$root -m Playground --port $port`; wait = false)
        try
            timedwait(600; pollint = 1) do
                process_exited(server) || try
                    HTTP.get(address; retry = false).status == 200
                catch
                    false
                end
            end
            @test process_running(server)

            landing = String(HTTP.get(address).body)
            @test occursin("href=\"./presentations/specimen/index.html\"", landing)
            @test occursin("href=\"./presentations/starter/index.html\"", landing)
            @test occursin("LCM presentation hostile specimen", landing)
            @test occursin("Presentation title", landing)

            specimen = String(HTTP.get("$address/presentations/specimen/").body)
            slides = split(specimen, "<section ")[2:end]
            slots = Dict(
                "full-canvas" => 1,
                "balanced" => 2,
                "feature-sidebar" => 2,
                "top-split" => 3,
                "dashboard-grid" => 4,
                "media-story" => 2
            )
            layouts = [layout.captures[1]
                       for slide in slides
                       for layout in eachmatch(r"data-lcm-layout=\"([^\"]+)\"", slide)]
            @test length(layouts) == 8
            @test Set(layouts) == Set(keys(slots))
            for slide in slides, layout in eachmatch(r"data-lcm-layout=\"([^\"]+)\"", slide)

                @test count(r"class=\"lcm-slot[ \"]", slide) == slots[layout.captures[1]]
            end
            notes = [note.captures[1]
                     for slide in slides
                     for note in eachmatch(r"data-lcm-math-target=\"([^\"]+)\"", slide)]
            @test Set(notes) == Set(["impedance-imag", "admittance-imag"])
            for slide in slides,
                note in eachmatch(r"data-lcm-math-target=\"([^\"]+)\"", slide)

                @test occursin("\\cssId{$(note.captures[1])}", slide)
            end

            cells = JSON.parse(String(HTTP.get("$address/presentations/specimen/cells.json").body))
            source = read(joinpath(root, "presentations", "specimen", "index.qmd"), String)
            sources = [cell.captures[1]
                       for cell in
                           eachmatch(r"```\{\.julia \.lcm-cell[^}]*\}\n(.*?)\n```"s, source)]
            @test length(sources) == 2
            @test [cell["source"] for cell in cells] == sources
            @test [cell["index"] for cell in cells] == 1:length(sources)
            frames = [frame.match
                      for frame in
                          eachmatch(r"<iframe [^>]*data-lcm-src=\"live/cells/\d+\"[^>]*>", specimen)]
            @test length(frames) == length(cells)
            @test !any(frame -> occursin(" src=", frame), frames)

            states = JSON.parse(String(HTTP.get("$address/api/presentations").body))
            @test states["specimen"] == Dict("state" => "stopped", "starts" => 0)
            @test states["starter"] == Dict("state" => "stopped", "starts" => 0)

            live = "$address/presentations/specimen/live"
            waiting = HTTP.get("$live/cells/1"; status_exception = false)
            @test waiting.status == 503
            @test occursin("lc-widget-shell", String(copy(waiting.body)))
            @test occursin("Julia starting", String(waiting.body))
            @test JSON.parse(String(HTTP.get("$address/api/presentations").body))["specimen"] ==
                  Dict("state" => "starting", "starts" => 1)
            timedwait(600; pollint = 1) do
                JSON.parse(String(HTTP.get("$address/api/presentations").body))["specimen"]["state"] !=
                "starting"
            end
            page = String(HTTP.get("$live/cells/1").body)
            @test occursin("lc-widget-shell", page)
            @test JSON.parse(String(HTTP.get("$address/api/presentations").body))["specimen"] ==
                  Dict("state" => "running", "starts" => 1)
            brand = match(r"href=\"(/presentations/specimen/live/assets/[0-9a-f]+-brand\.css)\"", page)
            @test HTTP.get("$address$(brand.captures[1])").status == 200
            @test HTTP.get("$address/presentations/unknown/live/cells/1"; status_exception = false).status ==
                  404
            closing = try
                HTTP.WebSockets.open(ws -> foreach(identity, ws),
                    "ws://127.0.0.1:$port/presentations/specimen/live/0")
                1000
            catch error
                error.message.code
            end
            @test closing == 4404

            @test HTTP.post("$address/api/presentations/specimen/restart").status == 202
            @test HTTP.post("$address/api/presentations/unknown/restart"; status_exception = false).status ==
                  404
            timedwait(600; pollint = 1) do
                JSON.parse(String(HTTP.get("$address/api/presentations").body))["specimen"]["state"] !=
                "starting"
            end
            @test JSON.parse(String(HTTP.get("$address/api/presentations").body))["specimen"] ==
                  Dict("state" => "running", "starts" => 2)

            worker = parse(Int, readchomp(`pgrep -P $(getpid(server))`))
            run(`kill -9 $worker`)
            timedwait(10; pollint = 0.1) do
                JSON.parse(String(HTTP.get("$address/api/presentations").body))["specimen"]["state"] ==
                "stopped"
            end
            @test JSON.parse(String(HTTP.get("$address/api/presentations").body))["specimen"]["state"] ==
                  "stopped"
            @test HTTP.get("$address/presentations/specimen/").status == 200
            stopped = HTTP.get("$live/cells/1"; status_exception = false)
            @test stopped.status == 503
            @test occursin("Julia stopped", String(stopped.body))
            @test JSON.parse(String(HTTP.get("$address/api/presentations").body))["specimen"] ==
                  Dict("state" => "stopped", "starts" => 2)

            HTTP.post("$address/api/presentations/specimen/restart")
            timedwait(600; pollint = 1) do
                JSON.parse(String(HTTP.get("$address/api/presentations").body))["specimen"]["state"] !=
                "starting"
            end
            @test JSON.parse(String(HTTP.get("$address/api/presentations").body))["specimen"]["state"] ==
                  "running"
            worker = parse(Int, readchomp(`pgrep -P $(getpid(server))`))
            kill(server)
            wait(server)
            timedwait(() -> !isdir("/proc/$worker"), 30; pollint = 0.1)
            @test !isdir("/proc/$worker")
        finally
            kill(server)
            wait(server)
        end
    end

    @testset "a presentation's process refuses a request without its secret" begin
        presentation = Playground.Presentation("specimen")
        try
            cell = "http://127.0.0.1:$(presentation.port)/cells/1"
            @test HTTP.get(cell; status_exception = false).status == 403
            @test HTTP.get(cell, ["Playground-Secret" => presentation.secret]).status == 200
        finally
            close(presentation)
        end
        @test !isopen(presentation)
    end

    @testset "julia -m Playground refuses a presentation name the relay cannot carry" begin
        mktempdir(joinpath(root, "presentations"); prefix = "Specimen_") do folder
            touch(joinpath(folder, "index.qmd"))
            output = IOBuffer()
            command = `$(Base.julia_cmd()) --project=$root -m Playground --no-render`
            process = run(pipeline(ignorestatus(command); stdout = output, stderr = output))
            @test process.exitcode != 0
            @test occursin(
                "the presentation '$(basename(folder))' needs a name of lowercase letters, digits and hyphens",
                String(take!(output)))
        end
    end

    @testset "deck.lua refuses a deck it cannot lay out" begin
        refusals = [
            "unknown-layout.qmd" => "lcm-deck: unknown layout 'mosaic'",
            "missing-anchor.qmd" => "lcm-deck: missing \\cssId anchor for 'impedance-imag'",
            "anchor-elsewhere.qmd" => "lcm-deck: math note and anchor 'impedance-imag' must be on the same slide",
            "wrong-slots.qmd" => "lcm-deck: layout 'balanced' requires 2 direct .lcm-slot children; found 1",
            "julia-engine.qmd" => "lcm-deck: slide 'A slide': a {julia} block runs in Quarto's Julia engine; write a cell as {.julia .lcm-cell}",
            "cell-attribute.qmd" => "lcm-deck: slide 'A slide': a cell takes the attributes autorun and echo; found 'eval'",
            "cell-value.qmd" => "lcm-deck: slide 'A slide': the attribute autorun of a cell is true or false; found 'yes'",
            "cell-outside-slot.qmd" => "lcm-deck: slide 'A slide': a cell belongs inside a .lcm-slot"
        ]
        mktempdir(@__DIR__) do folder
            cp(joinpath(root, "_extensions"), joinpath(folder, "_extensions"))
            cp(joinpath(root, "assets"), joinpath(folder, "assets"))
            for (fixture, sentence) in refusals
                cp(joinpath(@__DIR__, "fixtures", fixture), joinpath(folder, fixture))
                output = IOBuffer()
                render = `$(quarto()) render $(joinpath(folder, fixture)) -M engine:markdown`
                process = run(pipeline(ignorestatus(render); stdout = output, stderr = output))
                @test process.exitcode != 0
                @test occursin(sentence, String(take!(output)))
            end
        end
    end
end
