function captured_cli(arguments)
    return mktemp() do _, output
        result = redirect_stdout(output) do
            LineCableModelsPlayground.run_cli(arguments)
        end
        flush(output)
        seekstart(output)
        return result, read(output, String)
    end
end

@testset "normalized CLI" begin
    result, general = captured_cli(String[])
    @test isnothing(result)
    @test occursin("Usage: lcm <feature> <action> [options]", general)
    @test occursin("playground", general)
    @test occursin("worker", general)
    @test occursin("nats", general)
    @test occursin("container", general)
    @test !occursin("linecablemodels ", general)

    _, playground = captured_cli(["playground", "--help"])
    @test occursin("lcm playground build", playground)
    @test occursin("lcm playground start", playground)

    _, start_help = captured_cli(["playground", "start", "--help"])
    @test occursin("--proxy-url URL", start_help)
    @test occursin("--xray", start_help)

    options = LineCableModelsPlayground.parse_start_options([
        "--no-open",
        "--no-render",
        "--xray",
    ])
    @test options.xray
    @test !options.open_browser
    @test !options.render_before_start

    _, nats = captured_cli(["nats", "--help"])
    @test occursin("lcm nats init", nats)
    @test occursin("lcm nats status", nats)
    @test occursin("NATS_CONNECT_URL", nats)
    @test occursin("NATS_ADMIN_URL", nats)

    _, container = captured_cli(["container", "--help"])
    @test occursin("lcm container resolve", container)
    @test occursin("lcm container start", container)
    @test occursin("--runtime RUNTIME", container)
    @test occursin("docker", container)
    @test occursin("podman", container)

    @test_throws ArgumentError LineCableModelsPlayground.run_cli(["runtime", "status"])
    @test_throws ArgumentError LineCableModelsPlayground.run_cli(["unknown", "action"])
end
