@testitem "Core / test discovery rejects excluded setups malformed source and empty selection" tags=[:unit] begin
    helper = joinpath(pkgdir(LineCableModels), "test", "support", "runner.jl")
    project = dirname(Base.active_project())
    function invoke(root; selected=true)
        code = "include($(repr(helper))); ValidationTestRunner.run_tests($(repr(root)); filter=(_ -> $selected))"
        output = IOBuffer()
        process = run(pipeline(ignorestatus(`$(Base.julia_cmd()) --startup-file=no --project=$project -e $code`);
            stdout=output, stderr=output); wait=false)
        status = timedwait(() -> process_exited(process), 60.0; pollint=0.05)
        if status === :timed_out
            kill(process)
            wait(process)
            error("test-discovery subprocess exceeded 60 seconds")
        end
        wait(process)
        return success(process), String(take!(output))
    end
    mktempdir() do root
        write(joinpath(root, "Project.toml"), "[deps]\n")
        write(joinpath(root, "JuliaTestItems.toml"), "config-version = 1\nexclude = [\"excluded/**\"]\n")
        excluded = mkpath(joinpath(root, "excluded"))
        marker = joinpath(root, "unexpected-setup")
        write(joinpath(excluded, "setup.jl"), """
            @testmodule ChosenSetup begin
                write($(repr(marker)), "excluded")
                const answer = :wrong
            end
            @testmodule ExcludedOnly begin
                write($(repr(marker)), "excluded-only")
            end
            """)
        active = joinpath(root, "active.jl")
        write(active, """
            @testmodule ChosenSetup begin
                const answer = :current
            end
            @testitem "included contract" setup=[ChosenSetup] begin
                @test ChosenSetup.answer === :current
            end
            """)
        ok, output = invoke(root)
        @test ok
        @test occursin("Selected 1 maintained test items", output)
        @test !isfile(marker)

        ok, output = invoke(root; selected=false)
        @test !ok
        @test occursin("No test items selected", output)

        write(active, """
            @testitem "missing setup" setup=[ExcludedOnly] begin
                @test true
            end
            """)
        ok, output = invoke(root)
        @test !ok
        @test occursin("Test setup ExcludedOnly is not defined", output)
        @test !isfile(marker)

        write(active, "@testitem \"valid body\" begin\n@test true\nend\nend\n")
        ok, output = invoke(root)
        @test !ok
        @test occursin("ParseError", output)
        @test occursin("active.jl", output)

        write(joinpath(root, "JuliaTestItems.toml"), "config-version = 1\nexlcude = [\"excluded/**\"]\n")
        ok, output = invoke(root)
        @test !ok
        @test occursin("Invalid test discovery keys", output)
    end
end
