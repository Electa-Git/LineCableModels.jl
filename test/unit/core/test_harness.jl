@testitem "Core / test discovery rejects excluded setups malformed source and empty selection" tags=[:unit] begin
    helper = joinpath(pkgdir(LineCableModels), "test", "support", "runner.jl")
    include(helper)
    directory = joinpath(pkgdir(LineCableModels), "test")
    item(file, name, tags) = (; filename=joinpath(directory, file), name, tags)
    ordinary = item("unit/current.jl", "unresolved current contract", [:unit])
    sampling = item("integration/sampling.jl", "current sampling contract", [:integration])
    native = item("extensions/study.jl", "unresolved native study", [:extension, :fem_numerical])
    select(args) = ValidationTestRunner.selection(args, directory)
    @test select(String[])(ordinary)
    @test select(String[])(sampling)
    @test !select(["--list"])(native)
    @test select(["tag:integration"])(sampling)
    @test select(["tag:fem_numerical"])(native)
    @test select(["extensions/study"])(native)
    @test select(["tag:integration", "integration/"])(sampling)
    @test !select(["tag:fem_numerical", "integration/"])(native)
    @test select(["tag:unit", "tag:fem_numerical", "CURRENT", "STUDY"])(ordinary)
    @test !select(["tag:unit"])(merge(ordinary, (; filename=joinpath(dirname(directory), "outside.jl"))))
    @test_throws ErrorException select(["--unknown"])
    @test_throws ErrorException select(["tag:"])
    project = dirname(Base.active_project())
    function invoke(root; selected=true, list=false)
        code = "include($(repr(helper))); ValidationTestRunner.run_tests($(repr(root)); filter=(_ -> $selected), list=$list)"
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
        body_marker = joinpath(root, "body-executed")
        write(active, """
            @testmodule ChosenSetup begin
                const answer = :current
            end
            @testitem "included contract" setup=[ChosenSetup] begin
                write($(repr(body_marker)), "executed")
                @test ChosenSetup.answer === :current
            end
            """)
        ok, output = invoke(root; list=true)
        @test ok
        @test occursin("Listed 1 maintained test items in 1 files; no test bodies executed", output)
        @test !isfile(body_marker)
        @test !isfile(marker)

        ok, output = invoke(root)
        @test ok
        @test occursin("Selected 1 maintained test items", output)
        @test !isfile(marker)
        @test isfile(body_marker)
        @test occursin("Starting [", output)
        @test occursin("run completed", output)
        @test occursin("s elapsed", output)

        write(active, "@testitem \"failed contract\" begin\n@test false\nend\n")
        ok, output = invoke(root)
        @test !ok
        @test occursin("Selected 1 maintained test items in 1 files; run failed", output)

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
