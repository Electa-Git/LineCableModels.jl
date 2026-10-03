@testitem "Core / test discovery rejects excluded setups malformed source and empty selection" tags=[:unit, :units, :slow] begin
    helper = joinpath(pkgdir(LineCableModels), "test", "support", "runner.jl")
    include(helper)
    directory = joinpath(pkgdir(LineCableModels), "test")
    item(file, name, tags) = (; filename=joinpath(directory, file), name, tags)
    ordinary = item("unit/current.jl", "unresolved current item", [:unit])
    sampling = item("integration/sampling.jl", "current sampling item", [:integration])
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
            @testitem "included item" setup=[ChosenSetup] begin
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

        write(active, "@testitem \"failed item\" begin\n@test false\nend\n")
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

@testitem "Core / changed: runs quality in a fresh process" tags=[:unit, :units, :slow] begin
    root = pkgdir(LineCableModels)
    # Before quality runs, the parent process gains a method inside a package module,
    # defined from a test file (this one), as `test/integration/pscad/parser_tests.jl`
    # does. A2 reports such a method when it shares the process; in the fresh process
    # it passes.
    runner = joinpath(root, "test", "support", "runner.jl")
    planted = joinpath(root, "test", "unit", "core", "test_harness.jl")
    program = """
        using TestItemRunner
        include($(repr(runner)))
        import LineCableModels
        include_string(LineCableModels.PSCAD, "planted_probe(::Val{:planted}) = 1", $(repr(planted)))
        method = only(methods(LineCableModels.PSCAD.planted_probe))
        println("PLANTED ", String(method.file))
        ValidationTestRunner.run_changed($(repr(root)),
            Set([("test/unit/units/units.jl", "Units / locked public vocabulary")]);
            quality = ["Quality / architecture / A2 placement"])
        println("CHANGED RUN PASSED")
        """
    output = IOBuffer()
    command = `$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) -e $program`
    process = run(pipeline(ignorestatus(command); stdout = output, stderr = output); wait = false)
    status = timedwait(() -> process_exited(process), 600.0; pollint = 0.1)
    status === :timed_out && kill(process)
    wait(process)
    text = String(take!(output))
    @test status === :ok
    @test success(process)
    @test occursin("PLANTED $planted", text)
    @test count("Selected 1 maintained test items in 1 files; run completed", text) == 2
    @test occursin("Quality items in a fresh process: Quality / architecture / A2 placement", text)
    @test occursin("CHANGED RUN PASSED", text)
    success(process) || println(text)
end
