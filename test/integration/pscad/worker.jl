# Run in ext/LineCableModelsPSCADExt/remote, the worker's own PythonCall environment.
using Test, TOML, PythonCall

module Worker
include(joinpath(@__DIR__, "../../../ext/LineCableModelsPSCADExt/remote/runner.jl"))
end

const automation = pyimport("runpy").run_path(joinpath(@__DIR__, "worker_fixture.py"))

function run_fixture(check; configure = state -> nothing, input_change = input -> nothing,
        arguments_change = arguments -> nothing)
    mktempdir() do root
        output = joinpath(root, "outputs")
        project = joinpath(root, "generated.pscx")
        write(project, "synthetic project")
        input = Dict("schema_version" => 3,
            "solver" => Dict("version" => "5.1.0", "installation" => "fixture"),
            "native_settings" => Dict(
                "ground" => Dict("EarthForm" => Dict("value" => 2,
                    "readback" => "DIRECT_NUMERICAL_INTEGRATION")),
                "frequency" => Dict("Output" => Dict("value" => 1, "readback" => "YES")),
                "configuration" => Dict("Freq" => Dict("value" => 60.0, "readback" => 60.0))))
        input_change(input)
        open(io -> TOML.print(io, input), joinpath(root, "computation.toml"), "w")
        state = automation["install"](joinpath(root, "raw"))
        try
            configure(state)
            arguments = [project, output, "generated", "fixture", "default",
                "0.1", "10.0", "1", "5.1.0", "0"]
            arguments_change(arguments)
            check(state, output, arguments)
        finally
            automation["uninstall"](state)
        end
    end
end

@testset "PSCAD worker / native automation protocol" begin
    @test_throws ArgumentError Worker.main(String[])
    for (index, value) in ((10, "3"), (9, "5.0.0"), (4, "../escape"))
        run_fixture(; arguments_change = args -> (args[index] = value)) do state, output, args
            @test_throws ArgumentError Worker.main(args)
            @test pyconvert(Int, state.compiled) == 0
            @test isempty(pyconvert(Vector{String}, state.loaded))
        end
    end
    run_fixture(; input_change = input -> (input["schema_version"] = 2)) do state, output, args
        @test_throws ArgumentError Worker.main(args)
        @test pyconvert(Int, state.compiled) == 0
    end

    # Complete worker execution, including the fallback canvas lookup and cleanup
    # diagnostics. Readbacks and raw output bytes are checked independently.
    run_fixture(; configure = state -> begin
        state.canvas_fallback = true
        state.cleanup_error = true
    end, arguments_change = args -> (args[10] = "2")) do state, output, args
        @test Worker.main(args) === nothing
        settings = TOML.parsefile(joinpath(output, "native-settings.toml"))
        @test settings["configuration"]["Freq"] == 60.0
        @test settings["ground"]["EarthForm"] == "DIRECT_NUMERICAL_INTEGRATION"
        @test settings["frequency"]["Output"] == "YES"
        @test pyconvert(Dict, state.frequency.parameters())["Numf"] == 1
        @test pyconvert(Dict, state.line.parameters())["Name"] == "fixture"
        for name in ("zm", "zp", "ym", "yp")
            @test read(joinpath(output, "result_" * name * ".out"), String) ==
                "LOG10(FN) FN element\n-1 0.1 2\n1 10 3\n"
        end
        @test TOML.parsefile(joinpath(output, "solver.toml")) ==
            Dict("version" => "5.1.0", "installation" => "fixture")
        @test parse(Float64, read(joinpath(output, "timing.txt"), String)) >= 0
        @test pyconvert(Int, state.saved) == 1
        @test pyconvert(Int, state.compiled) == 1
        @test pyconvert(Int, state.unloaded) == 1
        @test pyconvert(Int, state.quit) == 1
        log = read(joinpath(output, "pscad-console.txt"), String)
        @test occursin("synthetic unload failure", log)
        @test occursin("synthetic quit failure", log)
        @test occursin("automation fixture output", log)
    end

    for (configure, expected, compiled, unloaded, quit) in (
        (s -> (s.automation_version = "0.0.0"), "mhi.pscad 3.1.2 is required", 0, 0, 0),
        (s -> (s.licensed = false), "refused the configured license", 0, 0, 1),
        (s -> (s.version = "5.0.0"), "unexpected version", 0, 0, 1),
        (s -> (s.change_identity_after = 0), "changed after input preparation", 0, 0, 1),
        (s -> (s.line_count = 0), "exactly one line-data row", 0, 1, 1),
        (s -> s.components.remove(s.ground), "0 components named", 0, 1, 1),
        (s -> s.frequency.values.pop("FS"), "has no field FS", 0, 1, 1),
        (s -> (s.reject_field = "Freq"), "expected readback 60.0", 0, 1, 1),
        (s -> s.components.remove(s.cable), "no master:Cable_Coax", 0, 1, 1),
        (s -> s.cable.values.__setitem__("elim2", "ELIMINATE"), "conductor elimination", 0, 1, 1),
        (s -> begin s.compile_error = true; s.diagnostics_error = true end,
            "synthetic compile failure", 1, 1, 1),
        (s -> (s.change_identity_after = 1), "changed during calculation", 1, 1, 1),
    )
        run_fixture(; configure) do state, output, args
            caught = try
                Worker.main(args)
                nothing
            catch error
                error
            end
            @test caught isa ErrorException
            @test occursin(expected, sprint(showerror, caught))
            @test pyconvert(Int, state.compiled) == compiled
            @test pyconvert(Int, state.unloaded) == unloaded
            @test pyconvert(Int, state.quit) == quit
            @test occursin(expected, read(joinpath(output, "pscad-console.txt"), String))
            @test !isfile(joinpath(output, "solver.toml"))
        end
    end
end
