@testitem "PSCAD / reader failures propagate and a failed diagnostics read keeps the failure" tags=[:integration, :pscad] setup=[PSCADStandIn] begin
    const P = Base.get_extension(LineCableModels, :LineCableModelsPSCADExt)
    const S = PSCADStandIn
    frequencies = 10.0 .^ range(-1, 5; length = 101)
    function staged(root)
        design = build(CableDesign, "wire", terminal(:core,
            core(Material(:conductor, 1.72e-8, 1.0); r = 0.004),
            insulation(Material(:insulator, 1e14, 2.3); t = 0.002)))
        system = build(LineCableSystem, [design], [Pose2(0, -1)];
            connections = [Dict(:core => 1)], system_id = "stand-in")
        problem = LineParametersProblem(system; earth_props = homogeneous(rho = 100.0), frequencies)
        project = joinpath(mkpath(joinpath(root, "case")), "stand_in.pscx")
        write(project, P._pscad_inputs(problem, Formulation(:pscad),
            LineCableModels.DataModel.flatten.(system.designs)).project)
        return project
    end
    settings = Dict(row => Dict{String, Any}() for row in ("ground", "frequency", "configuration"))
    function failure(root; timeout = 60)
        project = staged(root)
        return try
            P.run_remote_pscad(P.RemoteConfig(S.remote(root; timeout)), project,
                joinpath(dirname(project), "outputs"), settings, frequencies; line_name = "lcm_read")
            nothing
        catch error
            error
        end
    end
    # A detailed file that never completes fails the run, and the run's error states the phase
    # and the reader's failure.
    mktempdir() do root
        application = S.install(root)
        application.solves = false
        caught = failure(root; timeout = 1)
        @test caught isa ErrorException
        message = sprint(showerror, caught)
        @test occursin("during detailed-output completion", message)
        @test occursin("did not become complete within 1 seconds; last observed 0 of 101 rows. " *
                       "PSCAD writes the detailed files into its working directory, which must " *
                       "be the configuration's `shared_root`. A borrowed PSCAD must have been " *
                       "started there.", message)
        @test isempty(S.loaded(application))
        @test !S.silence(application)
    end
    # A failed read of PSCAD's diagnostics takes their place, and the failure stays the error.
    mktempdir() do root
        application = S.install(root)
        application.compile_error = "synthetic compile failure"
        application.output_error = "synthetic output failure"
        caught = failure(root)
        @test caught isa ErrorException
        message = sprint(showerror, caught)
        @test occursin("during line-constants computation: Python: RuntimeError: " *
                       "synthetic compile failure", message)
        @test occursin("PSCAD diagnostics are unavailable", message)
        @test occursin("synthetic output failure", message)
        console = read(joinpath(root, "case", "outputs", "pscad-console.txt"), String)
        @test occursin("synthetic compile failure", console)
        @test occursin("synthetic output failure", console)
        @test isempty(S.loaded(application))
        @test !S.silence(application)
    end
end
