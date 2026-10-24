@testitem "Extensions / PSCAD / unloaded backend is absent" tags=[:extension, :core_only, :importexport] begin
    @test Base.get_extension(LineCableModels, :LineCableModelsPSCADExt) === nothing
    @test_throws ArgumentError("unknown backend formulation :pscad. A backend held in a " *
        "package extension loads with its trigger package.") Formulation(:pscad)
    # The static project works without the extension.
    design = build(CableDesign, "wire", terminal(:core,
        core(Material(:conductor, 1.72e-8, 1.0); r = 0.004),
        insulation(Material(:insulator, 1e14, 2.3); t = 0.002)))
    system = build(LineCableSystem, [design], [Pose2(0, -1)];
        connections = [Dict(:core => 1)], system_id = "absent")
    mktempdir() do directory
        path = export_data(:pscad, system, homogeneous(rho = 100.0);
            file_name = joinpath(directory, "project.pscx"))
        @test path == joinpath(directory, "absent_project.pscx")
        @test occursin("<project name=\"absent\"", read(path, String))
    end
end
