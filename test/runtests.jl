using LineCableModels
using Test
using TestItemRunner

@testsnippet defaults begin
    const TEST_TOL = 1e-8
    using Measurements
    using Measurements: measurement, uncertainty, value
    using DataFrames
end

@testsnippet defs_materials begin
    materials = MaterialsLibrary(add_defaults=true)
    copper_props = Material(1.7241e-8, 1.0, 1.0, 20.0, 0.00393)
    aluminum_props = Material(2.8264e-8, 1.0, 1.0, 20.0, 0.00429)
    insulator_props = Material(1e14, 2.3, 1.0, 20.0, 0.0)
    semicon_props = Material(1000.0, 1000.0, 1.0, 20.0, 0.0)
end

@testsnippet cable_system_export begin

    cables_library = CablesLibrary()
    cables_library = load!(cables_library, file_name=joinpath(@__DIR__, "./cable_test.json"))

    # Retrieve the reloaded design
    cable_design = collect(values(cables_library.data))[1]
    x0, y0 = 0.0, -1.0
    xa, ya, xb, yb, xc, yc = trifoil_formation(x0, y0, 0.035);

    # Initialize the `LineCableSystem` with the first cable (phase A):
    cablepos = CablePosition(cable_design, xa, ya,
        Dict("core" => 1, "sheath" => 0, "jacket" => 0))
    cable_system = LineCableSystem("test_cable_sys", 1000.0, cablepos)

    # Add remaining cables (phases B and C):
    add!(cable_system, cable_design, xb, yb,
        Dict("core" => 2, "sheath" => 0, "jacket" => 0))
    add!(cable_system, cable_design, xc, yc,
        Dict("core" => 3, "sheath" => 0, "jacket" => 0))

    freqs = sort(abs.(randn(3)))
    earth_params_atp = EarthModel(freqs, 100.0, 10.0, 1.0)
    num_phases = cable_system.num_phases

    # Create minimal mock objects for the other required arguments
    problem_atp = LineParametersProblem(
        cable_system,
        temperature=20.0,  # Operating temperature
        earth_props=earth_params_atp,
        frequencies=freqs,   # Frequency for the analysis
        );

end

@run_package_tests verbose = true
