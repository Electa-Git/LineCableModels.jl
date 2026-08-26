@testmodule TestFixtures begin
    using LineCableModels
    using LineCableModels.DataModel: CablePosition, LineCableSystem, trifoil_formation
    using LineCableModels.EarthProps: EarthModel

    const FIXTURE_ROOT = normpath(joinpath(@__DIR__, "..", "fixtures"))
    const MV_CABLE_DESIGN_PATH = joinpath(FIXTURE_ROOT, "data", "mv_cable_design.json")
    const GOLDEN_ROOT = joinpath(FIXTURE_ROOT, "golden")

    function mv_cable_design()
        library = CablesLibrary()
        load!(library; file_name = MV_CABLE_DESIGN_PATH)
        return only(values(library.data))
    end

    materials_library() = LineCableModels.Materials.MaterialsLibrary(add_defaults = true)
    copper_material() = LineCableModels.Materials.Material(
        1.7241e-8,
        1.0,
        1.0,
        20.0,
        0.00393
    )
    aluminum_material() = LineCableModels.Materials.Material(
        2.8264e-8,
        1.0,
        1.0,
        20.0,
        0.00429
    )
    insulator_material() = LineCableModels.Materials.Material(1e14, 2.3, 1.0, 20.0, 0.0)
    semicon_material() = LineCableModels.Materials.Material(1000.0, 1000.0, 1.0, 20.0, 0.0)

    function three_phase_system(; line_length = 1_000.0, spacing = 0.035)
        design = mv_cable_design()
        xa, ya, xb, yb, xc, yc = trifoil_formation(0.0, -1.0, spacing)
        connections(phase) = Dict("core" => phase, "sheath" => 0, "jacket" => 0)
        position = CablePosition(design, xa, ya, connections(1))
        system = LineCableSystem("three-phase-system", line_length, position)
        add!(system, design, xb, yb, connections(2))
        add!(system, design, xc, yc, connections(3))
        return system
    end

    function line_parameters_problem(; frequencies = [50.0])
        system = three_phase_system()
        earth = EarthModel(100.0, 10.0, 1.0)
        return LineCableModels.Engine.LineParametersProblem(
            system;
            temperature = 20.0,
            earth_props = earth,
            frequencies
        )
    end

    function two_conductor_results(; frequencies = [50.0, 500.0])
        frequency = collect(frequencies)
        count = length(frequency)
        resistance = repeat([2.0e-4 0.5e-4; 0.5e-4 2.5e-4], 1, 1, count)
        inductance = repeat([1.5e-6 0.2e-6; 0.2e-6 1.8e-6], 1, 1, count)
        conductance = repeat([2.0e-9 -0.5e-9; -0.5e-9 2.5e-9], 1, 1, count)
        capacitance = repeat([8.0e-10 -2.0e-10; -2.0e-10 9.0e-10], 1, 1, count)
        angular_frequency = reshape(2π .* frequency, 1, 1, :)
        impedance = complex.(resistance, inductance .* angular_frequency)
        admittance = complex.(conductance, capacitance .* angular_frequency)
        return LineCableModels.Engine.LineParameters(
            impedance,
            admittance,
            frequency;
            basis = :pul
        )
    end

    function cable_monte_carlo_result()
        values = [1.0, 2.0, 3.0, 4.0]
        summary = LineCableModels.UQ.SampleSummary(values)
        histogram = LineCableModels.UQ.HistogramDensity(
            [1.0, 3.0, 5.0],
            [0.25, 0.25]
        )
        representation = LineCableModels.DataModel.CableConstants(2.5, 2.5, 2.5)
        statistics = (R = summary, L = summary, C = summary)
        samples = (
            R = copy(values),
            L = copy(values),
            C = copy(values)
        )
        histograms = (R = histogram, L = histogram, C = histogram)
        formulation = LineCableModels.MonteCarlo(
            LineCableModels.Formulation();
            trials = length(values),
            seed = 1,
            return_samples = true,
            return_histograms = true
        )
        return LineCableModels.MonteCarloResult(
            formulation,
            [representation],
            [statistics],
            [samples],
            [histograms],
            UInt64(1),
            UInt64[1],
            [length(values)]
        )
    end
end

@testsnippet MaterialFixtures begin
    materials = TestFixtures.materials_library()
    copper_props = TestFixtures.copper_material()
    aluminum_props = TestFixtures.aluminum_material()
    insulator_props = TestFixtures.insulator_material()
    semicon_props = TestFixtures.semicon_material()
end

@testsnippet CableSystemFixture begin
    cable_system = TestFixtures.three_phase_system()
    problem_atp = TestFixtures.line_parameters_problem()
    earth_props = problem_atp.earth_props
    freqs = problem_atp.frequencies
    num_phases = LineCableModels.nphases(cable_system)
end
