@testmodule TestFixtures begin
    using LineCableModels
    using LineCableModels.DataModel: CableDesign, LineCableSystem, trefoil_formation
    using LineCableModels.EarthProps: EarthModel

    const FIXTURE_ROOT = normpath(joinpath(@__DIR__, "..", "fixtures"))
    const MV_CABLE_DESIGN_PATH = joinpath(FIXTURE_ROOT, "data", "mv_cable_design.json")
    const GOLDEN_ROOT = joinpath(FIXTURE_ROOT, "golden")

    function mv_cable_design()
        aluminum = LineCableModels.Material(
            kind = :conductor,
            rho = 2.8264e-8,
            eps_r = 1.0,
            mu_r = 1.000022,
            T0 = 20.0,
            alpha = 0.00429
        )
        copper = LineCableModels.Material(
            kind = :conductor,
            rho = 1.7241e-8,
            eps_r = 1.0,
            mu_r = 0.999994,
            T0 = 20.0,
            alpha = 0.00393
        )
        semicon_outer = LineCableModels.Material(
            kind = :semicon, rho = 5300.0, eps_r = 32.3
        )
        semicon_inner = LineCableModels.Material(
            kind = :semicon, rho = 1000.0, eps_r = 1000.0
        )
        semicon_screen = LineCableModels.Material(
            kind = :semicon, rho = 500.0, eps_r = 1000.0
        )
        xlpe = LineCableModels.Material(
            kind = :insulator, rho = 197e12, eps_r = 2.3
        )

        parts = LineCableModels.AbstractCablePart[]
        inner_radii = (0.0, 0.00235, 0.00705, 0.01175, 0.01645)
        wire_counts = (1, 6, 12, 18, 24)
        lay_ratios = (0.0, 15.0, 13.5, 12.5, 11.0)
        for index in eachindex(wire_counts)
            wire = LineCableModels.Region(
                Symbol(:core_wire_, index),
                LineCableModels.DiskDefinition(0.00235),
                aluminum
            )
            path = iszero(lay_ratios[index]) ? nothing :
                   LineCableModels.Helix(
                LineCableModels.LayRatio(lay_ratios[index])
            )
            push!(parts,
                LineCableModels.Group(
                    :core,
                    wire;
                    pattern = LineCableModels.Ring(
                        wire_counts[index];
                        r = wire_counts[index] == 1 ? 0.0 :
                            inner_radii[index] + 0.00235
                    ),
                    path
                ))
        end
        for (tag, inner, outer, material) in (
            (:core_semicon_1, 0.02115, 0.02145, semicon_outer),
            (:core_semicon_2, 0.02145, 0.02205, semicon_inner),
            (:core_insulation, 0.02205, 0.03005, xlpe),
            (:core_semicon_3, 0.03005, 0.030350000000000002, semicon_screen),
            (:core_semicon_4, 0.030350000000000002,
            0.030650000000000004, semicon_outer)
        )
            push!(parts, LineCableModels.Region(
                tag,
                LineCableModels.AnnulusDefinition(inner, outer),
                material
            ))
        end

        screen_wire = LineCableModels.Region(
            :sheath_wire,
            LineCableModels.DiskDefinition(0.000475),
            copper
        )
        push!(parts,
            LineCableModels.Group(
                :sheath,
                screen_wire;
                pattern = LineCableModels.Ring(
                    49;
                    r = 0.030650000000000004 + 0.000475
                ),
                path = LineCableModels.Helix(LineCableModels.LayRatio(10.0))
            ))
        tape_inner = 0.0316
        tape_outer = 0.031700000000000006
        tape_width = 0.01
        push!(parts,
            LineCableModels.Group(
                :sheath,
                LineCableModels.Region(
                    :sheath_tape,
                    LineCableModels.SectorDefinition(
                        tape_inner,
                        tape_outer,
                        -tape_width / (tape_inner + tape_outer),
                        2tape_width / (tape_inner + tape_outer)
                    ),
                    copper
                );
                path = LineCableModels.Helix(LineCableModels.LayRatio(10.0))
            ))
        push!(parts,
            LineCableModels.Region(
                :sheath_semicon,
                LineCableModels.AnnulusDefinition(
                    0.031700000000000006,
                    0.03200000000000001
                ),
                semicon_outer
            ))
        push!(parts,
            LineCableModels.Group(
                :jacket,
                LineCableModels.Region(
                    :jacket_metal,
                    LineCableModels.AnnulusDefinition(
                        0.03200000000000001,
                        0.032150000000000005
                    ),
                    aluminum
                )
            ))
        push!(parts,
            LineCableModels.Region(
                :jacket_bedding,
                LineCableModels.AnnulusDefinition(
                    0.032150000000000005,
                    0.032200000000000006
                ),
                xlpe
            ))
        push!(parts,
            LineCableModels.Region(
                :jacket_insulation,
                LineCableModels.AnnulusDefinition(
                    0.032200000000000006,
                    0.034600000000000006
                ),
                xlpe
            ))

        return build(
            CableDesign,
            "test_cable",
            LineCableModels.Stack(parts)
        )
    end

    materials_library() = LineCableModels.Materials.MaterialsLibrary(add_defaults = true)
    copper_material() = LineCableModels.Materials.Material(
        :conductor,
        1.7241e-8,
        1.0,
        1.0,
        20.0,
        0.00393
    )
    aluminum_material() = LineCableModels.Materials.Material(
        :conductor,
        2.8264e-8,
        1.0,
        1.0,
        20.0,
        0.00429
    )
    insulator_material() = LineCableModels.Materials.Material(
        :insulator, 1e14, 2.3, 1.0, 20.0, 0.0
    )
    semicon_material() = LineCableModels.Materials.Material(
        :semicon, 1000.0, 1000.0, 1.0, 20.0, 0.0
    )

    function three_phase_system(; line_length = 1_000.0, spacing = 0.035)
        design = mv_cable_design()
        xa, ya, xb, yb, xc, yc = trefoil_formation(0.0, -1.0, spacing)
        connections(phase) = Dict("core" => phase, "sheath" => 0, "jacket" => 0)
        return build(
            LineCableSystem,
            CableDesign[design, design, design],
            [(xa, ya), (xb, yb), (xc, yc)];
            system_id = "three-phase-system",
            line_length,
            connections = [
                connections(1),
                connections(2),
                connections(3)
            ]
        )
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
