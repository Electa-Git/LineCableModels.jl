@testitem "PSCAD benchmark / 525 kV 1600 mm² bipole" tags=[:gauntlet] setup=[GauntletSupport] begin
    using Test
    using DataFrames
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport
    using .GauntletSupport.PSCADBenchmarks

    materials=MaterialsLibrary(add_defaults = true)
    copper=Material(materials, :copper)
    semicon1=Material(materials, :semicon1)
    semicon2=Material(materials, :semicon2)
    pe=Material(materials, :pe)
    polyacrylate=Material(materials, :polyacrylate)
    lead=Material(materials, :lead)
    pp=Material(materials, :pp)
    steel=Material(materials, :steel)

    num_core_wires=127
    num_armor_wires=68
    strand_diameter=3.6649e-3
    inner_semicon_thickness=2.0e-3
    insulation_thickness=26.0e-3
    outer_semicon_thickness=1.8e-3
    water_blocking_thickness=0.3e-3
    lead_screen_thickness=3.3e-3
    inner_sheath_thickness=3.0e-3
    bedding_thickness=3.0e-3
    armor_wire_diameter=5.827e-3
    jacket_thickness=10.0e-3

    strands_per_layer=6
    strand_layers=6
    @test 1 + strands_per_layer * sum(1:strand_layers) == num_core_wires
    core_parts=(
        Conductor.Stranded(
            :core;
            layers = strand_layers+1,
            wire_radius = strand_diameter/2,
            num_wires = strands_per_layer,
            lay_ratio = 11.0,
            material = copper
        ),
        Insulator.Semicon(
            :core; thickness = inner_semicon_thickness, material = semicon1
        ),
        Insulator.Tubular(
            :core; thickness = insulation_thickness, material = pe
        ),
        Insulator.Semicon(
            :core; thickness = outer_semicon_thickness, material = semicon2
        ),
        Insulator.Semicon(
            :core; thickness = water_blocking_thickness, material = polyacrylate
        )
    )
    sheath_parts=(
        Conductor.Tubular(
            :sheath; thickness = lead_screen_thickness, material = lead
        ),
        Insulator.Tubular(
            :sheath; thickness = inner_sheath_thickness, material = pe
        ),
        Insulator.Tubular(:sheath; thickness = bedding_thickness, material = pp)
    )
    armor_parts=(
        Conductor.Wires(
            :armor;
            wire_radius = armor_wire_diameter/2,
            num_wires = num_armor_wires,
            lay_ratio = 10.0,
            material = steel
        ),
        Insulator.Tubular(:armor; thickness = jacket_thickness, material = pp)
    )
    nominal_data=(
        designation_code = "(N)2XH(F)RK2Y",
        U0 = 500.0,
        U = 525.0,
        conductor_cross_section = 1600.0,
        screen_cross_section = 1000.0,
        resistance = nothing,
        capacitance = nothing,
        inductance = nothing
    )
    design=only(Gridspace(CableBuilder(
        "525kV_1600mm2",
        core_parts,
        sheath_parts,
        armor_parts;
        nominal = nominal_data
    )))

    earth=Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0)
    frequencies_value=collect(10.0 .^ range(0, stop = 6, length = 101))
    positions=(
        at(
            x = -0.5,
            y = -1.0,
            phases = (:core=>1, :sheath=>2, :armor=>3)
        ),
        at(
            x = 0.5,
            y = -1.0,
            phases = (:core=>4, :sheath=>5, :armor=>6)
        )
    )
    problem=only(Gridspace(SystemBuilder(
        "benchmark_525kV_1600mm2_bipole_pscad",
        design,
        positions;
        length = 1000.0,
        temperature = 20.0,
        earth,
        frequencies = frequencies_value
    )))
    reference_problem=LineParametersProblem(
        problem.system;
        temperature = problem.temperature,
        earth_props = problem.earth_props,
        frequencies = problem.frequencies
    )
    reference_formulation=Formulation(
        :pscad;
        earth_impedance = EarthImpedance.Wedepohl(),
        options = (output_stem = "525kV_bipole",)
    )
    formulation=Formulation(
        earth_impedance = EarthImpedance.Pollaczek(),
        earth_admittance = EarthAdmittance.IdealGround(),
        insulation_admittance = InsulationAdmittance.Lossless(),
        options = (
            kron_reduction = false,
            reduce_bundle = false,
            ideal_transposition = false
        )
    )
    tolerances=(
        reference = (
            Z = (absolute = 1.0e-6, relative = 5.0e-2),
            Y = (absolute = 1.0e-9, relative = 5.0e-2)
        ),
        regression = (
            Z = (absolute = 1.0e-9, relative = 1.0e-6),
            Y = (absolute = 1.0e-12, relative = 1.0e-6)
        ),
        performance = (
            median_time_ratio = 1.20,
            bytes_ratio = 1.05,
            allocations_ratio = 1.05
        )
    )
    case=GauntletCase(
        :benchmark_525kV_1600mm2_bipole_pscad,
        :pscad,
        @__FILE__,
        reference_problem,
        reference_formulation,
        problem,
        formulation,
        [
            "cable:1:core", "cable:1:sheath", "cable:1:armor",
            "cable:2:core", "cable:2:sheath", "cable:2:armor"
        ],
        (6, 6, length(frequencies_value)),
        tolerances
    )

    if !ISHEADLESS
        display(DataFrame(design, :components))
        display(DataFrame(problem.system))
        display(DataFrame(problem.earth_props))
    end

    inferred=@inferred compute(problem, formulation)
    @test size(Z(inferred)) == (6, 6, 101)
    execution_options=(verbosity = (default = 0, PSCAD = 2),)
    outcome=run_case(case; options = execution_options)
    @test outcome.reference isa LineParameters
    @test outcome.candidate isa LineParameters
    @test size(Z(outcome.reference)) == (6, 6, 101)
    @test size(Y(outcome.reference)) == (6, 6, 101)
    @test frequencies(outcome.reference) == frequencies_value

    if !ISHEADLESS
        errors=DataFrame(
            outcome.comparison;
            zero_atol = (
                Z = tolerances.reference.Z.absolute,
                Y = tolerances.reference.Y.absolute
            )
        )
        sort!(errors, [:quantity, :rms_absolute], rev = [false, true])
        display(errors)
    end

    if outcome.mode===:snapshot
        @test comparison_passes(outcome.regression.Z, tolerances.regression.Z)
        @test comparison_passes(outcome.regression.Y, tolerances.regression.Y)
        if outcome.performance.comparable
            @test outcome.performance.passes
        end
    else
        @test outcome.reference_execution.backend === :pscad
        @test outcome.reference_execution.exit_code == 0
        @test isfile(outcome.reference_execution.console_path)
        @test outcome.timings.reference >= 0
        @test outcome.timings.julia.samples == 10
    end
end
