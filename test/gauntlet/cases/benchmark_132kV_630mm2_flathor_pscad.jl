@testitem "PSCAD benchmark / 132 kV 630 mm² flat horizontal" tags=[:gauntlet] setup=[GauntletSupport] begin
    using Test
    using DataFrames
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport
    using .GauntletSupport.PSCADBenchmarks

    materials=MaterialsLibrary(add_defaults = true)
    aluminum=Material(materials, :aluminum)
    copper=Material(materials, :copper)
    pe=Material(materials, :pe)
    xlpe=Material(materials, :xlpe)
    semicon1=Material(materials, :semicon1)
    semicon2=Material(materials, :semicon2)
    polyacrylate=Material(materials, :polyacrylate)

    core_layers=5
    core_strand_diameter=3.6648500338805653e-3
    screen_wires=19
    screen_wire_diameter=2.5881867280128637e-3
    semicon_tape_thickness=0.3e-3
    inner_semicon_thickness=0.768e-3
    insulation_thickness=16.0e-3
    outer_semicon_thickness=0.472e-3
    copper_tape_thickness=0.1e-3
    copper_tape_width=10.0e-3
    water_blocking_thickness=0.94e-3
    aluminum_tape_thickness=0.15e-3
    jacket_thickness=5.05e-3

    core_parts=(
        Conductor.Stranded(
            :core;
            layers = core_layers,
            wire_radius = core_strand_diameter/2,
            num_wires = 6,
            lay_ratio = 13.0,
            material = copper
        ),
        Insulator.Semicon(
            :core; thickness = semicon_tape_thickness, material = polyacrylate
        ),
        Insulator.Semicon(
            :core; thickness = inner_semicon_thickness, material = semicon1
        ),
        Insulator.Tubular(
            :core; thickness = insulation_thickness, material = xlpe
        ),
        Insulator.Semicon(
            :core; thickness = outer_semicon_thickness, material = semicon2
        ),
        Insulator.Semicon(
            :core; thickness = semicon_tape_thickness, material = polyacrylate
        )
    )
    sheath_parts=(
        Conductor.Wires(
            :sheath;
            wire_radius = screen_wire_diameter/2,
            num_wires = screen_wires,
            lay_ratio = 10.0,
            material = copper
        ),
        Conductor.Strip(
            :sheath;
            thickness = copper_tape_thickness,
            width = copper_tape_width,
            lay_ratio = 10.0,
            material = copper
        ),
        Insulator.Semicon(
            :sheath; thickness = water_blocking_thickness, material = polyacrylate
        )
    )
    jacket_parts=(
        Conductor.Tubular(
            :jacket; thickness = aluminum_tape_thickness, material = aluminum
        ),
        Insulator.Tubular(:jacket; thickness = jacket_thickness, material = pe)
    )
    nominal_data=(
        designation_code = "A2XS(FL)2Y",
        U0 = 76.0,
        U = 132.0,
        conductor_cross_section = 630.0,
        screen_cross_section = 95.0,
        resistance = 0.0395,
        capacitance = 0.2,
        inductance = 0.12
    )
    design=only(CableBuilder(
        "132kV_630mm2",
        core_parts,
        sheath_parts,
        jacket_parts;
        nominal = nominal_data
    ))

    earth=Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0)
    frequencies_value=collect(10.0 .^ range(0, stop = 6, length = 101))
    positions=(
        at(
            x = -0.5,
            y = -1.0,
            phases = (:core=>1, :sheath=>2, :jacket=>3)
        ),
        at(
            x = 0.0,
            y = -1.0,
            phases = (:core=>4, :sheath=>5, :jacket=>6)
        ),
        at(
            x = 0.5,
            y = -1.0,
            phases = (:core=>7, :sheath=>8, :jacket=>9)
        )
    )
    problem=only(SystemBuilder(
        "benchmark_132kV_630mm2_flathor_pscad",
        design,
        positions;
        length = 1000.0,
        temperature = 20.0,
        earth,
        frequencies = frequencies_value
    ))
    reference_problem=LineParametersProblem(
        problem.system;
        temperature = problem.temperature,
        earth_props = problem.earth_props,
        frequencies = problem.frequencies
    )
    reference_formulation=Formulation(
        :pscad;
        earth_impedance = EarthImpedance.Wedepohl(),
        options = (output_stem = "132kV_630_flat",)
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
        :benchmark_132kV_630mm2_flathor_pscad,
        :pscad,
        @__FILE__,
        reference_problem,
        reference_formulation,
        problem,
        formulation,
        [
            "cable:1:core", "cable:1:sheath", "cable:1:jacket",
            "cable:2:core", "cable:2:sheath", "cable:2:jacket",
            "cable:3:core", "cable:3:sheath", "cable:3:jacket"
        ],
        (9, 9, length(frequencies_value)),
        tolerances
    )

    if !ISHEADLESS
        display(DataFrame(design, :components))
        display(DataFrame(problem.system))
        display(DataFrame(problem.earth_props))
    end

    inferred=@inferred compute(problem, formulation)
    @test size(Z(inferred)) == (9, 9, 101)
    execution_options=(verbosity = (default = 0, PSCAD = 2),)
    outcome=run_case(case; options = execution_options)
    @test outcome.reference isa LineParameters
    @test outcome.candidate isa LineParameters
    @test size(Z(outcome.reference)) == (9, 9, 101)
    @test size(Y(outcome.reference)) == (9, 9, 101)
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
