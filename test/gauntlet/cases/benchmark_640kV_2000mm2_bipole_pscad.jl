@testitem "PSCAD benchmark / 640 kV 2000 mm² rectangular-stranded bipole" tags=[:gauntlet] setup=[GauntletSupport] begin
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
    lead=Material(materials, :lead)

    core_strand_radius=6.543707496202785e-3/2
    core_strand_area=π*core_strand_radius^2
    ring_counts=(6, 12, 18, 24)
    semicon_tape_thickness=0.3e-3
    inner_semicon_thickness=0.768e-3
    insulation_thickness=32.0e-3
    outer_semicon_thickness=0.472e-3
    lead_screen_thickness=3.3e-3
    inner_sheath_thickness=3.0e-3
    aluminum_tape_thickness=0.15e-3
    jacket_thickness=6.05e-3

    const DM=LineCableModels.DataModel
    core_conductor=let group=DM.ConductorGroup(DM.Tubular(
            0.0, core_strand_radius, copper
        ))
        for count in ring_counts
            radius_in=group.r_ex
            strand_width=prevfloat(2π*radius_in/count)
            strand_thickness=core_strand_area/strand_width
            radius_ext=sqrt(
                radius_in^2+count*strand_width*strand_thickness/π
            )
            add!(
                group,
                DM.RectStrands(
                    radius_in,
                    radius_ext,
                    strand_thickness,
                    strand_width,
                    count,
                    13.0,
                    copper
                )
            )
        end
        group
    end
    @test core_conductor.cross_section ≈ 61core_strand_area

    r_in=core_conductor.r_ex
    core_insulation=DM.InsulatorGroup(DM.Semicon(
        r_in, r_in+semicon_tape_thickness, polyacrylate
    ))
    r_in=core_insulation.r_ex
    add!(core_insulation, DM.Semicon(r_in, r_in+inner_semicon_thickness, semicon1))
    r_in=core_insulation.r_ex
    add!(core_insulation, DM.Insulator(r_in, r_in+insulation_thickness, xlpe))
    r_in=core_insulation.r_ex
    add!(core_insulation, DM.Semicon(r_in, r_in+outer_semicon_thickness, semicon2))
    r_in=core_insulation.r_ex
    add!(core_insulation, DM.Semicon(
        r_in, r_in+semicon_tape_thickness, polyacrylate
    ))
    core_component=DM.CableComponent("core", core_conductor, core_insulation)

    r_in=core_insulation.r_ex
    sheath_conductor=DM.ConductorGroup(DM.Tubular(
        r_in, r_in+lead_screen_thickness, lead
    ))
    r_in=sheath_conductor.r_ex
    sheath_insulation=DM.InsulatorGroup(DM.Insulator(
        r_in, r_in+inner_sheath_thickness, pe
    ))
    sheath_component=DM.CableComponent("sheath", sheath_conductor, sheath_insulation)

    r_in=sheath_insulation.r_ex
    jacket_conductor=DM.ConductorGroup(DM.Tubular(
        r_in, r_in+aluminum_tape_thickness, aluminum
    ))
    r_in=jacket_conductor.r_ex
    jacket_insulation=DM.InsulatorGroup(DM.Insulator(
        r_in, r_in+jacket_thickness, pe
    ))
    jacket_component=DM.CableComponent("jacket", jacket_conductor, jacket_insulation)

    design=DM.CableDesign(
        "640kV_2000mm2",
        [core_component, sheath_component, jacket_component];
        nominal_data = DM.NominalData()
    )

    earth=Earth(rho = 100.0, eps_r = 10.0, mu_r = 1.0)
    frequencies_value=collect(10.0 .^ range(0, stop = 6, length = 101))
    positions=(
        at(
            x = -0.5,
            y = -1.0,
            phases = (:core=>1, :sheath=>2, :jacket=>3)
        ),
        at(
            x = 0.5,
            y = -1.0,
            phases = (:core=>4, :sheath=>5, :jacket=>6)
        )
    )
    problem=SystemBuilder(
        "benchmark_640kV_2000mm2_bipole_pscad",
        design,
        positions;
        length = 1000.0,
        temperature = 20.0,
        earth,
        frequencies = frequencies_value
    )
    reference_problem=LineParametersProblem(
        problem.system;
        temperature = problem.temperature,
        earth_props = problem.earth_props,
        frequencies = problem.frequencies
    )
    reference_formulation=Formulation(
        :pscad;
        earth_impedance = EarthImpedance.Wedepohl()
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
        :benchmark_640kV_2000mm2_bipole_pscad,
        :pscad,
        @__FILE__,
        reference_problem,
        reference_formulation,
        problem,
        formulation,
        [
            "cable:1:core", "cable:1:sheath", "cable:1:jacket",
            "cable:2:core", "cable:2:sheath", "cable:2:jacket"
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
    execution_options=(reference = (
        output_stem = "640kV_2000_bipole",
        verbosity = (default = 0, PSCAD = 2)
    ),)
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
