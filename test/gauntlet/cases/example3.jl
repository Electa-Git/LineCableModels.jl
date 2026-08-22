@testitem "PSCAD gauntlet / example 3" tags=[:gauntlet] setup=[GauntletSupport] begin
    using Test
    using DataFrames
    using LineCableModels
    using LineCableModels.Engine
    using .GauntletSupport

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
    design=only(CableBuilder(
        "525kV_1600mm2",
        core_parts,
        sheath_parts,
        armor_parts;
        nominal = nominal_data
    ))

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
    problem=only(SystemBuilder(
        "525kV_1600mm2_bipole",
        design,
        positions;
        length = 1000.0,
        temperature = 20.0,
        earth,
        frequencies = frequencies_value
    ))
    formulation=Formulation(
        earth_impedance = EarthImpedance.Pollaczek(),
        earth_admittance = EarthAdmittance.Pollaczek(),
        insulation_admittance = InsulationAdmittance.Lossless()
    )
    tolerances=(
        local_vs_pscad = (
            Z = (absolute = 1.0e-6, relative = 5.0e-2),
            Y = (absolute = 1.0e-9, relative = 5.0e-2)
        ),
        legacy_vs_current = (
            Z = (absolute = 1.0e-9, relative = 1.0e-4),
            Y = (absolute = 1.0e-12, relative = 1.0e-4)
        )
    )
    case=GauntletCase(
        :example3,
        @__FILE__,
        problem,
        formulation,
        :underground_wedepohl_pollaczek_lossless,
        ["positive:core", "negative:core"],
        true,
        (2, 2, length(frequencies_value)),
        tolerances,
        true
    )

    if !ISHEADLESS
        display(DataFrame(design, :components))
        display(DataFrame(problem.system))
        display(DataFrame(problem.earth_props))
    end

    outcome=run_case(case)
    @test outcome.reference isa LineParameters
    @test outcome.candidate isa LineParameters
    @test size(Z(outcome.reference)) == (2, 2, 101)
    @test size(Y(outcome.reference)) == (2, 2, 101)
    @test frequencies(outcome.reference) == frequencies_value
    @test outcome.comparison.Z.absolute <= tolerances.local_vs_pscad.Z.absolute ||
          outcome.comparison.Z.relative <= tolerances.local_vs_pscad.Z.relative
    @test outcome.comparison.Y.absolute <= tolerances.local_vs_pscad.Y.absolute ||
          outcome.comparison.Y.relative <= tolerances.local_vs_pscad.Y.relative

    if outcome.mode!==:snapshot
        @test outcome.legacy.reference isa LineParameters
        @test outcome.legacy.comparison.Z.absolute <=
              tolerances.legacy_vs_current.Z.absolute ||
              outcome.legacy.comparison.Z.relative <=
              tolerances.legacy_vs_current.Z.relative
        @test outcome.legacy.comparison.Y.absolute <=
              tolerances.legacy_vs_current.Y.absolute ||
              outcome.legacy.comparison.Y.relative <=
              tolerances.legacy_vs_current.Y.relative
        @test outcome.pscad.exit_code == 0
        @test isfile(outcome.pscad.console_path)
        @test outcome.timings.pscad >= 0
        @test outcome.timings.julia.samples == 10
    end
end
