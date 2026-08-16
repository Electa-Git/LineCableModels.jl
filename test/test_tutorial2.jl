@testitem "examples/tutorial2.jl public workflow" setup = [defaults] begin
    const LCM = LineCableModels
    const PB = LineCableModels.PlotBuilder

    materials = LCM.MaterialsLibrary(add_defaults=true)
    aluminum = LCM.Material(materials, :aluminum)
    copper = LCM.Material(materials, :copper)
    polyacrylate = LCM.Material(materials, :polyacrylate)
    semicon1 = LCM.Material(materials, :semicon1)
    semicon2 = LCM.Material(materials, :semicon2)
    pe = LCM.Material(materials, :pe)

    num_co_wires = 61
    num_sc_wires = 49
    d_w = 4.7e-3
    t_sc_in = 0.6e-3
    t_ins = 8e-3
    t_sc_out = 0.3e-3
    d_ws = 0.95e-3
    t_cut = 0.1e-3
    w_cut = 10e-3
    t_wbt = 0.3e-3
    t_sct = 0.3e-3
    t_alt = 0.15e-3
    t_pet = 0.05e-3
    t_jac = 2.4e-3

    core_wire_counts = (1, 6, 12, 18, 24)
    @test sum(core_wire_counts) == num_co_wires
    core_conductors = (
        LCM.Conductor.Wires(
            :core; wire_radius=d_w / 2, num_wires=1,
            lay_ratio=0.0, material=aluminum),
        LCM.Conductor.Wires(
            :core; wire_radius=d_w / 2, num_wires=6,
            lay_ratio=15.0, material=aluminum),
        LCM.Conductor.Wires(
            :core; wire_radius=d_w / 2, num_wires=12,
            lay_ratio=13.5, material=aluminum),
        LCM.Conductor.Wires(
            :core; wire_radius=d_w / 2, num_wires=18,
            lay_ratio=12.5, material=aluminum),
        LCM.Conductor.Wires(
            :core; wire_radius=d_w / 2, num_wires=24,
            lay_ratio=11.0, material=aluminum),
    )
    core_parts = (
        core_conductors,
        (
            LCM.Insulator.Semicon(
                :core; thickness=t_sct, material=polyacrylate),
            LCM.Insulator.Semicon(
                :core; thickness=t_sc_in, material=semicon1),
        ),
        LCM.Insulator.Tubular(:core; thickness=t_ins, material=pe),
        (
            LCM.Insulator.Semicon(
                :core; thickness=t_sc_out, material=semicon2),
            LCM.Insulator.Semicon(
                :core; thickness=t_sct, material=polyacrylate),
        ),
    )
    sheath_parts = (
        LCM.Conductor.Wires(
            :sheath;
            wire_radius=d_ws / 2,
            num_wires=num_sc_wires,
            lay_ratio=10.0,
            material=copper,
        ),
        LCM.Conductor.Strip(
            :sheath;
            thickness=t_cut,
            width=w_cut,
            lay_ratio=10.0,
            material=copper,
        ),
        LCM.Insulator.Semicon(
            :sheath; thickness=t_wbt, material=polyacrylate),
    )
    jacket_parts = (
        LCM.Conductor.Tubular(:jacket; thickness=t_alt, material=aluminum),
        LCM.Insulator.Tubular(:jacket; thickness=t_pet, material=pe),
        LCM.Insulator.Tubular(:jacket; thickness=t_jac, material=pe),
    )

    cable_id = "18kV_1000mm2"
    datasheet_info = (;
        designation_code="NA2XS(FL)2Y",
        U0=18.0,
        U=30.0,
        conductor_cross_section=1000.0,
        screen_cross_section=35.0,
        resistance=0.0291,
        capacitance=0.39,
        inductance=0.3,
    )

    core_design = only(LCM.CableBuilder(cable_id, core_parts; nominal=datasheet_info))
    screened_design = only(LCM.CableBuilder(
        cable_id, core_parts, sheath_parts; nominal=datasheet_info))
    cable_design = only(LCM.CableBuilder(
        cable_id, core_parts, sheath_parts, jacket_parts; nominal=datasheet_info))

    @test length(core_design.components) == 1
    @test length(screened_design.components) == 2
    @test length(cable_design.components) == 3
    @test [component.id for component in cable_design.components] ==
          ["core", "sheath", "jacket"]
    @test length(cable_design.components[1].conductor_group.layers) == 5
    @test length(cable_design.components[1].insulator_group.layers) == 5
    @test cable_design.nominal_data.designation_code == "NA2XS(FL)2Y"

    constants = LCM.compute!(LCM.CableConstantsProblem(cable_design), LCM.Formulation())
    @test isapprox(constants.R * 1e3, 0.0275677; atol=1e-5)
    @test isapprox(constants.L * 1e6, 0.287184; atol=1e-5)
    @test isapprox(constants.C * 1e9, 0.413357; atol=1e-5)
    @test DataFrame(constants) isa DataFrame
    @test DataFrame(cable_design, :components) isa DataFrame
    @test DataFrame(cable_design, :detailed) isa DataFrame

    cable_render = PB.make_render(
        LCM.DataModel.CablePreviewPlotSpec, cable_design)
    @test length(cable_render.figures) == 1

    mktempdir() do temp_dir
        library = LCM.CablesLibrary()
        add!(library, cable_design)
        library_file = joinpath(temp_dir, "cables_library.json")
        @test LCM.save(library; file_name=library_file) == abspath(library_file)

        loaded_library = LCM.CablesLibrary()
        LCM.load!(loaded_library; file_name=library_file)
        loaded_design = get(loaded_library, cable_id)
        @test loaded_design.cable_id == cable_id
        @test length(loaded_design.components) == 3
        @test DataFrame(loaded_library) isa DataFrame

        frequencies = collect(10.0 .^ range(0, stop=6, length=10))
        formation = LCM.trifoil(
            x=0.0,
            y=-1.0,
            spacing=70e-3,
            phases=(
                :core => (1, 2, 3),
                :sheath => 0,
                :jacket => 0,
            ),
        )
        problem = only(LCM.SystemBuilder(
            "18kV_1000mm2_trifoil",
            loaded_design,
            formation;
            length=1000.0,
            temperature=20.0,
            earth=LCM.Earth(rho=100.0, eps_r=10.0, mu_r=1.0),
            frequencies,
        ))
        system = problem.system
        earth_model = problem.earth_props

        @test length(system.cables) == 3
        @test DataFrame(system) isa DataFrame
        @test DataFrame(earth_model) isa DataFrame
        positions = [(cable.horz, cable.vert) for cable in system.cables]
        @test hypot(
            positions[1][1] - positions[2][1],
            positions[1][2] - positions[2][2],
        ) ≈ 70e-3

        system_render = PB.make_render(
            LCM.DataModel.SystemPreviewPlotSpec,
            system;
            earth_model,
            zoom_factor=2.0,
        )
        @test length(system_render.figures) == 1

        pscad = LCM.export_data(
            :pscad,
            system,
            earth_model;
            file_name=joinpath(temp_dir, "pscad_export.pscx"),
        )
        atp = LCM.export_data(
            :atp,
            system,
            earth_model;
            file_name=joinpath(temp_dir, "atp_export.xml"),
        )
        @test isfile(pscad)
        @test filesize(pscad) > 0
        @test isfile(atp)
        @test filesize(atp) > 0

        if get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false") == "true"
            cable_plot = LCM.preview(
                cable_design;
                display_plot=false,
                controls=false,
            )
            @test isempty(cable_plot.controls)
            system_plot = LCM.preview(
                system;
                earth_model,
                zoom_factor=2.0,
                display_plot=false,
                controls=false,
            )
            @test isempty(system_plot.controls)
        end
    end
end
