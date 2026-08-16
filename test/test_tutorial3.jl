@testitem "examples/tutorial3.jl public workflow" setup = [defaults] begin
    const LCM = LineCableModels
    const PB = LineCableModels.PlotBuilder

    materials = LCM.MaterialsLibrary(add_defaults=true)
    copper = LCM.Material(materials, :copper)
    semicon1 = LCM.Material(materials, :semicon1)
    semicon2 = LCM.Material(materials, :semicon2)
    pe = LCM.Material(materials, :pe)
    polyacrylate = LCM.Material(materials, :polyacrylate)
    lead = LCM.Material(materials, :lead)
    pp = LCM.Material(materials, :pp)
    steel = LCM.Material(materials, :steel)

    num_co_wires = 127
    num_ar_wires = 68
    d_w = 3.6649e-3
    t_sc_in = 2e-3
    t_ins = 26e-3
    t_sc_out = 1.8e-3
    t_wbt = 0.3e-3
    t_sc = 3.3e-3
    t_pe = 3e-3
    t_bed = 3e-3
    d_wa = 5.827e-3
    t_jac = 10e-3

    n_strands = 6
    n_layers = 6
    @test 1 + n_strands * sum(1:n_layers) == num_co_wires
    core_parts = (
        LCM.Conductor.Stranded(
            :core;
            layers=n_layers + 1,
            wire_radius=d_w / 2,
            num_wires=n_strands,
            lay_ratio=11.0,
            material=copper,
        ),
        LCM.Insulator.Semicon(:core; thickness=t_sc_in, material=semicon1),
        LCM.Insulator.Tubular(:core; thickness=t_ins, material=pe),
        LCM.Insulator.Semicon(:core; thickness=t_sc_out, material=semicon2),
        LCM.Insulator.Semicon(
            :core; thickness=t_wbt, material=polyacrylate),
    )
    sheath_parts = (
        LCM.Conductor.Tubular(:sheath; thickness=t_sc, material=lead),
        LCM.Insulator.Tubular(:sheath; thickness=t_pe, material=pe),
        LCM.Insulator.Tubular(:sheath; thickness=t_bed, material=pp),
    )
    armor_parts = (
        LCM.Conductor.Wires(
            :armor;
            wire_radius=d_wa / 2,
            num_wires=num_ar_wires,
            lay_ratio=10.0,
            material=steel,
        ),
        LCM.Insulator.Tubular(:armor; thickness=t_jac, material=pp),
    )
    datasheet_info = (;
        designation_code="(N)2XH(F)RK2Y",
        U0=500.0,
        U=525.0,
        conductor_cross_section=1600.0,
        screen_cross_section=1000.0,
        resistance=nothing,
        capacitance=nothing,
        inductance=nothing,
    )

    cable_id = "525kV_1600mm2"
    cable_design = only(LCM.CableBuilder(
        cable_id,
        core_parts,
        sheath_parts,
        armor_parts;
        nominal=datasheet_info,
    ))
    @test [component.id for component in cable_design.components] ==
          ["core", "sheath", "armor"]
    @test length(cable_design.components[1].conductor_group.layers) == 7
    @test length(cable_design.components[1].insulator_group.layers) == 4
    @test cable_design.nominal_data.conductor_cross_section == 1600.0

    constants = LCM.compute!(LCM.CableConstantsProblem(cable_design), LCM.Formulation())
    @test constants.R ≈ 1.337933694851607e-5
    @test constants.L ≈ 4.0671617988215346e-7
    @test constants.C ≈ 1.836871062680782e-10
    @test DataFrame(constants) isa DataFrame
    @test DataFrame(cable_design, :components) isa DataFrame

    cable_render = PB.make_render(
        LCM.DataModel.CablePreviewPlotSpec, cable_design)
    @test length(cable_render.figures) == 1

    mktempdir() do temp_dir
        library = LCM.CablesLibrary()
        add!(library, cable_design)
        library_file = joinpath(temp_dir, "cables_library.json")
        LCM.save(library; file_name=library_file)

        loaded_library = LCM.CablesLibrary()
        LCM.load!(loaded_library; file_name=library_file)
        loaded_design = get(loaded_library, cable_id)
        @test loaded_design.cable_id == cable_id

        frequencies = collect(10.0 .^ range(0, stop=6, length=61))
        positions = (
            LCM.at(
                x=-0.5,
                y=-1.0,
                phases=(:core => 1, :sheath => 0, :armor => 0),
            ),
            LCM.at(
                x=0.5,
                y=-1.0,
                phases=(:core => 2, :sheath => 0, :armor => 0),
            ),
        )
        problem = only(LCM.SystemBuilder(
            "525kV_1600mm2_bipole",
            loaded_design,
            positions;
            length=1000.0,
            temperature=20.0,
            earth=LCM.Earth(rho=100.0, eps_r=10.0, mu_r=1.0),
            frequencies,
        ))
        system = problem.system
        earth_model = problem.earth_props
        @test length(system.cables) == 2
        @test system.cables[1].horz == -0.5
        @test system.cables[2].horz == 0.5
        @test DataFrame(system) isa DataFrame
        @test DataFrame(earth_model) isa DataFrame

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
        atp_system = LCM.export_data(
            :atp,
            system,
            earth_model;
            file_name=joinpath(temp_dir, "atp_export.xml"),
        )
        @test isfile(pscad)
        @test filesize(pscad) > 0
        @test isfile(atp_system)
        @test filesize(atp_system) > 0

        formulation = LCM.Formulation()
        line_parameters = LCM.compute!(
            problem, formulation; options=(verbosity=0,))
        @test LCM.nfrequencies(line_parameters) == 61
        @test size(LCM.Z(line_parameters)) == (2, 2, 61)
        @test size(LCM.Y(line_parameters)) == (2, 2, 61)

        series_rl, shunt_gc = DataFrame(
            line_parameters, (LCM.R, LCM.L, LCM.G, LCM.C);
            length_unit=:kilo,
            tol=1e-9,
        )
        @test size(series_rl, 1) > 0
        @test size(shunt_gc, 1) > 0

        rlcg_render = PB.make_render(
            LCM.Engine.LineParameterPlotSpec,
            line_parameters;
            quantities=(LCM.R, LCM.L, LCM.G, LCM.C),
            xscale=:log10,
            length_unit=:kilo,
        )
        zy_render = PB.make_render(
            LCM.Engine.LineParameterPlotSpec,
            line_parameters;
            xscale=:log10,
            length_unit=:kilo,
        )
        series_zy_render = PB.make_render(
            LCM.Engine.LineParameterPlotSpec,
            LCM.Z(line_parameters);
            frequencies=LCM.frequencies(line_parameters),
            xscale=:log10,
            length_unit=:kilo,
        )
        shunt_zy_render = PB.make_render(
            LCM.Engine.LineParameterPlotSpec,
            LCM.Y(line_parameters);
            frequencies=LCM.frequencies(line_parameters),
            xscale=:log10,
            length_unit=:kilo,
        )
        @test length(rlcg_render.figures) == 2
        @test all(page -> length(page.views) == 2, rlcg_render.figures)
        @test length(zy_render.figures) == 2
        @test all(page -> length(page.views) == 2, zy_render.figures)
        @test all(
            view -> length(view.series) == LCM.nconductors(line_parameters)^2,
            Iterators.flatten(page.views for page in zy_render.figures),
        )
        @test getproperty.(first(zy_render.figures).views, :key) ==
              [(component=:Z_re,), (component=:Z_im,)]
        @test getproperty.(last(zy_render.figures).views, :key) ==
              [(component=:Y_re,), (component=:Y_im,)]
        @test length(series_zy_render.figures) == 1
        @test length(shunt_zy_render.figures) == 1
        @test length(only(series_zy_render.figures).views) == 2
        @test length(only(shunt_zy_render.figures).views) == 2
        @test getproperty.(only(series_zy_render.figures).views, :key) ==
              [(component=:Z_re,), (component=:Z_im,)]
        @test getproperty.(only(shunt_zy_render.figures).views, :key) ==
              [(component=:Y_re,), (component=:Y_im,)]
        @test_throws ArgumentError PB.make_render(
            LCM.Engine.LineParameterPlotSpec,
            line_parameters;
            mode=:ZY,
        )
        @test_throws ArgumentError PB.make_render(
            LCM.Engine.LineParameterPlotSpec,
            LCM.Z(line_parameters);
            frequencies=LCM.frequencies(line_parameters),
            quantities=(LCM.C,),
        )

        atp_zy = LCM.export_data(
            :atp,
            line_parameters;
            file_name=joinpath(temp_dir, "ZY_export.xml"),
            cable_system=system,
        )
        @test isfile(atp_zy)
        @test filesize(atp_zy) > 0

        transform, sequence_parameters = LCM.Fortescue(tol=1e-5)(line_parameters)
        @test size(transform) == (2, 2)
        @test LCM.nfrequencies(sequence_parameters) == 61
        series_zy, shunt_zy = DataFrame(
            sequence_parameters; length_unit=:kilo, tol=1e-9)
        @test size(series_zy, 1) > 0
        @test size(shunt_zy, 1) > 0

        sequence_render = PB.make_render(
            LCM.Engine.LineParameterPlotSpec,
            sequence_parameters;
            quantities=(LCM.R, LCM.L, LCM.G, LCM.C),
            xscale=:log10,
            length_unit=:kilo,
        )
        @test length(sequence_render.figures) == 2
        @test all(page -> length(page.views) == 2, sequence_render.figures)

        if get(ENV, "LINECABLEMODELS_TEST_PLOTTING", "false") == "true"
            extension = Base.get_extension(LCM, :LineCableModelsMakieExt)
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
            line_plots = extension.plot(
                line_parameters,
                (LCM.R, LCM.L, LCM.G, LCM.C);
                xscale=:log10,
                display_plot=false,
                controls=false,
            )
            @test !isempty(line_plots)
            @test all(handle -> isempty(handle.controls), line_plots)
        end
    end
end
