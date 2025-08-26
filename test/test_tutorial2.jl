@testitem "examples/tutorial2.jl tests" setup = [defaults] begin
    # Replicate the setup from the tutorial
    materials = MaterialsLibrary(add_defaults=true)

    # Cable dimensions from the tutorial
    num_co_wires = 61
    num_sc_wires = 49
    d_core = 38.1e-3
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

    # Test Core and Main Insulation construction
    @testset "core and main insulation" begin
        material_al = get(materials, "aluminum")
        @test material_al isa LineCableModels.Material
        core = ConductorGroup(WireArray(0.0, Diameter(d_w), 1, 0.0, material_al))
        add!(core, WireArray, Diameter(d_w), 6, 15.0, material_al)
        add!(core, WireArray, Diameter(d_w), 12, 13.5, material_al)
        add!(core, WireArray, Diameter(d_w), 18, 12.5, material_al)
        add!(core, WireArray, Diameter(d_w), 24, 11.0, material_al)

        @test length(core.layers) == 5
        @test isapprox(core.radius_ext * 2, 0.0423, atol=1e-4)

        material_poly = get(materials, "polyacrylate")
        material_sc1 = get(materials, "semicon1")
        material_pe = get(materials, "pe")
        material_sc2 = get(materials, "semicon2")

        main_insu = InsulatorGroup(Semicon(core, Thickness(t_sct), material_poly))
        add!(main_insu, Semicon, Thickness(t_sc_in), material_sc1)
        add!(main_insu, Insulator, Thickness(t_ins), material_pe)
        add!(main_insu, Semicon, Thickness(t_sc_out), material_sc2)
        add!(main_insu, Semicon, Thickness(t_sct), material_poly)

        @test length(main_insu.layers) == 5

        core_cc = CableComponent("core", core, main_insu)
        @test core_cc.id == "core"
    end

    # Build the full cable design step-by-step as in the tutorial
    # This also tests the constructors and `add!` methods implicitly
    material_al = get(materials, "aluminum")
    core = ConductorGroup(WireArray(0.0, Diameter(d_w), 1, 0.0, material_al))
    add!(core, WireArray, Diameter(d_w), 6, 15.0, material_al)
    add!(core, WireArray, Diameter(d_w), 12, 13.5, material_al)
    add!(core, WireArray, Diameter(d_w), 18, 12.5, material_al)
    add!(core, WireArray, Diameter(d_w), 24, 11.0, material_al)

    material_poly = get(materials, "polyacrylate")
    material_sc1 = get(materials, "semicon1")
    material_pe = get(materials, "pe")
    material_sc2 = get(materials, "semicon2")
    main_insu = InsulatorGroup(Semicon(core, Thickness(t_sct), material_poly))
    add!(main_insu, Semicon, Thickness(t_sc_in), material_sc1)
    add!(main_insu, Insulator, Thickness(t_ins), material_pe)
    add!(main_insu, Semicon, Thickness(t_sc_out), material_sc2)
    add!(main_insu, Semicon, Thickness(t_sct), material_poly)

    core_cc = CableComponent("core", core, main_insu)

    cable_id = "18kV_1000mm2"
    datasheet_info = NominalData(
        designation_code="NA2XS(FL)2Y", U0=18.0, U=30.0,
        conductor_cross_section=1000.0, screen_cross_section=35.0,
        resistance=0.0291, capacitance=0.39, inductance=0.3
    )
    cable_design = CableDesign(cable_id, core_cc, nominal_data=datasheet_info)

    @test length(cable_design.components) == 1
    @test cable_design.cable_id == cable_id

    material_cu = get(materials, "copper")
    lay_ratio = 10.0
    screen_con = ConductorGroup(WireArray(main_insu, Diameter(d_ws), num_sc_wires, lay_ratio, material_cu))
    add!(screen_con, Strip, Thickness(t_cut), w_cut, lay_ratio, material_cu)
    screen_insu = InsulatorGroup(Semicon(screen_con, Thickness(t_wbt), material_poly))
    sheath_cc = CableComponent("sheath", screen_con, screen_insu)
    add!(cable_design, sheath_cc)

    @test length(cable_design.components) == 2
    @test cable_design.components[2].id == "sheath"

    jacket_con = ConductorGroup(Tubular(screen_insu, Thickness(t_alt), material_al))
    jacket_insu = InsulatorGroup(Insulator(jacket_con, Thickness(t_pet), material_pe))
    add!(jacket_insu, Insulator, Thickness(t_jac), material_pe)
    add!(cable_design, "jacket", jacket_con, jacket_insu)

    @test length(cable_design.components) == 3
    @test cable_design.components[3].id == "jacket"

    @testset "calculated parameters vs hard-coded values" begin
        core_df = DataFrame(cable_design, :baseparams)

        # Hard-coded values from the tutorial
        expected_R = 0.0275677
        expected_L = 0.287184
        expected_C = 0.413357

        # Test R
        computed_R = core_df[core_df.parameter.=="R [Ω/km]", :computed][1]
        @test isapprox(computed_R, expected_R, atol=1e-5)

        # Test L
        computed_L = core_df[core_df.parameter.=="L [mH/km]", :computed][1]
        @test isapprox(computed_L, expected_L, atol=1e-5)

        # Test C
        computed_C = core_df[core_df.parameter.=="C [μF/km]", :computed][1]
        @test isapprox(computed_C, expected_C, atol=1e-5)
    end

    @testset "dataframes and library" begin
        # Test that DataFrame constructors do not throw errors
        @test DataFrame(cable_design, :components) isa DataFrame
        @test DataFrame(cable_design, :detailed) isa DataFrame

        # Test CablesLibrary functionality
        library = CablesLibrary()
        add!(library, cable_design)
        @test length(library) == 1
        @test DataFrame(library) isa DataFrame

        # Test saving and loading
        mktempdir(joinpath(@__DIR__)) do temp_dir
            output_file = joinpath(temp_dir, "cables_library.json")
            save(library, file_name=output_file)
            @test isfile(output_file)

            loaded_library = CablesLibrary()
            load!(loaded_library, file_name=output_file)
            @test length(loaded_library) == 1
            @test loaded_library.data[cable_id].cable_id == cable_id
        end
    end

    @testset "cable system and export" begin
        f = 10.0 .^ range(0, stop=6, length=10)
        earth_params = EarthModel(f, 100.0, 10.0, 1.0)
        @test DataFrame(earth_params) isa DataFrame

        x0, y0 = 0.0, -1.0
        xa, ya, xb, yb, xc, yc = trifoil_formation(x0, y0, 0.035)

        cablepos = CablePosition(cable_design, xa, ya, Dict("core" => 1, "sheath" => 0, "jacket" => 0))
        cable_system = LineCableSystem("18kV_1000mm2_trifoil", 1000.0, cablepos)
        add!(cable_system, cable_design, xb, yb, Dict("core" => 2, "sheath" => 0, "jacket" => 0))
        add!(cable_system, cable_design, xc, yc, Dict("core" => 3, "sheath" => 0, "jacket" => 0))

        @test length(cable_system.cables) == 3
        @test DataFrame(cable_system) isa DataFrame

        # Test PSCAD export
        mktempdir(joinpath(@__DIR__)) do temp_dir
            output_file = joinpath(temp_dir, "$(cable_system.system_id)_export.pscx")
            export_data(:pscad, cable_system, earth_params, file_name=output_file)
            @test isfile(output_file)
            # Check if file has content
            @test filesize(output_file) > 0
        end
    end

    @testset "preview functions" begin
        # Test that preview functions execute without error
        # Note: This does not check the plot content, only that they don't crash.
        @test preview(cable_design, display_plot=false) isa Any

        f = 10.0 .^ range(0, stop=6, length=10)
        earth_params = EarthModel(f, 100.0, 10.0, 1.0)
        x0, y0 = 0.0, -1.0
        xa, ya, xb, yb, xc, yc = trifoil_formation(x0, y0, 0.035)
        cablepos = CablePosition(cable_design, xa, ya, Dict("core" => 1, "sheath" => 0, "jacket" => 0))
        cable_system = LineCableSystem("18kV_1000mm2_trifoil", 1000.0, cablepos)
        add!(cable_system, cable_design, xb, yb, Dict("core" => 2, "sheath" => 0, "jacket" => 0))
        add!(cable_system, cable_design, xc, yc, Dict("core" => 3, "sheath" => 0, "jacket" => 0))

    end

    @testset "user error handling and robustness" begin
        # Test invalid material request
        @test get(materials, "unobtanium") === nothing

        # Test invalid geometric parameters
        @test_throws ArgumentError WireArray(0.0, Diameter(-1.0), 1, 0.0, material_al)
        @test_throws ArgumentError Insulator(core, Thickness(-1.0), material_pe)
        @test_throws ArgumentError WireArray(core, Diameter(d_w), 1, -1.0, material_al) # Negative lay ratio
        @test_throws ArgumentError Strip(core, Thickness(-0.1), w_cut, lay_ratio, material_cu)

        # Test empty object creation
        @test_throws ArgumentError ConductorGroup()
        @test_throws ArgumentError InsulatorGroup()
        @test_throws MethodError CableDesign("empty_cable")
        @test_throws MethodError LineCableSystem("empty_system", 1000.0)

        # Test trifoil formation with negative radius
        @test_throws AssertionError trifoil_formation(0.0, -1.0, -1.0)

        # Test adding a cable at an overlapping position in LineCableSystem
        x0, y0 = 0.0, -1.0
        xa, ya, xb, yb, xc, yc = trifoil_formation(x0, y0, 0.1) # Use a valid distance
        cablepos = CablePosition(cable_design, xa, ya, Dict("core" => 1))
        cable_system = LineCableSystem("overlap_test", 1000.0, cablepos)
        @test_throws ArgumentError add!(cable_system, cable_design, xa, ya, Dict("core" => 2))

        # Test conductor at interface
        @test_throws AssertionError CablePosition(cable_design, 0.0, 0.0, Dict("core" => 1))

        # Test invalid phase mapping
        @test_throws ArgumentError CablePosition(cable_design, 1.0, 1.0, Dict("non_existent_component" => 1))


        # Test exporting a system where some components are grounded (valid case)
        f = 10.0 .^ range(0, stop=6, length=10)
        earth_params = EarthModel(f, 100.0, 10.0, 1.0)
        cablepos_partially_grounded = CablePosition(cable_design, xa, ya, Dict("core" => 1, "sheath" => 0, "jacket" => 0))
        system_partially_grounded = LineCableSystem("partially_grounded_system", 1000.0, cablepos_partially_grounded)
        mktempdir(joinpath(@__DIR__)) do temp_dir
            output_file = joinpath(temp_dir, "partially_grounded_export.pscx")
            # This should run without error
            export_data(:pscad, system_partially_grounded, earth_params, file_name=output_file)
            @test isfile(output_file)
            @test filesize(output_file) > 0
        end
    end
end