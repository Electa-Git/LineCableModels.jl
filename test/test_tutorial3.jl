@testitem "examples/tutorial3.jl tests" setup = [defaults] begin


    mktempdir(joinpath(@__DIR__)) do tmpdir
        # Materials
        materials = MaterialsLibrary(add_defaults=true)
        lead = Material(21.4e-8, 1.0, 0.999983, 20.0, 0.00400)
        add!(materials, "lead", lead)
        steel = Material(13.8e-8, 1.0, 300.0, 20.0, 0.00450)
        add!(materials, "steel", steel)
        pp = Material(1e15, 2.8, 1.0, 20.0, 0.0)
        add!(materials, "pp", pp)

        @test haskey(materials, "lead")
        @test haskey(materials, "steel")
        @test haskey(materials, "pp")

        # Cable dimensions
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

        # Core and main insulation
        material_cu = get(materials, "copper")
        n = 6
        core = ConductorGroup(WireArray(0.0, Diameter(d_w), 1, 0.0, material_cu))
        add!(core, WireArray, Diameter(d_w), 1 * n, 11.0, material_cu)
        add!(core, WireArray, Diameter(d_w), 2 * n, 11.0, material_cu)
        add!(core, WireArray, Diameter(d_w), 3 * n, 11.0, material_cu)
        add!(core, WireArray, Diameter(d_w), 4 * n, 11.0, material_cu)
        add!(core, WireArray, Diameter(d_w), 5 * n, 11.0, material_cu)
        add!(core, WireArray, Diameter(d_w), 6 * n, 11.0, material_cu)

        material_sc1 = get(materials, "semicon1")
        main_insu = InsulatorGroup(Semicon(core, Thickness(t_sc_in), material_sc1))
        material_pe = get(materials, "pe")
        add!(main_insu, Insulator, Thickness(t_ins), material_pe)
        material_sc2 = get(materials, "semicon2")
        add!(main_insu, Semicon, Thickness(t_sc_out), material_sc2)
        material_pa = get(materials, "polyacrylate")
        add!(main_insu, Semicon, Thickness(t_wbt), material_pa)

        core_cc = CableComponent("core", core, main_insu)
        cable_id = "525kV_1600mm2"
        datasheet_info = NominalData(U=525.0, conductor_cross_section=1600.0)
        cable_design = CableDesign(cable_id, core_cc, nominal_data=datasheet_info)

        @test length(cable_design.components) == 1
        @test cable_design.components[1].id == "core"

        # Lead screen/sheath
        material_lead = get(materials, "lead")
        screen_con = ConductorGroup(Tubular(main_insu, Thickness(t_sc), material_lead))
        material_pe_sheath = get(materials, "pe")
        screen_insu = InsulatorGroup(Insulator(screen_con, Thickness(t_pe), material_pe_sheath))
        material_pp_bedding = get(materials, "pp")
        add!(screen_insu, Insulator, Thickness(t_bed), material_pp_bedding)
        sheath_cc = CableComponent("sheath", screen_con, screen_insu)
        add!(cable_design, sheath_cc)

        @test length(cable_design.components) == 2
        @test cable_design.components[2].id == "sheath"

        # Armor and outer jacket components
        lay_ratio = 10.0
        material_steel = get(materials, "steel")
        armor_con = ConductorGroup(WireArray(screen_insu, Diameter(d_wa), num_ar_wires, lay_ratio, material_steel))
        material_pp_jacket = get(materials, "pp")
        armor_insu = InsulatorGroup(Insulator(armor_con, Thickness(t_jac), material_pp_jacket))
        add!(cable_design, "armor", armor_con, armor_insu)

        @test length(cable_design.components) == 3
        @test cable_design.components[3].id == "armor"

        # Saving the cable design
        library = CablesLibrary()
        library_file = joinpath(tmpdir, "cables_library.json")
        add!(library, cable_design)
        save(library, file_name=library_file)

        loaded_library = CablesLibrary()
        load!(loaded_library, file_name=library_file)
        @test haskey(loaded_library, cable_id)
        reloaded_design = get(loaded_library, cable_id)
        @test reloaded_design.cable_id == cable_design.cable_id
        @test length(reloaded_design.components) == length(cable_design.components)

        # Defining a cable system
        f = 1e-3
        earth_params = EarthModel([f], 100.0, 10.0, 1.0)
        xp = -0.5
        xn = 0.5
        y0 = -1.0
        cablepos = CablePosition(cable_design, xp, y0, Dict("core" => 1, "sheath" => 0, "armor" => 0))
        cable_system = LineCableSystem("525kV_1600mm2_bipole", 1000.0, cablepos)
        add!(cable_system, cable_design, xn, y0, Dict("core" => 2, "sheath" => 0, "armor" => 0))

        @test length(cable_system.cables) == 2

        # FEM calculations
        problem = LineParametersProblem(cable_system, temperature=20.0, earth_props=earth_params, frequencies=[f])
        rho_g = earth_params.layers[end].rho_g[1]
        mu_g = earth_params.layers[end].mu_g[1]
        skin_depth_earth = abs(sqrt(rho_g / (1im * (2 * pi * f) * mu_g)))
        domain_radius = clamp(skin_depth_earth, 5.0, 5000.0)

        opts = (force_remesh=true, force_overwrite=true, plot_field_maps=false, mesh_only=false, save_path=joinpath(tmpdir, "fem_output"), keep_run_files=false, verbosity=0)

        formulation = FormulationSet(:FEM,
            impedance=Darwin(),
            admittance=Electrodynamics(),
            domain_radius=domain_radius,
            domain_radius_inf=domain_radius * 1.25,
            elements_per_length_conductor=1,
            elements_per_length_insulator=2,
            elements_per_length_semicon=1,
            elements_per_length_interfaces=5,
            points_per_circumference=16,
            mesh_size_min=1e-6,
            mesh_size_max=domain_radius / 5,
            mesh_size_default=domain_radius / 10,
            mesh_algorithm=5,
            mesh_max_retries=20,
            materials=materials,
            options=opts
        )

        workspace, line_params = compute!(problem, formulation)

        @test line_params isa LineParameters
        @test size(line_params.Z) == (6, 6, 1)
        @test size(line_params.Y) == (6, 6, 1)

        R = real(line_params.Z[1, 1, 1]) * 1000
        L = imag(line_params.Z[1, 1, 1]) / (2π * f) * 1e6
        C = imag(line_params.Y[1, 1, 1]) / (2π * f) * 1e9

        # Check if the results match hard-coded benchmarks
        @test isapprox(R, 0.01303, atol=1e-5)
        @test isapprox(L, 2.7606, atol=1e-4)
        @test isapprox(C, 0.1837, atol=1e-4)
    end
end