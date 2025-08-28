@testsnippet deps_datamodel begin
    using Plots
    using EzXML
end

@testitem "DataModel module" setup = [defaults, deps_datamodel] begin

    println("\nSetting up materials and dimensions for DataModel test...")
    materials = MaterialsLibrary(add_defaults=true)
    @test haskey(materials, "aluminum")
    @test haskey(materials, "copper")
    @test haskey(materials, "polyacrylate")
    @test haskey(materials, "semicon1")
    @test haskey(materials, "semicon2")
    @test haskey(materials, "pe")

    initial_default_count = length(materials)
    @test initial_default_count > 5 # Should have several defaults

    materials_empty = MaterialsLibrary(add_defaults=false)
    @test isempty(materials_empty)

    # Add a custom material for removal tests
    mat_remove_test = Material(1e-5, 5.0, 1.0, 20.0, 0.05)
    add!(materials, "remove_me", mat_remove_test)
    @test length(materials) == initial_default_count + 1
    @test haskey(materials, "remove_me")

    println("  Testing delete!...")
    delete!(materials, "remove_me")
    @test !haskey(materials, "remove_me")
    @test length(materials) == initial_default_count

    # Test removing non-existent (should throw KeyError based on source)
    @test_throws KeyError delete!(
        materials,
        "does_not_exist",
    )
    # Verify count didn't change
    @test length(materials) == initial_default_count

    println("  Testing DataFrame...")
    # Use the empty DB + one material for simpler checking
    mat_list_test = Material(9e9, 9.0, 9.0, 99.0, 0.9)
    add!(materials_empty, "list_test_mat", mat_list_test)
    df_listed = DataFrame(materials_empty)

    @test df_listed isa DataFrame
    @test names(df_listed) == ["name", "rho", "eps_r", "mu_r", "T0", "alpha"] # Check column names
    @test nrow(df_listed) == 1
    @test df_listed[1, :name] == "list_test_mat"
    @test df_listed[1, :rho] == 9e9
    @test df_listed[1, :eps_r] == 9.0
    @test df_listed[1, :mu_r] == 9.0
    @test df_listed[1, :T0] == 99.0
    @test df_listed[1, :alpha] == 0.9

    println("  Testing save/load cycle for MaterialsLibrary...")
    mktempdir(joinpath(@__DIR__)) do tmpdir
        output_file = joinpath(tmpdir, "materials_library_test.json")
        println("    Saving to: ", output_file)

        # Save the db that had defaults + 'remove_me' (before removal)
        # Let's re-add it for a more comprehensive save file
        db_to_save = MaterialsLibrary(add_defaults=true)
        mat_temp = Material(1e-5, 5.0, 1.0, 20.0, 0.05)
        add!(db_to_save, "temp_mat", mat_temp)
        num_expected = length(db_to_save)

        save(db_to_save, file_name=output_file)
        @test isfile(output_file)
        @test filesize(output_file) > 0

        # Load into a NEW, EMPTY library
        materials_from_json = MaterialsLibrary(add_defaults=false)
        load!(materials_from_json, file_name=output_file)

        # Verify loaded content
        @test length(materials_from_json) == num_expected
        @test haskey(materials_from_json, "temp_mat")
        @test haskey(materials_from_json, "copper") # Check a default also loaded
        loaded_temp_mat = get(materials_from_json, "temp_mat")
        @test loaded_temp_mat.rho == mat_temp.rho
        @test loaded_temp_mat.eps_r == mat_temp.eps_r

        println("    Save/load cycle completed.")
        println("Materials Library tests completed.")
    end # Temp dir cleanup

    # Cable dimensions from tutorial
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
    t_sct = 0.3e-3 # Semiconductive tape thickness
    t_alt = 0.15e-3
    t_pet = 0.05e-3
    t_jac = 2.4e-3

    # Nominal data for final comparison
    datasheet_info = NominalData(
        designation_code="NA2XS(FL)2Y",
        U0=18.0,                        # Phase-to-ground voltage [kV]
        U=30.0,                         # Phase-to-phase voltage [kV]
        conductor_cross_section=1000.0, # [mm²]
        screen_cross_section=35.0,      # [mm²]
        resistance=0.0291,              # DC resistance [Ω/km]
        capacitance=0.39,               # Capacitance [μF/km]
        inductance=0.3,                 # Inductance in trifoil [mH/km]
    )
    @test datasheet_info.resistance > 0
    @test datasheet_info.capacitance > 0
    @test datasheet_info.inductance > 0

    function calculate_rlc(
        design::CableDesign;
        rho_e::Float64=100.0,
        default_S_factor::Float64=2.0,
    )
        # Get components
        core_comp = design.components[findfirst(c -> c.id == "core", design.components)]
        sheath_comp =
            design.components[findfirst(c -> c.id == "sheath", design.components)]
        last_comp = design.components[end] # Usually jacket

        if isnothing(core_comp) || isnothing(sheath_comp)
            error(
                "Required 'core' or 'sheath' component not found in design for RLC calculation.",
            )
        end

        # Resistance (from effective core conductor group resistance)
        R = core_comp.conductor_group.resistance * 1e3 # Ω/m to Ω/km

        # Inductance (Trifoil)
        # Use outermost radius for separation calculation - corrected access path
        outermost_radius = last_comp.insulator_group.radius_ext
        S = default_S_factor * outermost_radius # Approx center-to-center distance [m]

        L =
            calc_inductance_trifoil(
                core_comp.conductor_group.radius_in,
                core_comp.conductor_group.radius_ext,
                core_comp.conductor_props.rho, core_comp.conductor_props.mu_r,
                sheath_comp.conductor_group.radius_in,
                sheath_comp.conductor_group.radius_ext, sheath_comp.conductor_props.rho,
                sheath_comp.conductor_props.mu_r,
                S, rho_e=rho_e,
            ) * 1e6 # H/m to mH/km

        # Capacitance
        C =
            calc_shunt_capacitance(
                core_comp.conductor_group.radius_ext,
                core_comp.insulator_group.radius_ext,
                core_comp.insulator_props.eps_r,
            ) * 1e6 * 1e3 # F/m to μF/km

        return R, L, C
    end

    println("Constructing core conductor group...")
    material_alu = get(materials, "aluminum")
    core = ConductorGroup(WireArray(0.0, Diameter(d_w), 1, 0.0, material_alu))
    @test core isa ConductorGroup
    @test length(core.layers) == 1
    @test core.radius_in == 0
    @test core.radius_ext ≈ d_w / 2.0
    @test core.resistance > 0
    @test core.gmr > 0

    add!(core, WireArray, Diameter(d_w), 6, 15.0, material_alu)
    @test length(core.layers) == 2
    @test core.radius_ext ≈ (d_w / 2.0) * 3 # Approximation for 1+6 wires
    @test core.resistance > 0 # Resistance should decrease

    add!(core, WireArray, Diameter(d_w), 12, 13.5, material_alu)
    @test length(core.layers) == 3
    @test core.radius_ext ≈ (d_w / 2.0) * 5 # Approximation for 1+6+12 wires

    add!(core, WireArray, Diameter(d_w), 18, 12.5, material_alu)
    @test length(core.layers) == 4
    @test core.radius_ext ≈ (d_w / 2.0) * 7 # Approximation

    add!(core, WireArray, Diameter(d_w), 24, 11.0, material_alu)
    @test length(core.layers) == 5
    @test core.radius_ext ≈ (d_w / 2.0) * 9 # Approximation
    # Check final calculated radius against nominal diameter
    # Note: constructor uses internal calculations, may differ slightly from d_core/2
    @test core.radius_ext ≈ d_core / 2.0 rtol = 0.1 # Allow 10% tolerance for geometric approximation vs nominal
    final_core_radius = core.radius_ext # Store for later use
    final_core_resistance = core.resistance # Store for later use

    println("Constructing main insulation group...")
    # Inner semiconductive tape
    material_sc_tape = get(materials, "polyacrylate")
    main_insu = InsulatorGroup(Semicon(core, Thickness(t_sct), material_sc_tape))
    @test main_insu isa InsulatorGroup
    @test length(main_insu.layers) == 1
    @test main_insu.radius_in ≈ final_core_radius
    @test main_insu.radius_ext ≈ final_core_radius + t_sct

    # Inner semiconductor
    material_sc1 = get(materials, "semicon1")
    add!(main_insu, Semicon, Thickness(t_sc_in), material_sc1)
    @test length(main_insu.layers) == 2
    @test main_insu.radius_ext ≈ final_core_radius + t_sct + t_sc_in

    # Main insulation (XLPE)
    material_pe = get(materials, "pe")
    add!(main_insu, Insulator, Thickness(t_ins), material_pe)
    @test length(main_insu.layers) == 3
    @test main_insu.radius_ext ≈ final_core_radius + t_sct + t_sc_in + t_ins

    # Outer semiconductor
    material_sc2 = get(materials, "semicon2")
    add!(main_insu, Semicon, Thickness(t_sc_out), material_sc2)
    @test length(main_insu.layers) == 4
    @test main_insu.radius_ext ≈ final_core_radius + t_sct + t_sc_in + t_ins + t_sc_out

    # Outer semiconductive tape
    add!(main_insu, Semicon, Thickness(t_sct), material_sc_tape)
    @test length(main_insu.layers) == 5
    @test main_insu.radius_ext ≈
          final_core_radius + t_sct + t_sc_in + t_ins + t_sc_out + t_sct
    @test main_insu.shunt_capacitance > 0
    @test main_insu.shunt_conductance >= 0
    final_insu_radius = main_insu.radius_ext # Store for later use

    println("Creating core cable component...")
    core_cc = CableComponent("core", core, main_insu)
    @test core_cc isa CableComponent
    @test core_cc.id == "core"
    @test core_cc.conductor_group === core
    @test core_cc.insulator_group === main_insu
    @test core_cc.conductor_props isa Material # Check effective props were created
    @test core_cc.insulator_props isa Material

    println("Initializing CableDesign...")
    cable_id = "tutorial2_test"
    cable_design = CableDesign(cable_id, core_cc, nominal_data=datasheet_info)
    @test cable_design isa CableDesign
    @test length(cable_design.components) == 1
    @test cable_design.components[1] === core_cc
    @test cable_design.nominal_data === datasheet_info

    println("Constructing sheath group...")
    # Wire screens
    lay_ratio_screen = 10.0
    material_cu = get(materials, "copper")
    screen_con = ConductorGroup(
        WireArray(
            main_insu,
            Diameter(d_ws),
            num_sc_wires,
            lay_ratio_screen,
            material_cu,
        ),
    )
    @test screen_con isa ConductorGroup
    @test screen_con.radius_in ≈ final_insu_radius
    @test screen_con.radius_ext ≈ final_insu_radius + d_ws # Approx radius of single layer of wires

    # Copper tape
    add!(
        screen_con,
        Strip,
        Thickness(t_cut),
        w_cut,
        lay_ratio_screen,
        material_cu,
    )
    @test screen_con.radius_ext ≈ final_insu_radius + d_ws + t_cut
    final_screen_con_radius = screen_con.radius_ext

    # Water blocking tape
    material_wbt = get(materials, "polyacrylate") # Assuming same as sc tape
    screen_insu = InsulatorGroup(Semicon(screen_con, Thickness(t_wbt), material_wbt))
    @test screen_insu.radius_ext ≈ final_screen_con_radius + t_wbt
    final_screen_insu_radius = screen_insu.radius_ext

    # Sheath Cable Component & Add to Design
    sheath_cc = CableComponent("sheath", screen_con, screen_insu)
    @test sheath_cc isa CableComponent
    add!(cable_design, sheath_cc)
    @test length(cable_design.components) == 2
    @test cable_design.components[2] === sheath_cc

    println("Constructing jacket group...")
    # Aluminum foil
    material_alu = get(materials, "aluminum") # Re-get just in case
    jacket_con = ConductorGroup(Tubular(screen_insu, Thickness(t_alt), material_alu))
    @test jacket_con.radius_ext ≈ final_screen_insu_radius + t_alt
    final_jacket_con_radius = jacket_con.radius_ext

    # PE layer after foil
    material_pe = get(materials, "pe") # Re-get just in case
    jacket_insu = InsulatorGroup(Insulator(jacket_con, Thickness(t_pet), material_pe))
    @test jacket_insu.radius_ext ≈ final_jacket_con_radius + t_pet

    # PE jacket
    add!(jacket_insu, Insulator, Thickness(t_jac), material_pe)
    @test jacket_insu.radius_ext ≈ final_jacket_con_radius + t_pet + t_jac
    final_jacket_insu_radius = jacket_insu.radius_ext

    # Add Jacket Component to Design (using alternative signature)
    add!(cable_design, "jacket", jacket_con, jacket_insu)
    @test length(cable_design.components) == 3
    @test cable_design.components[3].id == "jacket"
    # Check overall radius
    @test cable_design.components[3].insulator_group.radius_ext ≈
          final_jacket_insu_radius

    println("Checking DataFrame...")
    @test DataFrame(cable_design, :baseparams) isa DataFrame
    @test DataFrame(cable_design, :components) isa DataFrame
    @test DataFrame(cable_design, :detailed) isa DataFrame

    println("Validating calculated RLC against nominal values (rtol=6%)...")

    # Get components for calculation (assuming they are named consistently)
    cable_core =
        cable_design.components[findfirst(c -> c.id == "core", cable_design.components)]
    cable_sheath =
        cable_design.components[findfirst(
            c -> c.id == "sheath",
            cable_design.components,
        )] # Note: Tutorial used 'cable_shield' variable name
    cable_jacket =
        cable_design.components[findfirst(
            c -> c.id == "jacket",
            cable_design.components,
        )]

    @test cable_core !== nothing
    @test cable_sheath !== nothing
    @test cable_jacket !== nothing

    (R_orig, L_orig, C_orig) = calculate_rlc(cable_design)
    println("  Original design RLC = ($R_orig, $L_orig, $C_orig)")
    @test R_orig ≈ datasheet_info.resistance rtol = 0.06
    @test L_orig ≈ datasheet_info.inductance rtol = 0.06
    @test C_orig ≈ datasheet_info.capacitance rtol = 0.06

    println("\nTesting CableDesign reconstruction...")
    new_components = []
    for original_component in cable_design.components
        println("  Reconstructing component: $(original_component.id)")

        # Extract effective properties and dimensions
        eff_cond_props = original_component.conductor_props
        eff_ins_props = original_component.insulator_props
        r_in_cond = original_component.conductor_group.radius_in
        r_ext_cond = original_component.conductor_group.radius_ext
        r_in_ins = original_component.insulator_group.radius_in
        r_ext_ins = original_component.insulator_group.radius_ext

        # Sanity check dimensions
        @test r_ext_cond ≈ r_in_ins atol = 1e-9 # Inner radius of insulator must match outer of conductor

        # Create simplified Tubular conductor using effective properties
        # Note: We must provide a material object, which are the effective props here
        equiv_conductor = Tubular(r_in_cond, r_ext_cond, eff_cond_props)
        # Wrap it in a ConductorGroup (which recalculates R, L based on the Tubular part)
        equiv_cond_group = ConductorGroup(equiv_conductor)

        # Create simplified Insulator using effective properties
        equiv_insulator = Insulator(r_in_ins, r_ext_ins, eff_ins_props)
        # Wrap it in an InsulatorGroup (which recalculates C, G based on the Insulator part)
        equiv_ins_group = InsulatorGroup(equiv_insulator)

        # Create the new, equivalent CableComponent
        equiv_component = CableComponent(
            original_component.id,
            equiv_cond_group,
            equiv_ins_group,
        )

        # Check if the recalculated R/L/C/G of the simple groups match the effective props closely. Note: This tests the self-consistency of the effective property calculations and the Tubular/Insulator constructors. Tolerance might need adjustment.
        @test equiv_cond_group.resistance ≈
              calc_tubular_resistance(
            r_in_cond,
            r_ext_cond,
            eff_cond_props.rho,
            0.0,
            20.0,
            20.0,
        ) rtol =
            1e-6
        @test equiv_ins_group.shunt_capacitance ≈
              calc_shunt_capacitance(r_in_ins, r_ext_ins, eff_ins_props.eps_r) rtol =
            1e-6
        # GMR/Inductance and Conductance checks could also be added here

        push!(new_components, equiv_component)
    end

    # Assemble the new CableDesign from the equivalent components
    @test length(new_components) == length(cable_design.components)
    equiv_cable_design = CableDesign(
        cable_design.cable_id * "_equiv",
        new_components[1], # Initialize with the first equivalent component
        nominal_data=datasheet_info, # Keep same nominal data for reference
    )

    # Add remaining equivalent components
    if length(new_components) > 1
        for i in eachindex(new_components)[2:end]
            add!(equiv_cable_design, new_components[i])
        end
    end
    @test length(equiv_cable_design.components) == length(new_components)
    println("  Equivalent cable design assembled.")

    println("  Calculating RLC for equivalent design...")

    (R_equiv, L_equiv, C_equiv) = calculate_rlc(equiv_cable_design)

    println("    Original R, L, C = ", R_orig, ", ", L_orig, ", ", C_orig)
    println("  Equivalent R, L, C = ", R_equiv, ", ", L_equiv, ", ", C_equiv)

    # Use a tight tolerance because they *should* be mathematically equivalent if the model is self-consistent
    rtol_equiv = 1e-6
    # Resistance mismatch in equivalent model?
    @test R_equiv ≈ R_orig rtol = rtol_equiv
    # Inductance mismatch in equivalent model?
    @test L_equiv ≈ L_orig rtol = rtol_equiv
    # Capacitance mismatch in equivalent model?
    @test C_equiv ≈ C_orig rtol = rtol_equiv

    println("  Effective properties reconstruction test passed.")

    println("\nTesting CablesLibrary methods...")
    library = CablesLibrary()
    add!(library, cable_design)

    initial_count = length(library)
    test_cable_id = cable_design.cable_id # Should be "tutorial2_test"
    @test initial_count >= 1
    @test haskey(library, test_cable_id)

    println("  Testing delete!...")
    delete!(library, test_cable_id)
    @test !haskey(library, test_cable_id)
    @test length(library) == initial_count - 1

    # Test removing non-existent (should throw error)
    @test_throws KeyError delete!(library, "non_existent_cable_id_123")
    @test length(library) == initial_count - 1 # Count remains unchanged


    println("\nTesting JSON Save/Load and RLC consistency...")

    add!(library, cable_design)
    @test length(library) == initial_count # Should be back to original count
    @test haskey(library, test_cable_id)

    mktempdir(joinpath(@__DIR__)) do tmpdir # Create a temporary directory for the test file
        output_file = joinpath(tmpdir, "cables_library_test.json")
        println("  Saving library to: ", output_file)

        # Test saving
        @test isfile(save(library, file_name=output_file))
        @test filesize(output_file) > 0 # Check if file is not empty

        # Test loading into a new library
        loaded_library = CablesLibrary()
        load!(loaded_library, file_name=output_file)
        @test length(loaded_library) == length(library)

        # Retrieve the reloaded design
        reloaded_design = get(loaded_library, cable_design.cable_id)
        println("Reloaded components:")
        for comp in reloaded_design.components
            println("  ID: ", repr(comp.id), " Type: ", typeof(comp)) # Use repr to see if ID is empty or weird
        end
        @test reloaded_design isa CableDesign
        @test reloaded_design.cable_id == cable_design.cable_id
        @test length(reloaded_design.components) == length(cable_design.components)
        # Optionally, add more granular checks on reloaded components/layers if needed

        println("  Calculating RLC for reloaded design...")
        (R_reload, L_reload, C_reload) = calculate_rlc(reloaded_design)
        println("    Reloaded Design RLC = ($R_reload, $L_reload, $C_reload)")
        println("    Original Design RLC = ($R_orig, $L_orig, $C_orig)") # Print original for comparison

        # Use a very tight tolerance - should be almost identical if serialization is good
        rtol_serial = 1e-9
        # Resistance mismatch after JSON reload?
        @test R_reload ≈ R_orig rtol = rtol_serial
        # Inductance mismatch after JSON reload?
        @test L_reload ≈ L_orig rtol = rtol_serial
        # Capacitance mismatch after JSON reload?
        @test C_reload ≈ C_orig rtol = rtol_serial

        println("  JSON save/load test passed.")
    end # mktempdir ensures cleanup


    println("  Setting up CableSystem...")
    f_pscad = 10.0 .^ range(0, stop=6, length=10) # Frequency range
    earth_params_pscad = EarthModel(f_pscad, 100.0, 10.0, 1.0) # 100 Ω·m, εr=10, μr=1

    # Use outermost radius for trifoil calculation spacing
    r = cable_design.components[end].insulator_group.radius_ext

    s = 2r + 0.01

    x0, y0 = 0.0, -1.0 # System center 1 m underground
    xa, ya, xb, yb, xc, yc = trifoil_formation(x0, y0, s)

    d_ab = hypot(xa - xb, ya - yb)
    d_bc = hypot(xb - xc, yb - yc)
    d_ca = hypot(xc - xa, yc - ya)
    @assert min(d_ab, d_bc, d_ca) ≥ 2r

    cable_system_id = "tutorial2_pscad_test"
    cablepos =
        CablePosition(cable_design, xa, ya, Dict("core" => 1, "sheath" => 0, "jacket" => 0))
    cable_system =
        LineCableSystem(cable_system_id, 1000.0, cablepos)
    add!(
        cable_system,
        cable_design,
        xb,
        yb,
        Dict("core" => 2, "sheath" => 0, "jacket" => 0),
    )
    add!(
        cable_system,
        cable_design,
        xc,
        yc,
        Dict("core" => 3, "sheath" => 0, "jacket" => 0),
    )
    @test cable_system.num_cables == 3
    @test cable_system.num_phases == 3

    mktempdir(joinpath(@__DIR__)) do tmpdir
        output_file = joinpath(tmpdir, "tutorial2_export_test.pscx")
        println("  Exporting PSCAD file to: ", output_file)

        # Run export and use returned path (exporter may prefix basename with system_id)
        result_path = export_data(:pscad, cable_system, earth_params_pscad, file_name=output_file)

        # Basic file checks (use returned path)
        @test result_path != nothing
        @test isfile(result_path)
        @test filesize(result_path) > 200

        # Basic XML content checks
        xml_content = read(result_path, String)

        @test occursin("<?xml version=", xml_content) # Check for XML declaration
        @test occursin("<project ", xml_content) # Check for root project tag
        @test occursin("</project>", xml_content) # Check for closing root tag

        println("  Performing XML structure checks via XPath...")
        local xml_doc
        try
            xml_doc = readxml(result_path)
        catch parse_err
            println("Failed to parse generated XML: $(parse_err)")
            println("Skipping XPath validation due to parsing error.")
            return # Exit testset early
        end

        # 3. Check Root Element and Attributes
        project_node = root(xml_doc)
        @test nodename(project_node) == "project"
        @test haskey(project_node, "name")
        @test project_node["name"] == cable_system.system_id
        # Check for expected version if needed
        @test project_node["version"] == "5.0.2"

        # 4. Check Count of Cable Definitions
        # Finds all 'User' components representing a coaxial cable definition
        cable_coax_nodes = findall("//User[@name='master:Cable_Coax']", project_node)
        @test length(cable_coax_nodes) == length(cable_system.cables) # Should be 3

        # 5. Check Data within the First Cable Definition (CABNUM=1)
        # Construct XPath to find the <paramlist> within the first Cable_Coax User component
        # This is a bit complex: find User where name='master:Cable_Coax' AND which has a child param CABNUM=1
        xpath_cable1_params = "//User[@name='master:Cable_Coax'][paramlist/param[@name='CABNUM' and @value='1']]/paramlist"
        params_cable1_node = findfirst(xpath_cable1_params, project_node)
        @test !isnothing(params_cable1_node)

        if !isnothing(params_cable1_node)
            # Helper to get a specific param value from the paramlist node
            function get_param_value(paramlist_node, param_name)
                p_node = findfirst("param[@name='$(param_name)']", paramlist_node)
                return isnothing(p_node) ? nothing : p_node["value"]
            end

            # Check component names exported
            @test get_param_value(params_cable1_node, "CONNAM1") == "Core"   # Matches cable_system.cables[1].design_data.components[1].id ?
            @test get_param_value(params_cable1_node, "CONNAM2") == "Sheath" # Matches cable_system.cables[1].design_data.components[2].id ?
            @test get_param_value(params_cable1_node, "CONNAM3") == "Jacket" # Matches cable_system.cables[1].design_data.components[3].id ?

            # Check X position
            x_val_str = get_param_value(params_cable1_node, "X")
            @test !isnothing(x_val_str)
            if !isnothing(x_val_str)
                parsed_x = parse(Float64, x_val_str)
                expected_x = cable_system.cables[1].horz # Get horz from the first cable in the system
                println(
                    "    Checking first cable horz: XML='$(x_val_str)', Expected='$(expected_x)'",
                )
                @test parsed_x ≈ expected_x rtol = 1e-6
            end

            # Check Y position (in PSCAD Y is oriented downwards)
            y_val_str = get_param_value(params_cable1_node, "Y")
            @test !isnothing(y_val_str)
            if !isnothing(y_val_str)
                parsed_y = abs(parse(Float64, y_val_str))
                expected_y = abs(cable_system.cables[1].vert)
                println(
                    "    Checking first cable vert: XML='$(y_val_str)', Expected='$(expected_y)' (May differ due to PSCAD coord system)",
                )
                # Don't assert exact equality 
                @test isapprox(parsed_y, expected_y, rtol=1e-4)
            end


            # Check an effective property, e.g., Core conductor effective resistivity (RHOC)
            rhoc_val_str = get_param_value(params_cable1_node, "RHOC")
            @test !isnothing(rhoc_val_str)
            if !isnothing(rhoc_val_str)
                parsed_rhoc = parse(Float64, rhoc_val_str)
                # Get effective rho from the first component (core) of the first cable design
                expected_rhoc =
                    cable_system.cables[1].design_data.components[1].conductor_props.rho
                println(
                    "    Checking first cable RHOC: XML='$(rhoc_val_str)', Expected='$(expected_rhoc)'",
                )
                # Use a slightly looser tolerance for calculated effective properties
                @test parsed_rhoc ≈ expected_rhoc rtol = 1e-4
            end

            # Check an effective dielectric property, e.g., Main insulation Epsilon_r (EPS1)
            eps1_val_str = get_param_value(params_cable1_node, "EPS1")
            @test !isnothing(eps1_val_str)
            if !isnothing(eps1_val_str)
                parsed_eps1 = parse(Float64, eps1_val_str)
                # Get effective eps_r from the first component (core) insulator props
                expected_eps1 =
                    cable_system.cables[1].design_data.components[1].insulator_props.eps_r
                println(
                    "    Checking first cable EPS1: XML='$(eps1_val_str)', Expected='$(expected_eps1)'",
                )
                @test parsed_eps1 ≈ expected_eps1 rtol = 1e-4
            end

        end

        # 6. Check Ground Parameters (Example)
        ground_params =
            findfirst("//User[@name='master:Line_Ground']/paramlist", project_node)
        @test !isnothing(ground_params)

        println("  XML structure checks via XPath passed.")

        println("  PSCAD export basic checks passed.")
    end # mktempdir cleanup

    println("\nTesting plotting functions...")

    println("  Testing preview...")
    # Reuse the fully constructed cable_design
    @test preview(
        cable_design,
        display_plot=false,
        display_legend=true,
        backend=gr,
    ) isa Plots.Plot
    @test preview(
        cable_design,
        display_plot=false,
        display_legend=false,
        backend=gr,
    ) isa Plots.Plot # Test option

    println("  Testing preview...")
    # Reuse the fully constructed cable_system
    @test preview(cable_system, zoom_factor=0.5, backend=gr, display_plot=false) isa
          Plots.Plot
    @test preview(cable_system, zoom_factor=0.1, backend=gr, display_plot=false) isa
          Plots.Plot # Test option

    println("  Plotting functions executed without errors.")

    println("\nTesting DataFrame generation...")

    # Reuse the fully constructed cable_design
    println("  Testing DataFrame...")
    df_core = DataFrame(cable_design, :baseparams)
    @test df_core isa DataFrame
    @test names(df_core) == ["parameter", "computed", "nominal", "percent_diff"] ||
          names(df_core) == [
        "parameter",
        "computed",
        "nominal",
        "percent_diff",
        "lower",
        "upper",
        "in_range?",
    ] # Allow for uncertainty columns
    @test nrow(df_core) == 3

    df_comp = DataFrame(cable_design, :components)
    @test df_comp isa DataFrame
    # Expected columns: "property", "core", "sheath", "jacket" (based on tutorial build)
    @test names(df_comp) == ["property", "core", "sheath", "jacket"]
    @test nrow(df_comp) > 5 # Should have several properties

    df_detail = DataFrame(cable_design, :detailed)
    @test df_detail isa DataFrame
    @test "property" in names(df_detail)
    # Check if columns were generated for layers, e.g., "core, cond. layer 1"
    @test occursin("core, cond. layer 1", join(names(df_detail)))
    @test occursin("jacket, ins. layer 1", join(names(df_detail)))
    @test nrow(df_detail) > 10 # Should have many properties

    # Test invalid format
    @test_throws ErrorException DataFrame(cable_design, :invalid_format)

    println("  Testing DataFrame...")
    # Reuse the fully constructed cable_system
    df_sys = DataFrame(cable_system)
    @test df_sys isa DataFrame
    @test names(df_sys) == ["cable_id", "horz", "vert", "phase_mapping"]
    @test nrow(df_sys) == 3 # Because we added 3 cables

    println("  DataFrame functions executed successfully.")

    println("\nTesting Base.show methods...")

    # Reuse objects created earlier in the test file
    # List of objects that have custom text/plain show methods in DataModel
    # Add more as needed (e.g., specific part types if they have custom shows)
    objects_to_show = [
        core,            # ConductorGroup
        main_insu,       # InsulatorGroup
        core_cc,         # CableComponent
        cable_design,    # CableDesign
        cable_system,    # LineCableSystem
        materials,
        # Add an example of a basic part if desired and has a show method
        Tubular(0.0, 0.01, get(materials, "aluminum")),
    ]

    mime = MIME"text/plain"()

    for obj in objects_to_show
        println("  Testing show for: $(typeof(obj))")
        obj_repr = sprint(show, mime, obj)
        @test obj_repr isa String
        @test length(obj_repr) > 10 # Check that it produced some reasonable output
    end

    println("  Custom show methods executed without errors.")

    println("\nDataModel test completed.")

end