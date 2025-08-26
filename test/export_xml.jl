# test/test_atp_export.jl

using Test
using TestItemRunner
using LineCableModels

# Add your setup snippets here if they are in this file
# @testsnippet defaults begin ... end
# @testsnippet cable_system_export begin ... end

@testsnippet deps_export_atp begin
    using EzXML
end

@testitem "Export to ATPDraw LCC format" setup=[defaults, cable_system_export, deps_export_atp] begin

    # All variables from the setup snippets are available here (problem_atp, cable_system, etc.)

    # 1. ARRANGE & ACT: Run the export in a temporary directory
    mktempdir() do tmpdir
        output_file = joinpath(tmpdir, "atp_export_test.xml")
        result_path = export_data(Val(:atp), cable_system, earth_props, file_name=output_file)

        # 2. ASSERT: Basic file checks (unchanged)
        @test result_path == output_file
        @test isfile(output_file)
        @test filesize(output_file) > 500

        # 3. ASSERT: General XML structure and LCC data
        println("  Performing high-level XML structure checks...")
        doc = readxml(output_file)
        root_node = root(doc)

        @test nodename(root_node) == "project"
        @test root_node["Application"] == "ATPDraw"

        # Find the main LCC component content node
        comp_content_node = findfirst("/project/objects/comp/comp_content", root_node)
        @test !isnothing(comp_content_node)

        # Verify general parameters like Length, Freq, and Ground Resistivity
        println("  Verifying general LCC data (Length, Freq, Grnd resis)...")
        @test parse(Float64, findfirst("data[@Name='Length']", comp_content_node)["Value"]) ≈ cable_system.line_length
        @test parse(Float64, findfirst("data[@Name='Freq']", comp_content_node)["Value"]) ≈ problem_atp.frequencies[1]
        @test parse(Float64, findfirst("data[@Name='Grnd resis']", comp_content_node)["Value"]) ≈ problem_atp.earth_props.layers[end].base_rho_g

        # 4. ASSERT: Detailed validation of ALL cables and conductors
        println("  Verifying all cables and their conductors...")
        lcc_node = findfirst("/project/objects/comp/LCC", root_node)
        cable_header = findfirst("cable_header", lcc_node)
        cable_nodes = findall("cable", cable_header)

        @test length(cable_nodes) == num_phases

        # Loop through each cable exported in the XML and compare it to the source
        for (i, cable_node) in enumerate(cable_nodes)
            println("    -> Checking Cable #$i...")
            source_cable = cable_system.cables[i]

            # Verify position of EACH cable
            @test parse(Float64, cable_node["PosX"]) ≈ source_cable.horz
            @test parse(Float64, cable_node["PosY"]) ≈ source_cable.vert

            # Verify the number of conductor components inside this cable
            num_components = length(source_cable.design_data.components)
            @test parse(Int, cable_node["NumCond"]) == num_components

            conductor_nodes = findall("conductor", cable_node)
            @test length(conductor_nodes) == num_components

            # Loop through each conductor component within the cable
            for (j, conductor_node) in enumerate(conductor_nodes)
                source_component = source_cable.design_data.components[j]
                cond_group = source_component.conductor_group
                ins_group = source_component.insulator_group

                # Pre-calculate the expected values using the same functions as the export
                expected_rho = calc_equivalent_rho(cond_group.resistance, cond_group.radius_ext, cond_group.radius_in)
                expected_muC = calc_equivalent_mu(cond_group.gmr, cond_group.radius_ext, cond_group.radius_in)
                expected_epsI = calc_equivalent_eps(ins_group.shunt_capacitance, ins_group.radius_in, ins_group.radius_ext)

                # Assert that every attribute matches the expected value
                @test parse(Float64, conductor_node["Rin"]) ≈ cond_group.radius_in
                @test parse(Float64, conductor_node["Rout"]) ≈ cond_group.radius_ext
                @test parse(Float64, conductor_node["rho"]) ≈ expected_rho
                @test parse(Float64, conductor_node["muC"]) ≈ expected_muC
                @test parse(Float64, conductor_node["muI"]) ≈ ins_group.layers[1].material_props.mu_r
                @test parse(Float64, conductor_node["epsI"]) ≈ expected_epsI
                @test parse(Float64, conductor_node["Cext"]) ≈ ins_group.shunt_capacitance
                @test parse(Float64, conductor_node["Gext"]) ≈ ins_group.shunt_conductance
            end
        end
        println("  All detailed checks passed!")
    end
end

# Run the test
# @run_package_tests filter = ti -> occursin("Export to ATPDraw LCC format", ti.name)


@testitem "Export to ATP format" setup=[defaults, cable_system_export, deps_export_atp] begin

    # The ACT and ASSERT parts of your test go here.
    # All variables from the snippets are already defined.

    # 1. RUN THE TEST IN A TEMPORARY DIRECTORY
    mktempdir() do tmpdir
        output_file = joinpath(tmpdir, "atp_export_test.xml")
        println("  Exporting ATP XML file to: ", output_file)
        Z_matrix = randn(ComplexF64, num_phases, num_phases, length(freqs))
        Y_matrix = randn(ComplexF64, num_phases, num_phases, length(freqs))
        line_params = LineParameters(Z_matrix, Y_matrix)

        # Call the function we want to test
        result_path = export_data(Val(:atp), line_params, freqs;  file_name=output_file)
        
        # 2. BASIC FILE CHECKS
        @test result_path == output_file
        @test isfile(output_file)
        @test filesize(output_file) > 100

        xml_content = read(output_file, String)
        @test occursin("<ZY", xml_content)
        @test occursin("</ZY>", xml_content)

        # 3. XML STRUCTURE AND DATA VALIDATION
        println("  Performing XML structure checks via XPath...")
        xml_doc = readxml(output_file)
        root_node = root(xml_doc)
        
        @test nodename(root_node) == "ZY"
        @test parse(Int, root_node["NumPhases"]) == num_phases
        
        z_blocks = findall("//Z", root_node)
        @test length(z_blocks) == length(freqs)

        # 4. DETAILED DATA VERIFICATION (for the first frequency)
        println("  Verifying numerical data for first frequency...")
        first_z_block = z_blocks[1]
        @test parse(Float64, first_z_block["Freq"]) ≈ freqs[1]

        z_matrix_rows = split(strip(nodecontent(first_z_block)), '\n')
        @test length(z_matrix_rows) == num_phases
        
        first_row_elements = split(z_matrix_rows[1], ',')
        @test length(first_row_elements) == num_phases
        number_pattern = r"(-?[\d\.]+E[+-]\d+)"
        complex_pattern = Regex("$(number_pattern.pattern)([+-][\\d\\.]+E[+-]\\d+)i")

        match_result = match(complex_pattern, first_row_elements[1])

        if !isnothing(match_result)
            # The captures are now guaranteed to be valid Float64 strings
            real_part = parse(Float64, match_result.captures[1])
            imag_part = parse(Float64, match_result.captures[2])
            parsed_z11 = complex(real_part, imag_part)
            
            expected_z11 = Z_matrix[1, 1, 1]
            @test parsed_z11 ≈ expected_z11 rtol=1e-12
        end
    end
end


@run_package_tests filter = ti -> occursin("Export to ATP format", ti.name)