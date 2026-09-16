@testsnippet deps_export_atp begin
    using EzXML
end

# TODO: test if serialization works properly if uncertain types are used (Measurements)

@testitem "ImportExport / ATP / LineCableSystem LCC export" tags=[:integration] setup=[
    UseImportExportSupport, TestNumerics,
    TestFixtures, CableSystemFixture, deps_export_atp] begin

    # 1. ARRANGE & ACT: Run the export in a temporary directory
    mktempdir() do tmpdir
        output_file=joinpath(tmpdir, "atp_export_test.xml")
        result_path=export_data(:atp, cable_system, earth_props, file_name = output_file)
        expected_file=joinpath(
            dirname(output_file),
            "$(cable_system.system_id)_$(basename(output_file))"
        )

        # 2. ASSERT: Basic file checks (exporter prefixes basename with system_id)
        @test result_path == expected_file
        @test isfile(expected_file)

        # 3. ASSERT: General XML structure and LCC data
        doc=readxml(expected_file)
        root_node=root(doc)

        @test nodename(root_node) == "project"
        @test root_node["Application"] == "ATPDraw"

        # Find the main LCC component content node
        comp_content_node=findfirst("/project/objects/comp/comp_content", root_node)
        @test !isnothing(comp_content_node)

        # Verify general parameters like Length, Freq, and Ground Resistivity
        @test parse(
            Float64,
            findfirst("data[@Name='Length']", comp_content_node)["Value"]
        ) ≈ cable_system.line_length
        @test parse(Float64, findfirst("data[@Name='Freq']", comp_content_node)["Value"]) ≈
              problem_atp.frequencies[1]
        @test parse(
            Float64,
            findfirst("data[@Name='Grnd resis']", comp_content_node)["Value"]
        ) ≈ problem_atp.earth_props.layers[end].rho

        # 4. ASSERT: Detailed validation of ALL cables and conductors
        lcc_node=findfirst("/project/objects/comp/LCC", root_node)
        cable_header=findfirst("cable_header", lcc_node)
        cable_nodes=findall("cable", cable_header)

        @test length(cable_nodes) == num_phases

        # Loop through each cable exported in the XML and compare it to the source
        for (i, cable_node) in enumerate(cable_nodes)
            source_design=cable_system.designs[i]
            source_position=cable_system.positions[i]
            source_components=LineCableModels.DataModel.flatten(
                source_design, problem_atp.frequencies[1]
            )

            # Verify position of EACH cable
            @test parse(Float64, cable_node["PosX"]) ≈ source_position.x
            @test parse(Float64, cable_node["PosY"]) ≈ source_position.y

            # Verify the number of conductor components inside this cable
            num_components=length(source_components)
            @test parse(Int, cable_node["NumCond"]) == num_components

            conductor_nodes=findall("conductor", cable_node)
            @test length(conductor_nodes) == num_components

            # Loop through each conductor component within the cable
            for (j, conductor_node) in enumerate(conductor_nodes)
                source_component=source_components[j]
                conductor=source_component.conductor
                dielectric=source_component.dielectric

                expected_radius_in=conductor.r_in
                expected_radius_ext=conductor.r_ex
                expected_rho=conductor.material.rho
                expected_muC=conductor.material.mu_r
                expected_epsI=dielectric.material.eps_r
                expected_muI=dielectric.material.mu_r
                expected_Cext=dielectric.shunt_capacitance
                expected_Gext=dielectric.shunt_conductance

                # Assert that every attribute matches the expected value
                @test parse(Float64, conductor_node["Rin"]) ≈ expected_radius_in
                @test parse(Float64, conductor_node["Rout"]) ≈ expected_radius_ext
                @test parse(Float64, conductor_node["rho"]) ≈ expected_rho
                @test parse(Float64, conductor_node["muC"]) ≈ expected_muC
                @test parse(Float64, conductor_node["muI"]) ≈ expected_muI
                @test parse(Float64, conductor_node["epsI"]) ≈ expected_epsI
                @test parse(Float64, conductor_node["Cext"]) ≈ expected_Cext
                @test parse(Float64, conductor_node["Gext"]) ≈ expected_Gext
            end
        end
    end
end

@testitem "ImportExport / ATP / LineParameters ZY export" tags=[:integration] setup=[
    TestFixtures, deps_export_atp] begin
    source=TestFixtures.two_conductor_results()
    z,y=copy(Z(source)),copy(Y(source))
    # Distinct self/mutual entries and both signs expose ordering and parsing errors.
    z[1,2,:].*=-1
    y[2,1,:].*=-1
    parameters=LineParameters(z,y,frequencies(source))
    mktempdir() do directory
        path=joinpath(directory,"matrices.xml")
        @test export_data(:atp,parameters;file_name=path)==path
        document=readxml(path)
        node=root(document)
        @test nodename(node)=="ZY"
        @test parse(Int,node["NumPhases"])==2
        @test node["ZFmt"]=="R+Xi"
        @test node["YFmt"]=="G+Bi"
        pattern=r"^([+-]?[0-9.]+E[+-][0-9]+)([+-][0-9.]+E[+-][0-9]+)i$"
        for (name,expected) in (("Z",z),("Y",y))
            blocks=findall(name,node)
            @test length(blocks)==length(frequencies(source))
            for (k,block) in enumerate(blocks)
                @test parse(Float64,block["Freq"])==frequencies(source)[k]
                rows=split(strip(nodecontent(block)),'\n')
                @test length(rows)==2
                for (i,row) in enumerate(rows)
                    entries=split(row,',')
                    @test length(entries)==2
                    for (j,entry) in enumerate(entries)
                        parsed=match(pattern,strip(entry))
                        @test parsed!==nothing
                        parsed===nothing && error("Malformed ATP $name entry [$i,$j,$k]: $entry")
                        @test parse(Float64,parsed.captures[1])==real(expected[i,j,k])
                        @test parse(Float64,parsed.captures[2])==imag(expected[i,j,k])
                    end
                end
            end
        end
    end
end
