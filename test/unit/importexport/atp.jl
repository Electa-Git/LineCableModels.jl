@testsnippet deps_export_atp begin
    using EzXML
end

# TODO: test if serialization works properly if uncertain types are used (Measurements)

@testitem "ImportExport / ATP / LineCableSystem LCC export" tags=[:integration, :importexport] setup=[
    UseImportExportSupport, TestNumerics,
    TestFixtures, CableSystemFixture, deps_export_atp] begin
    # The lossless equivalent permittivity of a conductor's dielectric group: its layers in
    # series, eps_eq = log(b/a) / sum(log(r_ex/r_in) / eps_r).
    function lossless_permittivity(blueprint, index)
        layers=blueprint[index].layers
        isempty(layers) && return 0.0
        log(last(layers).r_ex/first(layers).r_in)/
            sum(layer->log(layer.r_ex/layer.r_in)/layer.material.eps_r, layers)
    end

    # 1. ARRANGE & ACT: run the export in a temporary directory
    mktempdir() do tmpdir
        output_file=joinpath(tmpdir, "atp_export_test.xml")
        result_path=export_data(:atp, cable_system, earth_props, file_name = output_file)
        expected_file=joinpath(
            dirname(output_file),
            "$(cable_system.system_id)_$(basename(output_file))"
        )

        # 2. ASSERT: basic file checks (exporter prefixes basename with system_id)
        @test result_path == expected_file
        @test isfile(expected_file)

        # 3. ASSERT: general XML structure and LCC data
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

        # 4. ASSERT: detailed validation of ALL cables and conductors
        lcc_node=findfirst("/project/objects/comp/LCC", root_node)
        cable_header=findfirst("cable_header", lcc_node)
        cable_nodes=findall("cable", cable_header)

        @test length(cable_nodes) == num_phases

        # Loop through each cable exported in the XML and compare it to the source
        for (i, cable_node) in enumerate(cable_nodes)
            source_design=cable_system.designs[i]
            source_position=cable_system.positions[i]
            blueprint=LineCableModels.DataModel.flatten(source_design)

            # Verify position of EACH cable
            @test parse(Float64, cable_node["PosX"]) ≈ source_position.x
            @test parse(Float64, cable_node["PosY"]) ≈ source_position.y

            # Verify the number of conductor components inside this cable
            num_components=length(blueprint.conductors)
            @test parse(Int, cable_node["NumCond"]) == num_components

            conductor_nodes=findall("conductor", cable_node)
            @test length(conductor_nodes) == num_components

            # The default formulation is lossless: the groups have no conductance, and
            # ATPDraw takes the capacitance from epsI alone.
            for (j, conductor_node) in enumerate(conductor_nodes)
                conductor=blueprint.conductors[j]
                @test parse(Float64, conductor_node["Rin"]) ≈ conductor.r_in
                @test parse(Float64, conductor_node["Rout"]) ≈ conductor.r_ex
                @test parse(Float64, conductor_node["rho"]) ≈ conductor.material.rho
                @test parse(Float64, conductor_node["muC"]) ≈ conductor.material.mu_r
                @test parse(Float64, conductor_node["muI"]) ≈ blueprint.mu_r_ins[j]
                @test parse(Float64, conductor_node["epsI"]) ≈
                      lossless_permittivity(blueprint, j) rtol=1e-12
                @test parse(Float64, conductor_node["Cext"]) == 0
                @test parse(Float64, conductor_node["Gext"]) == 0
            end
        end
    end
end

@testitem "ImportExport / ATP / LineParameters ZY export" tags=[:integration, :importexport] setup=[
    TestFixtures, deps_export_atp] begin
    source=TestFixtures.two_conductor_results()
    z,y=copy(Z(source)),copy(Y(source))
    # Distinct self and mutual entries and both signs expose ordering and parsing errors.
    z[1,2,:].*=-1
    y[2,1,:].*=-1
    parameters=LineParameters(z,y,frequencies(source); line_length=1.0)
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
