@testitem "ImportExport / PSCAD / system structure and physical values" tags=[:integration] setup=[
    ImportExportTestSupport,
    UseImportExportSupport,
    TestFixtures
] begin
    using EzXML

    function parameters(node)
        Dict(parameter["name"]=>parameter["value"]
        for parameter in findall("./paramlist/param", node))
    end

    system=TestFixtures.three_phase_system()
    earth=EarthModel(100.0, 10.0, 1.0)

    mktempdir() do directory
        requested=joinpath(directory, "system.pscx")
        output=export_data(:pscad, system, earth; file_name = requested)
        @test output == joinpath(directory, "$(system.system_id)_system.pscx")
        @test isfile(output)
        @test filesize(output) > 200

        document=readxml(output)
        project=root(document)
        @test nodename(project) == "project"
        @test project["name"] == system.system_id

        definitions=findall("./definitions/Definition", project)
        @test length(definitions) == 3
        station=only(filter(node->node["name"]=="DS", definitions))
        main=only(filter(node->node["name"]=="Main", definitions))
        cable_definition=only(filter(
            node->node["name"]=="CableSystem", definitions))
        @test station["classid"] == "StationDefn"
        @test main["classid"] == "UserCmpDefn"
        @test cable_definition["classid"] == "RowDefn"

        components=findall("./schematic/User", cable_definition)
        @test length(components) == length(system.cables) + 2
        @test getindex.(components, "id") ==
              string.(100000010:(100000011 + length(system.cables)))

        frequency=only(filter(
            node->node["defn"]==
            "master:Line_FrePhase_Options", components))
        cables=filter(node->node["defn"]=="master:Cable_Coax", components)
        ground=only(filter(node->node["defn"]=="master:Line_Ground", components))
        @test length(cables) == length(system.cables)
        @test parameters(frequency) == Dict(
            "FS" => "0.5",
            "FE" => "1.0E6",
            "Numf" => "100",
            "DCCOR" => "1",
            "CPASS" => "0"
        )

        for (index, (node, position)) in enumerate(zip(cables, system.cables))
            values=parameters(node)
            @test values["CABNUM"] == string(index)
            @test values["Name"] == position.design_data.cable_id
            @test parse(Float64, values["X"]) ≈ position.horz rtol=1e-5
            @test parse(Float64, values["Y"]) ≈ abs(position.vert) rtol=1e-5
            @test values["OHC"] == "0"
            @test values["LL"] == string(2length(position.design_data.components) - 1)
            @test values["CONNAM1"] ==
                  uppercasefirst(first(position.design_data.components).id)
            @test parse(Float64, values["R1"]) ≈
                  first(position.design_data.components).conductor_group.r_in rtol=1e-5
        end

        ground_values=parameters(ground)
        @test ground_values["GRRES"] == "100.0"
        @test ground_values["GPERM"] == "1.0"
        @test ground_values["GRP"] == "10.0"

        @test only(findall("./schematic", cable_definition))["classid"] == "RowCanvas"

        station_wire=only(findall("./schematic/Wire", station))
        @test station_wire["classid"] == "Branch"
        @test station_wire["name"] == "Main"
        @test station_wire["defn"] == "Main"
        main_reference=only(findall("./User", station_wire))
        @test main_reference["name"] == "$(system.system_id):Main"
        @test main_reference["defn"] == "$(system.system_id):Main"

        main_canvas=only(findall("./schematic", main))
        @test main_canvas["classid"] == "UserCanvas"
        cable_wire=only(findall("./Wire", main_canvas))
        @test cable_wire["classid"] == "Cable"
        @test cable_wire["name"] == "$(system.system_id):CableSystem"
        cable_reference=only(findall("./User", cable_wire))
        @test cable_reference["defn"] == "$(system.system_id):CableSystem"
        instance_values=parameters(cable_reference)
        @test instance_values["Name"] == system.system_id
        @test parse(Float64, instance_values["Freq"]) == 50.0
        @test parse(Float64, instance_values["Length"]) == system.line_length / 1000

        hierarchy=only(findall("./hierarchy/call", project))
        @test hierarchy["link"] == "100000001"
        @test hierarchy["name"] == "$(system.system_id):DS"
        main_call=only(findall("./call", hierarchy))
        @test main_call["link"] == "100000003"
        cable_call=only(findall("./call", main_call))
        @test cable_call["link"] == "100000008"
        @test cable_call["view"] == "true"
        @test length(findall("//Wire", document)) == 2
        @test occursin("LineCableModels.jl", read(output, String))
    end
end
