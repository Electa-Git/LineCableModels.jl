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

        blocked=joinpath(directory, "$(system.system_id)_blocked.pscx")
        mkdir(blocked)
        @test_throws EzXML.XMLError export_data(
            :pscad,
            system,
            earth;
            file_name = joinpath(directory, "blocked.pscx")
        )

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

        @test any(
            length(component.conductor_group.layers) > 1 ||
                length(component.insulator_group.layers) > 1
        for position in system.cables for component in position.design_data.components
        )
        outer_radii=("R2", "R4", "R6", "R8")
        insulation_radii=("R3", "R5", "R7", "R9")
        resistivities=("RHOC", "RHOS", "RHOA", "RHOO")
        conductor_mus=("PERMC", "PERMS", "PERMA", "PERMO")
        permittivities=("EPS1", "EPS2", "EPS3", "EPS4")
        insulation_mus=("PERM1", "PERM2", "PERM3", "PERM4")
        losses=("LT1", "LT2", "LT3", "LT4")
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
            for (component_index, component) in enumerate(position.design_data.components)
                @test parse(Float64, values[outer_radii[component_index]]) ≈
                      component.conductor_group.r_ex rtol=1e-5
                @test parse(Float64, values[insulation_radii[component_index]]) ≈
                      component.insulator_group.r_ex rtol=1e-5
                @test parse(Float64, values[resistivities[component_index]]) ≈
                      component.conductor_props.rho rtol=1e-5
                @test parse(Float64, values[conductor_mus[component_index]]) ≈
                      component.conductor_props.mu_r rtol=1e-5
                @test parse(Float64, values[permittivities[component_index]]) ≈
                      component.insulator_props.eps_r rtol=1e-5
                @test parse(Float64, values[insulation_mus[component_index]]) ≈
                      component.insulator_props.mu_r rtol=1e-5
                capacitance=component.insulator_group.shunt_capacitance
                expected_loss=iszero(capacitance) ? 0.0 :
                              min(
                    component.insulator_group.shunt_conductance/
                    (2π*50.0*capacitance),
                    10.0
                )
                @test parse(Float64, values[losses[component_index]])≈
                expected_loss rtol=1e-5 atol=1e-12
            end
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

@testitem "ImportExport / PSCAD / four-component limit" tags=[:integration] setup=[
    ImportExportTestSupport,
    UseImportExportSupport
] begin
    using LineCableModels
    import LineCableModels.DataModel

    metal=Material(:conductor, 1.7e-8, 1.0, 1.0, 20.0, 0.0039)
    dielectric=Material(:insulator, 1.0e14, 2.3, 1.0, 20.0, 0.0)
    components=let radius=0.0
        values=DataModel.CableComponent[]
        for index in 1:5
            conductor=DataModel.Tubular(radius, radius+1.0e-3, metal)
            insulation=DataModel.Insulator(
                radius+1.0e-3,
                radius+2.0e-3,
                dielectric
            )
            push!(values,
                DataModel.CableComponent(
                    "component$index",
                    DataModel.ConductorGroup(conductor),
                    DataModel.InsulatorGroup(insulation)
                ))
            radius+=2.0e-3
        end
        values
    end
    design=DataModel.CableDesign("five-components", components)
    connections=Dict(component.id=>index for (index, component) in enumerate(components))
    system=DataModel.LineCableSystem(
        "five-components",
        1_000.0,
        DataModel.CablePosition(design, 0.0, -1.0, connections)
    )
    earth=EarthModel(100.0, 1.0, 1.0)
    mktempdir() do directory
        @test_throws ArgumentError export_data(
            :pscad,
            system,
            earth;
            file_name = joinpath(directory, "invalid.pscx")
        )
        @test isempty(readdir(directory))
    end
end

@testitem "ImportExport / PSCAD / materialized round trip" tags=[:integration] setup=[
    ImportExportTestSupport,
    UseImportExportSupport,
    TestFixtures
] begin
    system=TestFixtures.three_phase_system()
    earth=EarthModel(100.0, 10.0, 1.0)

    mktempdir() do directory
        output=export_data(
            :pscad,
            system,
            earth;
            file_name = joinpath(directory, "system.pscx")
        )
        imported_earth, imported_system=import_data(:pscad, output)

        @test imported_earth isa EarthModel{Float64}
        @test imported_system isa LineCableSystem{Float64}
        @test imported_system.system_id == system.system_id
        @test imported_system.line_length ≈ system.line_length
        @test ncables(imported_system) == ncables(system)
        @test nphases(imported_system) == nphases(system)
        @test last(imported_earth.layers).rho ≈ last(earth.layers).rho
        @test last(imported_earth.layers).eps_r ≈ last(earth.layers).eps_r
        @test last(imported_earth.layers).mu_r ≈ last(earth.layers).mu_r

        for (imported, original) in zip(imported_system.cables, system.cables)
            @test imported.horz ≈ original.horz rtol=1e-5
            @test imported.vert ≈ original.vert rtol=1e-5
            @test imported.conn == original.conn
            @test imported.design_data.cable_id == original.design_data.cable_id
            @test getproperty.(imported.design_data.components, :id) ==
                  getproperty.(original.design_data.components, :id)
            for (actual, expected) in zip(
                imported.design_data.components,
                original.design_data.components
            )
                @test actual.conductor_group.r_in ≈ expected.conductor_group.r_in rtol=1e-5
                @test actual.conductor_group.r_ex ≈ expected.conductor_group.r_ex rtol=1e-5
                @test actual.insulator_group.r_ex ≈ expected.insulator_group.r_ex rtol=1e-5
                @test actual.conductor_props.rho ≈ expected.conductor_props.rho rtol=1e-5
                @test actual.conductor_props.mu_r ≈ expected.conductor_props.mu_r rtol=1e-5
                @test actual.insulator_props.eps_r ≈ expected.insulator_props.eps_r rtol=1e-5
                dielectric=only(actual.insulator_group.layers).material_props
                @test dielectric.rho > 0
                @test dielectric.mu_r ≈ expected.insulator_props.mu_r rtol=1e-5
            end
        end

        @test_throws ArgumentError import_data(:pscad, output * ".xml")
        @test_throws ArgumentError import_data(
            :pscad,
            joinpath(directory, "missing.pscx")
        )
    end
end
