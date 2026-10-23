@testitem "ImportExport / PSCAD / system structure and physical values" tags=[:integration, :pscad] setup=[
    UseImportExportSupport,
    TestFixtures
] begin
    using EzXML

    function parameters(node)
        Dict(parameter["name"]=>parameter["value"]
        for parameter in findall("./paramlist/param", node))
    end
    # The shunt admittance of a conductor's dielectric group at the reference temperature:
    # its layers' admittances in series, each from the law of its own class.
    function group_admittance(blueprint, index, formulation, frequency)
        layers=blueprint[index].layers
        isempty(layers) && return 0.0im
        laws=(formulation.formulas.insulation_admittance,
            formulation.formulas.semicon_admittance)
        inv(sum(layers) do layer
            kappa=constitutive(laws, layer.material, frequency, layer.material.T0)
            log(layer.r_ex/layer.r_in)/(2π*kappa)
        end)
    end

    # A split dielectric with distinct laws exercises the actual multilayer
    # exporter. A generic coaxial recipe need not include this extra scenario.
    metal=TestFixtures.conductor_material()
    first_layer=Material(kind=:insulator,rho=1e8,eps_r=3.0)
    second_layer=Material(kind=:insulator,rho=2e8,eps_r=5.0)
    design=build(CableDesign,"export-layers",Stack(
        terminal(:core,Region(:metal,Disk(.005),metal)),
        Region(:inner_dielectric,Shell(.002),first_layer),
        Region(:outer_dielectric,Shell(.003),second_layer),
        terminal(:sheath,Region(:return,Shell(.001),metal)),
        Region(:cover,Shell(.001),first_layer)))
    system=build(LineCableSystem,[design,design,design],
        [(-.08,-1.0),(0.,-1.08),(.08,-1.0)];
        connections=[Dict(:core=>i,:sheath=>0) for i in 1:3],system_id="export-current")
    earth=EarthModel(100.0, 10.0, 1.0)

    mktempdir() do directory
        requested=joinpath(directory, "system.pscx")
        lossy=Formulation(insulation_admittance=:lossy, semicon_admittance=:lossy)
        output=export_data(:pscad, system, earth; file_name = requested,
            formulation = lossy)
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
        # The namespace keeps only PSCAD's characters of the system id. The line keeps the id.
        @test project["name"] == "export_current"

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
        @test length(components) == length(system.designs) + 2
        @test getindex.(components, "id") ==
              string.(100000010:(100000011 + length(system.designs)))

        frequency=only(filter(
            node->node["defn"]==
            "master:Line_FrePhase_Options", components))
        cables=filter(node->node["defn"]=="master:Cable_Coax", components)
        ground=only(filter(node->node["defn"]=="master:Line_Ground", components))
        @test length(cables) == length(system.designs)
        @test parameters(frequency) == Dict(
            "FS" => "0.5",
            "FE" => "1.0E6",
            "Numf" => "100",
            "DCCOR" => "1",
            "CPASS" => "0",
            "enablf" => "1",
            "shntcab" => "0.0"
        )

        blueprints=[LineCableModels.DataModel.flatten(design) for design in system.designs]
        @test any(length(range) > 1
        for blueprint in blueprints for range in blueprint.dielectric_ranges)
        outer_radii=("R2", "R4", "R6", "R8")
        insulation_radii=("R3", "R5", "R7", "R9")
        resistivities=("RHOC", "RHOS", "RHOA", "RHOO")
        conductor_mus=("PERMC", "PERMS", "PERMA", "PERMO")
        permittivities=("EPS1", "EPS2", "EPS3", "EPS4")
        insulation_mus=("PERM1", "PERM2", "PERM3", "PERM4")
        losses=("LT1", "LT2", "LT3", "LT4")
        omega=2π*50.0
        for (index, (node, design, position, blueprint)) in enumerate(zip(
            cables,
            system.designs,
            system.positions,
            blueprints
        ))
            values=parameters(node)
            conductors=blueprint.conductors
            @test values["CABNUM"] == string(index)
            @test values["Name"] == design.cable_id
            @test parse(Float64, values["X"]) ≈ position.x rtol=1e-5
            @test parse(Float64, values["Y"]) ≈ abs(position.y) rtol=1e-5
            @test values["OHC"] == "0"
            @test values["LL"] == string(2length(conductors) - 1)
            @test values["CONNAM1"] == uppercasefirst(String(first(conductors).terminal))
            @test parse(Float64, values["R1"]) ≈ first(conductors).r_in rtol=1e-5
            for (component_index, conductor) in enumerate(conductors)
                layers=blueprint[component_index].layers
                outer=isempty(layers) ? conductor.r_ex : last(layers).r_ex
                admittance=group_admittance(blueprint, component_index, lossy, 50.0)
                capacitance=imag(admittance)/omega
                @test parse(Float64, values[outer_radii[component_index]]) ≈
                      conductor.r_ex rtol=1e-5
                @test parse(Float64, values[insulation_radii[component_index]]) ≈
                      outer rtol=1e-5
                @test parse(Float64, values[resistivities[component_index]]) ≈
                      conductor.material.rho rtol=1e-5
                @test parse(Float64, values[conductor_mus[component_index]]) ≈
                      conductor.material.mu_r rtol=1e-5
                @test parse(Float64, values[permittivities[component_index]]) ≈
                      LineCableModels.DataModel.BaseParams.equivalent_eps(capacitance, outer,
                          conductor.r_ex) rtol=1e-5
                @test parse(Float64, values[insulation_mus[component_index]]) ≈
                      blueprint.mu_r_ins[component_index] rtol=1e-5
                expected_loss=iszero(capacitance) ? 0.0 :
                              min(real(admittance)/(omega*capacitance), 10.0)
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
        @test main_reference["name"] == "export_current:Main"
        @test main_reference["defn"] == "export_current:Main"

        main_canvas=only(findall("./schematic", main))
        @test main_canvas["classid"] == "UserCanvas"
        cable_wire=only(findall("./Wire", main_canvas))
        @test cable_wire["classid"] == "Cable"
        @test cable_wire["name"] == "export_current:CableSystem"
        cable_reference=only(findall("./User", cable_wire))
        @test cable_reference["defn"] == "export_current:CableSystem"
        instance_values=parameters(cable_reference)
        @test instance_values["Name"] == system.system_id
        @test parse(Float64, instance_values["Freq"]) == 50.0
        @test parse(Float64, instance_values["Length"]) == system.line_length / 1000

        hierarchy=only(findall("./hierarchy/call", project))
        @test hierarchy["link"] == "100000001"
        @test hierarchy["name"] == "export_current:DS"
        main_call=only(findall("./call", hierarchy))
        @test main_call["link"] == "100000003"
        cable_call=only(findall("./call", main_call))
        @test cable_call["link"] == "100000008"
        @test cable_call["view"] == "true"
        @test length(findall("//Wire", document)) == 2
        @test occursin("LineCableModels.jl", read(output, String))
    end
end

@testitem "ImportExport / PSCAD / four-component limit" tags=[:integration, :pscad] setup=[
    UseImportExportSupport
] begin
    using LineCableModels
    import LineCableModels.DataModel

    metal=Material(:conductor, 1.7e-8, 1.0, 1.0, 20.0, 0.0039)
    dielectric=Material(:insulator, 1.0e14, 2.3, 1.0, 20.0, 0.0)
    parts=let radius=0.0
        values=DataModel.AbstractCablePart[]
        for index in 1:5
            terminal=Symbol("component$index")
            push!(values,
                DataModel.Group(
                    terminal,
                    DataModel.Region(
                        Symbol(terminal, :_conductor),
                        DataModel.Annulus(radius, radius+1.0e-3),
                        metal
                    )
                ))
            push!(values,
                DataModel.Region(
                    Symbol(terminal, :_insulation),
                    DataModel.Annulus(radius+1.0e-3, radius+2.0e-3),
                    dielectric
                ))
            radius+=2.0e-3
        end
        values
    end
    design=build(
        CableDesign,
        "five-components",
        DataModel.Stack(parts)
    )
    connections=Dict(name=>index for (index, name) in enumerate(design.terminal_order))
    system=build(
        LineCableSystem,
        design,
        DataModel.Pose2(0.0, -1.0, 0.0);
        connections,
        system_id = "five-components",
        line_length = 1_000.0
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

@testitem "ImportExport / PSCAD / materialized round trip" tags=[:integration, :pscad] setup=[
    UseImportExportSupport,
    TestFixtures
] begin
    # The shunt admittance of a conductor's dielectric group at the reference temperature:
    # its layers' admittances in series, each from the law of its own class.
    function group_admittance(blueprint, index, formulation, frequency)
        layers=blueprint[index].layers
        isempty(layers) && return 0.0im
        laws=(formulation.formulas.insulation_admittance,
            formulation.formulas.semicon_admittance)
        inv(sum(layers) do layer
            kappa=constitutive(laws, layer.material, frequency, layer.material.T0)
            log(layer.r_ex/layer.r_in)/(2π*kappa)
        end)
    end
    # A split dielectric with distinct laws exercises the actual multilayer
    # exporter. A generic coaxial recipe need not include this extra scenario.
    metal=TestFixtures.conductor_material()
    first_layer=Material(kind=:insulator,rho=1e8,eps_r=3.0)
    second_layer=Material(kind=:insulator,rho=2e8,eps_r=5.0)
    design=build(CableDesign,"export-layers",Stack(
        terminal(:core,Region(:metal,Disk(.005),metal)),
        Region(:inner_dielectric,Shell(.002),first_layer),
        Region(:outer_dielectric,Shell(.003),second_layer),
        terminal(:sheath,Region(:return,Shell(.001),metal)),
        Region(:cover,Shell(.001),first_layer)))
    system=build(LineCableSystem,[design,design,design],
        [(-.08,-1.0),(0.,-1.08),(.08,-1.0)];
        connections=[Dict(:core=>i,:sheath=>0) for i in 1:3],system_id="export-current")
    earth=EarthModel(100.0, 10.0, 1.0)

    mktempdir() do directory
        lossy=Formulation(insulation_admittance=:lossy, semicon_admittance=:lossy)
        output=export_data(
            :pscad,
            system,
            earth;
            formulation = lossy,
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

        for (imported_design, original_design, imported_position,
            original_position, imported_connections, original_connections) in zip(
            imported_system.designs,
            system.designs,
            imported_system.positions,
            system.positions,
            imported_system.connections,
            system.connections
        )
            @test imported_position.x ≈ original_position.x rtol=1e-5
            @test imported_position.y ≈ original_position.y rtol=1e-5
            @test imported_connections == original_connections
            @test imported_design.cable_id == original_design.cable_id
            @test imported_design.terminal_order == original_design.terminal_order
            # The import keeps the capacitance and conductance of each group at 50 Hz,
            # under the laws of the export.
            actual=LineCableModels.DataModel.flatten(imported_design)
            expected=LineCableModels.DataModel.flatten(original_design)
            for (index, (a, e)) in enumerate(zip(actual.conductors, expected.conductors))
                @test a.r_in ≈ e.r_in rtol=1e-5
                @test a.r_ex ≈ e.r_ex rtol=1e-5
                @test a.material.rho ≈ e.material.rho rtol=1e-5
                @test a.material.mu_r ≈ e.material.mu_r rtol=1e-5
                @test last(actual[index].layers).r_ex ≈ last(expected[index].layers).r_ex rtol=1e-5
                @test group_admittance(actual, index, lossy, 50.0) ≈
                      group_admittance(expected, index, lossy, 50.0) rtol=1e-5
                @test actual.mu_r_ins[index] ≈ expected.mu_r_ins[index] rtol=1e-5
            end
        end

        @test_throws ArgumentError import_data(:pscad, output * ".xml")
        @test_throws ArgumentError import_data(
            :pscad,
            joinpath(directory, "missing.pscx")
        )
    end
end
