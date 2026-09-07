@testitem "ImportExport / PSCAD / simplified physical inputs and circuits" tags=[:integration] begin
    using EzXML
    const DM = LineCableModels.DataModel
    epsilon0 = 8.8541878128e-12
    mktempdir() do directory
        # A minimal imported project, not an exported file used as its own
        # oracle. PSCAD resistance and capacitance inputs use Ω/km and μF/km.
        function write_project(overrides=Dict{String,String}(); binding="master:Cable_CoaxSimpl")
            document = parsexml("""
                <project name="Imported">
                  <definitions>
                    <Definition name="CableSystem"><schematic>
                      <User defn="master:Line_FrePhase_Options"><paramlist/></User>
                      <User defn="master:Line_Ground"><paramlist>
                        <param name="GRRES" value="100"/>
                        <param name="GRP" value="10"/>
                        <param name="GPERM" value="1"/>
                      </paramlist></User>
                      <User defn="physical"><paramlist/></User>
                    </schematic></Definition>
                    <Definition name="Main"><schematic>
                      <User defn="Imported:CableSystem"><paramlist>
                        <param name="Name" value="six_cables"/>
                        <param name="Length" value="3.25 [km]"/>
                      </paramlist></User>
                    </schematic></Definition>
                  </definitions>
                </project>
                """)
            physical = only(findall("//User[@defn='physical']", document))
            physical["defn"] = binding
            values = merge(Dict(
                "Name"=>"", "CABNUM"=>"1", "NC"=>"2", "D"=>"0.15",
                "X"=>"-0.5", "Y"=>"1.2", "FLT"=>"6D1 [Hz]", "LL"=>"3",
                "RorT"=>"0", "SemiCL"=>"0", "R1"=>"0", "R2"=>"0.004",
                "R3"=>"0.006", "R4"=>"0.0065", "R5"=>"0.0075",
                "DTC"=>"0", "DCRC"=>"0.35 [ohm/km]", "RHOC"=>"2.82e-8",
                "DTS"=>"1", "DCRS"=>"0.9", "RHOS"=>"2.2e-7",
                "PERMC"=>"1", "PERMS"=>"1.1", "DTI1"=>"0", "CI1"=>"0.3 [uF/km]",
                "EPS1"=>"2.3", "DTI2"=>"1", "CI2"=>"0.7", "EPS2"=>"3.2",
                "mu_r1"=>"1.2", "mu_r2"=>"1", "LT1"=>"0.02", "LT2"=>"0",
                "OHC"=>"1", "Y2"=>"2.0",
            ), overrides)
            parameters = only(findall("./paramlist", physical))
            for name in sort!(collect(keys(values)))
                node = addelement!(parameters, "param")
                node["name"] = name
                node["value"] = values[name]
            end
            path = joinpath(directory, "input.pscx")
            write(path, document)
            return path
        end

        earth, system = import_data(:pscad, write_project())
        @test ncables(system) == 6
        @test nphases(system) == 12
        @test system.line_length == 3250.0
        @test system.system_id == "six_cables"
        @test last(earth.layers).rho == 100.0
        @test [design.cable_id for design in system.designs] == ["cable$i" for i in 1:6]
        @test [position.x for position in system.positions] ≈ [-0.5 + 0.15i for i in 0:5]
        @test all(position -> position.y == -1.2, system.positions)
        @test collect(Iterators.flatten(system.connections)) == collect(1:12)
        for design in system.designs
            core, insulation, sheath, jacket = design.geometry.regions
            @test design.terminal_order == [:core, :sheath]
            @test DM.r_ex(core.primitive) == 0.004
            @test DM.r_in(sheath.primitive) == 0.006
            @test DM.r_ex(jacket.primitive) == 0.0075
            @test core.source.material.rho ≈ 0.35 / 1000 * pi * 0.004^2
            @test sheath.source.material.rho == 2.2e-7
            expected_epsilon = 0.3e-9 * log(0.006 / 0.004) / (2pi * epsilon0)
            @test insulation.source.material.eps_r ≈ expected_epsilon
            @test insulation.source.material.mu_r == 1.2
            @test inv(insulation.source.material.rho) ≈ 2pi * 60 * epsilon0 * expected_epsilon * 0.02
            @test jacket.source.material.eps_r == 3.2
            @test isinf(jacket.source.material.rho)
            @test all(region -> region.source.material.T0 == 20.0 &&
                iszero(region.source.material.alpha), design.geometry.regions)
        end
        # Exercise the complementary direct-permittivity / DC-resistance inputs.
        _, complementary = import_data(:pscad, write_project(Dict(
            "DTC"=>"1", "DTS"=>"0", "DTI1"=>"1", "DTI2"=>"0")))
        core, insulation, sheath, jacket = first(complementary.designs).geometry.regions
        @test core.source.material.rho == 2.82e-8
        @test insulation.source.material.eps_r == 2.3
        @test sheath.source.material.rho ≈ 0.9 / 1000 * pi * (0.0065^2 - 0.006^2)
        @test jacket.source.material.eps_r ≈ 0.7e-9 * log(0.0075 / 0.0065) / (2pi * epsilon0)
        _, selected_fields = import_data(:pscad, write_project(Dict(
            "RHOC"=>"unused", "EPS1"=>"unused", "DCRS"=>"unused", "CI2"=>"unused")))
        @test [region.source.material for region in first(selected_fields.designs).geometry.regions] ==
            [region.source.material for region in first(system.designs).geometry.regions]

        # LL=0 is the detailed bare-conductor form, not an insulated component.
        _, bare = import_data(:pscad, write_project(Dict("LL"=>"0"); binding="master:Cable_Coax"))
        @test ncables(bare) == nphases(bare) == 1
        @test first(bare.designs).terminal_order == [:conductor]
        @test length(first(bare.designs).geometry.regions) == 1
        @test first(bare.positions).y == 2.0

        for (field, value, error_type) in (
            ("NC", "0", DomainError), ("NC", "1.5", ArgumentError),
            ("RorT", "1", ArgumentError), ("SemiCL", "1", ArgumentError),
            ("LL", "2", ArgumentError), ("LL", "7", ArgumentError),
            ("DTC", "2", ArgumentError), ("DTI1", "2", ArgumentError),
            ("LT1", "-0.1", DomainError), ("LT1", "10.1", DomainError),
            ("FLT", "0", DomainError), ("FLT", "-60", DomainError),
            ("FLT", "Inf", DomainError), ("FLT", "text", ArgumentError),
            ("CABNUM", "2", ArgumentError),
        )
            @testset "$field=$value" begin
                @test_throws error_type import_data(:pscad, write_project(Dict(field=>value)))
            end
        end
        for binding in ("master:Cable_PipeType", "master:Line_Tower_3Phase")
            @test_throws ArgumentError import_data(:pscad, write_project(; binding))
        end
        for binding in ("master:Cable_Coax", "master:Cable_CoaxSimpl"), invalid in ("0", "-60", "Inf")
            path = write_project(Dict("FLT"=>invalid); binding)
            failure = try
                import_data(:pscad, path)
            catch exception
                exception
            end
            @test failure isa DomainError
            @test occursin("reference frequency must be positive and finite", sprint(showerror, failure))
        end

        # Multiple row definitions must not silently change which physical
        # system is imported. Explicit selection wins over the output flags.
        path = write_project()
        document = readxml(path)
        definitions = only(findall("/project/definitions", document))
        other = addelement!(definitions, "Definition")
        other["name"] = "Other"
        schematic = addelement!(other, "schematic")
        for binding in ("master:Line_FrePhase_Options", "master:Line_Ground")
            node = addelement!(schematic, "User")
            node["defn"] = binding
            addelement!(node, "paramlist")
        end
        write(path, document)
        @test_throws ArgumentError import_data(:pscad, path)
        for definition in ("CableSystem", "Imported:CableSystem")
            _, selected = import_data(:pscad, path; definition)
            @test ncables(selected) == 6
            @test selected.system_id == system.system_id
        end
        @test_throws ArgumentError import_data(:pscad, path; definition="Unknown")
        output = addelement!(only(findall(
            "//Definition[@name='CableSystem']/schematic/User[@defn='master:Line_FrePhase_Options']/paramlist",
            document)), "param")
        output["name"] = "Output"
        for enabled in ("1", " YES ", "enabled", "True")
            output["value"] = enabled
            write(path, document)
            @test ncables(last(import_data(:pscad, path))) == 6
        end
        other_output = addelement!(only(findall(
            "./schematic/User[@defn='master:Line_FrePhase_Options']/paramlist", other)), "param")
        other_output["name"] = "Output"
        other_output["value"] = "1"
        write(path, document)
        @test_throws ArgumentError import_data(:pscad, path)
        @test ncables(last(import_data(:pscad, path; definition="CableSystem"))) == 6
        other["name"] = "CableSystem"
        write(path, document)
        @test_throws ArgumentError import_data(:pscad, path; definition="CableSystem")
    end
end
