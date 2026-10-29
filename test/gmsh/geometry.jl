function geometry_entities(path)
    points = Dict{Int, NTuple{3, Float64}}()
    entities = Dict(name => Dict{Int, Vector{Int}}()
    for name in ("Line", "Circle", "Curve Loop", "Plane Surface"))
    patches = 0
    for line in eachline(path)
        entity = match(r"^(Point|Line|Circle|Curve Loop|Plane Surface)\((\d+)\)\s*=\s*\{([^}]*)\};", line)
        if entity !== nothing
            name, tag, entries = entity.captures
            values = split(entries, ',')
            if name == "Point"
                points[parse(Int, tag)] = Tuple(parse.(Float64, values[1:3]))
            else
                entities[name][parse(Int, tag)] = parse.(Int, values)
            end
        end
        count = match(r"^FEMAnnularPatches = (\d+);", line)
        count === nothing || (patches = parse(Int, count[1]))
    end
    return (; points, entities, patches)
end

@testset "Gmsh geometry" begin
    extension = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    Geometry, FieldModel = extension.Geometry, extension.FieldModel
    options = ComputationOptions(GmshBackend)
    @testset "$case" for (case, reference) in REFERENCES
        problem = LineParametersProblem(reference.system, reference.earth;
            temperature = reference.temperature,
            frequencies = reference.frequencies[1:10:end])
        frozen = geometry_entities(joinpath(
            @__DIR__, "..", "references", case, "physical.geo"))
        gmsh.initialize()
        try
            geometry = build(Geometry, FieldModel(problem, FEM()), options)
            @test length(geometry.surfaces) == length(problem.system.regions)
            @test length(geometry.outlines) == length(problem.system.designs)
            @test length(geometry.patches) == frozen.patches
            mktempdir() do directory
                path = joinpath(directory, "physical.geo_unrolled")
                gmsh.write(path)
                actual = geometry_entities(path)
                @test length(actual.points) == length(frozen.points)
                @test keys(actual.points) == keys(frozen.points)
                difference = maximum(maximum(abs.(actual.points[tag] .- xyz))
                for (tag, xyz) in frozen.points)
                @test difference <= 1e-12
                for name in ("Line", "Circle", "Curve Loop", "Plane Surface")
                    @test length(actual.entities[name]) == length(frozen.entities[name])
                    @test actual.entities[name] == frozen.entities[name]
                end
                written = joinpath(directory, "written.geo")
                write(written, geometry)
                @test !isfile(written * "_unrolled")
                actual_lines = readlines(written)
                expected_lines = readlines(joinpath(
                    @__DIR__, "..", "references", case, "physical.geo"))
                @test length(actual_lines) == length(expected_lines)
                number = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?"
                largest = 0.0
                terminal = 0
                cable = 0
                differing = Tuple{Int, String, String}[]
                for (line, (actual_line, expected_line)) in
                    enumerate(zip(actual_lines, expected_lines))
                    if startswith(expected_line, "Physical Surface(\"LCM/terminal/")
                        terminal += 1
                        pattern = r"^Physical Surface\(\"([^\"]+)\",(\d+)\) = \{([^}]*)\};$"
                        name = "LCM/terminal/$(lpad(terminal, 4, '0'))/$(geometry.model.label[terminal])"
                    elseif startswith(expected_line, "Physical Curve(\"LCM/cable_contour/")
                        cable += 1
                        pattern = r"^Physical Curve\(\"([^\"]+)\",(\d+)\) = \{([^}]*)\};$"
                        name = "LCM/cable_contour/$(lpad(cable, 4, '0'))/$(geometry.model.problem.system.designs[cable].id)"
                    else
                        actual_numbers = [parse(Float64, value.match)
                                          for value in eachmatch(number, actual_line)]
                        expected_numbers = [parse(Float64, value.match)
                                            for value in eachmatch(number, expected_line)]
                        if replace(actual_line, number => "") !=
                           replace(expected_line, number => "") ||
                           length(actual_numbers) != length(expected_numbers)
                            push!(differing, (line, actual_line, expected_line))
                        end
                        if length(actual_numbers) == length(expected_numbers)
                            largest = max(largest,
                                maximum(abs.(actual_numbers .- expected_numbers); init = 0.0))
                        end
                        continue
                    end
                    actual_group = match(pattern, actual_line)
                    expected_group = match(pattern, expected_line)
                    if actual_group === nothing ||
                       actual_group[2] != expected_group[2] ||
                       parse.(Int, split(actual_group[3], ',')) !=
                       parse.(Int, split(expected_group[3], ',')) ||
                       actual_group[1] != name
                        push!(differing, (line, actual_line, expected_line))
                    end
                end
                @test isempty(differing)
                @test largest <= 1e-12
                for (name, values) in (("FEMReceiverX", [p.x for p in geometry.anchors]),
                    ("FEMReceiverY", [p.y for p in geometry.anchors]),
                    ("FEMReceiverMetalDimension", [p.dimension for p in geometry.anchors]))
                    line = only(filter(line -> startswith(line, name * "()"), expected_lines))
                    frozen_values = parse.(Float64, split(match(r"\{([^}]*)\}", line)[1], ','))
                    @test length(values) == length(frozen_values)
                    @test all(abs.(values .- frozen_values) .<= 1e-12)
                end
            end
        finally
            gmsh.model.remove()
            gmsh.finalize()
        end
    end

    @testset "Registry loops" begin
        @testset "Bare wires" begin
            reference = last(REFERENCES[1])
            problem = LineParametersProblem(reference.system, reference.earth;
                temperature = reference.temperature,
                frequencies = reference.frequencies[1:10:end])
            gmsh.initialize()
            try
                geometry = build(Geometry, FieldModel(problem, FEM()), options)
                # Calls to draw visit two regions and two outlines: four paths. Each wire's
                # circle serves its region and its outline.
                @test length(geometry.loops) == 2
            finally
                gmsh.model.remove()
                gmsh.finalize()
            end
        end
        @testset "Insulated wires" begin
            reference = last(REFERENCES[2])
            problem = LineParametersProblem(reference.system, reference.earth;
                temperature = reference.temperature,
                frequencies = reference.frequencies[1:10:end])
            gmsh.initialize()
            try
                geometry = build(Geometry, FieldModel(problem, FEM()), options)
                # Calls to draw visit two conductors and two outlines: four paths. Each of
                # the eight patches makes its own loop outside the registry.
                @test length(geometry.loops) == 4
            finally
                gmsh.model.remove()
                gmsh.finalize()
            end
        end
        @testset "Insulated wires without strips" begin
            reference = last(REFERENCES[2])
            problem = LineParametersProblem(reference.system, reference.earth;
                temperature = reference.temperature,
                frequencies = reference.frequencies[1:10:end])
            gmsh.initialize()
            try
                geometry = build(Geometry, FieldModel(problem, FEM()),
                    ComputationOptions(GmshBackend; overrides = (ThinLayerStrips = 0,)))
                # Calls to draw visit two conductors and two insulations with an outer path
                # and a hole each. With the two outlines, there are eight paths.
                # Each conductor's loop serves its region and its
                # insulation's hole. Each insulation's outer circle serves the outline.
                @test length(geometry.loops) == 4
                @test isempty(geometry.patches)
                @test all(==(1), length.(geometry.surfaces))
                @test length(gmsh.model.getEntities(2)) == 4
            finally
                gmsh.model.remove()
                gmsh.finalize()
            end
        end
    end

    @testset "Touching trefoil" begin
        reference = last(REFERENCES[2])
        design = first(reference.system.designs)
        R = 0.0425 + 1.0e-3
        placements = [Placement(-R, -1), Placement(R, -1),
            Placement(0, -1 + sqrt(3) * R)]
        system = build(CableSystem, fill(design, 3), placements;
            connections = [(core = 1,), (core = 2,), (core = 3,)],
            id = "touching_trefoil", line_length = 1.0)
        problem = LineParametersProblem(system, reference.earth;
            temperature = reference.temperature,
            frequencies = reference.frequencies[1:10:end])
        gmsh.initialize()
        try
            geometry = build(Geometry, FieldModel(problem, FEM()), options)
            @test isempty(geometry.holes.air)
            @test isempty(geometry.pockets.air)
            @test length(geometry.holes.earth) == 1
            @test length(geometry.pockets.earth) == 1
            mktempdir() do directory
                path = joinpath(directory, "trefoil.geo")
                write(path, geometry)
                lines = readlines(path)
                earth = only(filter(line -> startswith(line, "FEMEarthPockets()"), lines))
                air = only(filter(line -> startswith(line, "FEMAirPockets()"), lines))
                @test parse.(Int, split(match(r"\{([^}]*)\}", earth)[1], ',')) ==
                      geometry.pockets.earth
                @test isempty(match(r"\{([^}]*)\}", air)[1])
            end
        finally
            gmsh.model.remove()
            gmsh.finalize()
        end
    end

    @testset "Flat row" begin
        reference = last(REFERENCES[2])
        design = first(reference.system.designs)
        R = 0.0425 + 1.0e-3
        placements = [Placement(-2R, -1), Placement(0, -1), Placement(2R, -1)]
        system = build(CableSystem, fill(design, 3), placements;
            connections = [(core = 1,), (core = 2,), (core = 3,)],
            id = "touching_row", line_length = 1.0)
        problem = LineParametersProblem(system, reference.earth;
            temperature = reference.temperature,
            frequencies = reference.frequencies[1:10:end])
        gmsh.initialize()
        try
            geometry = build(Geometry, FieldModel(problem, FEM()), options)
            @test isempty(geometry.holes.air)
            @test isempty(geometry.pockets.air)
            @test length(geometry.holes.earth) == 1
            @test isempty(geometry.pockets.earth)
            mktempdir() do directory
                path = joinpath(directory, "row.geo")
                write(path, geometry)
                lines = readlines(path)
                for medium in ("Air", "Earth")
                    line = only(filter(
                        line -> startswith(line, "FEM" * medium *
                                                 "Pockets()"), lines))
                    @test isempty(match(r"\{([^}]*)\}", line)[1])
                end
            end
        finally
            gmsh.model.remove()
            gmsh.finalize()
        end
    end
end
