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
        finally
            gmsh.model.remove()
            gmsh.finalize()
        end
    end
end
