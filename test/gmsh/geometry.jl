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
    @testset "the circle rule counts segments" begin
        for (τ, count) in ((1e-3, 96), (1e-2, 48), (1e-4, 384), (2.0, 12))
            @test extension.segments(τ) == count
        end
    end
    @testset "section records of an ellipse polygon and contour offset" begin
        ellipse = LineCableModels.Ellipse(3.0e-3, 2.0e-3)
        polygon = LineCableModels.Polygon([
            (0.0, 0.0), (0.004, 0.0), (0.005, 0.003), (0.002, 0.005), (-0.001, 0.003)])
        offset = LineCableModels.resolve(Rectangle(4.0e-3, 2.0e-3), Offset(0.5e-3))
        @test extension.section(ellipse) == (kind = 1, a = 2.0e-3, b = 0.0)
        @test extension.section(polygon).kind == 3
        @test collect(extension.section(polygon))[2:3] ≈ [0.006, 0.005] rtol = 1e-12
        @test extension.section(offset).kind == 3
        @test collect(extension.section(offset))[2:3] ≈ [0.005, 0.003] rtol = 1e-12
        @test extension.section(Offset(0.5e-3), offset).kind == 2
        @test collect(extension.section(Offset(0.5e-3), offset))[2:3] ≈ [0.5e-3, 0] rtol = 1e-12
        fill = LineCableModels.PerforatedShape(offset,
            LineCableModels.AbstractShape{Float64}[Circle(0.2e-3, Placement(0, 1.25e-3))])
        @test extension.section(Offset(0.5e-3), fill).kind == 2
        @test collect(extension.section(Offset(0.5e-3), fill))[2:3] ≈ [0.5e-3, 0] rtol = 1e-12
    end

    # The record of a band takes its thickness from the band itself, the same at every
    # placement of the band.
    @testset "a band's record is its thickness at every placement" begin
        band = LineCableModels.resolve(Rectangle(10.0e-3, 4.0e-3), Offset(1.0e-3))
        for place in (Placement(0, 0), Placement(0, 10.0), Placement(0.3, -1.2, 0.7))
            record = extension.section(Offset(1.0e-3), LineCableModels.resolve(place, band))
            @test record.kind == 2
            @test record.a - record.b ≈ 1.0e-3 rtol = 1e-12
        end
    end

    @testset "rectangle vertices split neighboring edges" begin
        reference = last(first(REFERENCES))
        model = FieldModel(reference, FEM())
        gmsh.initialize()
        try
            gmsh.model.add("rectangles")
            geometry = Geometry(model, options)
            @test geometry.segments == 96
            # Draw the smaller neighbor first to register its vertices on the edge.
            neighbor = Contour(
                Curve{Float64}[
                Segment((2.0, 0.5), (3.0, 0.5)),
                Segment((3.0, 0.5), (3.0, 1.5)),
                Segment((3.0, 1.5), (2.0, 1.5)),
                Segment((2.0, 1.5), (2.0, 0.5))],
                Vector{Curve{Float64}}[])
            rectangle = Contour(
                Curve{Float64}[
                Segment((0.0, 0.0), (2.0, 0.0)),
                Segment((2.0, 0.0), (2.0, 2.0)),
                Segment((2.0, 2.0), (0.0, 2.0)),
                Segment((0.0, 2.0), (0.0, 0.0))],
                Vector{Curve{Float64}}[])
            neighbor = draw(geometry, neighbor.outer)
            rectangle = draw(geometry, rectangle.outer)
            @test length(rectangle.curves) == 6
            @test length(intersect(abs.(rectangle.curves), abs.(neighbor.curves))) == 1
            @test Set(minmax(geometry.positions[p], geometry.positions[q])
            for (p, q) in keys(geometry.lines)
            if geometry.positions[p][1] == geometry.positions[q][1] == 2.0) ==
                  Set([((2.0, 0.0), (2.0, 0.5)),
                ((2.0, 0.5), (2.0, 1.5)), ((2.0, 1.5), (2.0, 2.0))])
        finally
            gmsh.model.remove()
            gmsh.finalize()
        end
    end

    @testset "the parallel of an elliptical duct is drawn within the tolerance" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        air = Material(:insulator, Inf)
        pvc = Material(:insulator, Inf, 3)
        design = build(CableDesign, "duct_wall",
            Stack(
                Enclosure(Region(:duct, LineCableModels.Ellipse(30e-3, 20e-3), air),
                    Terminal(:core, Region(:core, Circle(0.01), copper))),
                Region(:wall, Offset(3e-3), pvc)))
        system = build(CableSystem, [design], [Placement(0, -1)];
            connections = [(core = 1,)], id = "duct_wall", line_length = 1.0)
        record = extension.section(last(system.regions))
        @test record.kind == 2
        @test record.a - record.b ≈ 3e-3 rtol = 1e-12
        arc = LineCableModels.DataModel.EllipticArc(
            Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, 0.5e-3)
        @test LineCableModels.DataModel.turn(arc).fraction ≈ 1 / 4 rtol = 1e-12
        @test LineCableModels.DataModel.turn(arc).length ≈
              0.003 * 1.3221199657742158 + 0.5e-3 * π / 2 rtol = 1e-12
        @test LineCableModels.DataModel.turn(LineCableModels.DataModel.EllipticArc(
            Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, 0)).length ≈
              0.003 * 1.3221199657742158 rtol = 1e-12
        problem = LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
            frequencies = [0.1])
        gmsh.initialize()
        try
            geometry = build(Geometry, FieldModel(problem, FEM()), options)
            gmsh.model.geo.synchronize()
            points = Tuple{Float64, Float64}[]
            for tag in geometry.outlines[1].curves
                low, high = gmsh.model.getParametrizationBounds(1, abs(tag))
                parameters = collect(range(only(low), only(high); length = 4096))
                tag < 0 && reverse!(parameters)
                for parameter in parameters
                    xyz = gmsh.model.getValue(1, abs(tag), [parameter])
                    push!(points, (xyz[1], xyz[2]))
                end
                @test geometry.curves[abs(tag)] isa LineCableModels.DataModel.EllipticArc
                @test geometry.curves[abs(tag)].offset == 3e-3
            end
            A_d = sum(points[k][1] * points[mod1(k + 1, length(points))][2] -
                      points[k][2] * points[mod1(k + 1, length(points))][1]
            for k in eachindex(points)) / 2
            @test abs((A_d - π * 30e-3 * 20e-3) /
                      (0.12 * 1.3221199657742158 * 3e-3 + π * 3e-3^2) - 1) <= 1e-3
            repeated = draw(geometry, Contour(LineCableModels.resolve(Placement(0, -1), design.outer)).outer)
            @test repeated.curves == geometry.outlines[1].curves
        finally
            gmsh.model.remove()
            gmsh.finalize()
        end
    end

    @testset "an annular fill around round and rectangular courses" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        polymer = Material(:insulator, Inf)
        design = build(CableDesign, "mixed_courses",
            Terminal(:core, Region(:core_metal, Circle(0.01), copper)),
            Terminal(:screen,
                Enclosure(Region(:fill, Offset(0.004), polymer),
                    Stack(
                        Repeat(CircularPattern(6; r = Offset(0)),
                            Region(:wires, Circle(0.001), copper)),
                        Repeat(CircularPattern(3; r = Offset(0)),
                            Region(:strands, Rectangle(0.004, 0.002), copper))))))
        system = build(CableSystem, [design], [Placement(0.0, -1.0)];
            connections = [(core = 1, screen = 2)], id = "mixed_courses", line_length = 1.0)
        problem = LineParametersProblem(system, homogeneous(rho = 100); frequencies = [50.0])
        plan = extension.compartments(last(system.regions).shape)
        @test (
            plan.center, plan.inner_radius, plan.outer_radius, plan.shell_outer_radius) ==
              ((0.0, -1.0), 0.01, 0.012, 0.014)
        @test plan.angles ≈ collect(0:5) .* (π / 3)
        @test length(plan.holes) == 6
        @test length(plan.outer_spans) == 3
        @test first.(plan.outer_spans) ≈ [2 / 13 + k * 2π / 3 for k in 0:2]
        @test last.(plan.outer_spans) ≈ fill(2π / 3 - 4 / 13, 3)
        @test extension.section(system.regions[8]) == (kind = 3, a = 0.004, b = 0.002)
        @test extension.section(last(system.regions)) == (kind = 2, a = 0.014, b = 0.01)
        sector = system.regions[8].shape
        @test collect(extension.lowest(Contour(sector))) ≈
              [0.014cos(2 / 13), -1.0 - 0.014sin(2 / 13)]
        gmsh.initialize()
        try
            geometry = build(Geometry, FieldModel(problem, FEM()), options)
            @test length.(geometry.surfaces) == [fill(1, 10); 9]
            @test length(gmsh.model.getEntities(2)) == 19
        finally
            gmsh.model.remove()
            gmsh.finalize()
        end
    end

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
