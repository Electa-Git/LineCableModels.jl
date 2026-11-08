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
    @testset "the mesh parameters of round, annular and strip sections" begin
        dimension = LineCableModels.BaseGeometry.dimension
        controls = merge(options.overrides, (ConductorThicknessElements = 3,))
        for (x, radius) in ((LineCableModels.Ellipse(3.0e-3, 2.0e-3), 2.0e-3^2 / 3.0e-3),
            (LineCableModels.Ellipse(9.0e-3, 3.0e-3), 1.0e-3),
            (LineCableModels.Ellipse(3.0e-3, 3.0e-3), 3.0e-3), (Circle(2.0e-3), 2.0e-3))
            for (strand, target) in ((false, radius / 5), (true, radius))
                s = extension.section(x, x, strand, false, controls)
                @test s.target ≈ target atol = 1e-15
                @test s.width ≈ radius atol = 1e-15
                @test (s.divisions, s.fraction, s.depth) == (5, 0.5, 0)
            end
        end
        annulus = Annulus(1.0e-3, 3.0e-3)
        for layer in (false, true)
            s = extension.section(annulus, annulus, false, layer, controls)
            @test s.target ≈ 1.0e-3 rtol = 1e-12
            @test s.width ≈ 2.0e-3 rtol = 1e-12
            @test (s.divisions, s.fraction, s.depth) == (3, 0.45, 0)
        end
        # The insulated wire's layer, 1 mm thick at the radius 43.5 mm, region 2 of its model.
        insulated = FieldModel(Dict(REFERENCES)["two_insulated_wires"], FEM())
        for (patches, τ, target) in ((0, 1e-3, 0.5e-3), (1, 1e-3, 2π * 43.5e-3 / 96),
            (1, 1e-2, 2π * 43.5e-3 / 48))
            override = (ThinLayerStrips = patches, ConductorGeometryTolerance = τ)
            s = extension.section(insulated, 2, merge(controls, override))
            @test s.target ≈ target rtol = 1e-12
            @test s.width ≈ 1.0e-3 rtol = 1e-12
            @test (s.divisions, s.fraction, s.depth) == (3, 0.45, 0)
        end
        band = LineCableModels.resolve(Rectangle(4.0e-3, 2.0e-3), Offset(0.5e-3))
        fill = LineCableModels.PerforatedShape(band,
            LineCableModels.AbstractShape{Float64}[Circle(0.2e-3, Placement(0, 1.25e-3))])
        for shape in (band, fill)
            s = extension.section(Offset(0.5e-3), shape, false, false, controls)
            @test s.target ≈ 0.25e-3 rtol = 1e-12
            @test s.width ≈ 0.5e-3 rtol = 1e-12
            @test (s.divisions, s.fraction, s.depth) == (3, 0.45, 0)
        end
        # A strip's target reads its declared sides, and its width the `dimension` of its
        # resolved shape. The band's outer is 1 mm across, the band 0.5 mm thick. At the
        # placement 1.2 m away, the rectangle's `dimension` is 2.9e-15 m off its 2 mm.
        sector = LineCableModels.Sector(2.0e-3, 5.0e-3, π / 3, 0, 0, Placement(0, 0))
        polygon = LineCableModels.Polygon([
            (0.0, 0.0), (0.004, 0.0), (0.005, 0.003), (0.002, 0.005), (-0.001, 0.003)])
        square = LineCableModels.Polygon(
            [(-5.0e-3, -5.0e-3), (5.0e-3, -5.0e-3), (5.0e-3, 5.0e-3), (-5.0e-3, 5.0e-3)])
        for (x, side, width) in ((Rectangle(4.0e-3, 2.0e-3), 2.0e-3, 2.0e-3),
            (sector, 3.0e-3, dimension(sector)), (square, 10.0e-3, 10.0e-3),
            (polygon, dimension(polygon), dimension(polygon)), (band, 1.0e-3, 0.5e-3))
            for at in (Placement(0, 0), Placement(0, 10.0), Placement(0.3, -1.2, 0.7))
                shape = LineCableModels.resolve(at, x)
                for (strand, n) in ((false, 5), (true, 2))
                    s = extension.section(x, shape, strand, false, controls)
                    @test s.target ≈ side / n rtol = 1e-12
                    @test s.width === dimension(shape)
                    @test s.width ≈ width atol = 1e-14
                    @test s.depth === 0.45 * dimension(shape)
                    @test (s.divisions, s.fraction) == (5, 0.45)
                end
            end
        end
    end

    @testset "the model data carries each region's mesh parameters" begin
        model = FieldModel(Dict(REFERENCES)["two_insulated_wires"], FEM())
        mktempdir() do directory
            path = write(joinpath(directory, "model_data.pro"), model, options)
            lines = readlines(path)
            for i in eachindex(model.problem.system.regions)
                s = extension.section(model, i, options.overrides)
                for (name, number) in
                    (("RegionTarget", s.target), ("RegionWidth", s.width),
                    ("RegionDivisions", s.divisions), ("RegionFraction", s.fraction))
                    line = only(filter(startswith("$name~{$(i - 1)} = "), lines))
                    @test parse(Float64, chop(last(split(line, " = ")))) === Float64(number)
                end
            end
            @test !any(contains(r"^Region(Dimension|Repeated)"), lines)
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

    @testset "a hole 1 µm inside a circle lies in its face" begin
        reference = last(first(REFERENCES))
        model = FieldModel(reference, FEM())
        gmsh.initialize()
        try
            gmsh.model.add("walk")
            geometry = Geometry(model, options)
            # The hole passes 1 µm inside the circle at the angle 17π/64.
            at = Placement(8.999e-3 * cos(17π / 64), 8.999e-3 * sin(17π / 64))
            circle = draw(geometry, Contour(Circle(10e-3)).outer)
            hole = draw(geometry, Contour(Circle(1e-3, at)).outer)
            apart = draw(geometry, Contour(Circle(1e-3, Placement(20e-3, 0.0))).outer)
            surfaces, outside = extension.walk!(geometry, (circle.curves, -hole.curves))
            @test length(surfaces) == 1
            @test isempty(outside)
            surfaces, outside = extension.walk!(geometry, (circle.curves, -apart.curves))
            @test length(surfaces) == 1
            @test Set(only(outside)) == Set(-apart.curves)
        finally
            gmsh.model.remove()
            gmsh.finalize()
        end
    end

    @testset "an elliptic arc is drawn once and reused in either sense" begin
        reference = last(first(REFERENCES))
        model = FieldModel(reference, FEM())
        gmsh.initialize()
        try
            gmsh.model.add("elliptic")
            geometry = Geometry(model, options)
            for offset in (0.0, 0.5e-3)
                arc = LineCableModels.BaseGeometry.EllipticArc(
                    Placement(0.1, 0.2, π / 6), 3.0e-3, 2.0e-3, π / 8, 3π / 4, offset)
                reversed = LineCableModels.BaseGeometry.EllipticArc(
                    arc.at, arc.a, arc.b, arc.φ0 + arc.span, -arc.span, offset)
                drawn = draw(geometry, arc)
                stored = length(geometry.curves)
                @test draw(geometry, reversed) == -reverse(drawn)
                @test draw(geometry, arc) == drawn
                @test length(geometry.curves) == stored
            end
        finally
            gmsh.model.remove()
            gmsh.finalize()
        end
    end

    @testset "the chord of an elliptic arc of equal axes is the circular arc's" begin
        reference = last(first(REFERENCES))
        model = FieldModel(reference, FEM())
        chord = LineCableModels.BaseGeometry.direction
        gmsh.initialize()
        try
            for offset in (0.0, 0.4e-3), (φ0, span) in ((π / 8, π / 4), (3π / 8, -π / 4))

                gmsh.model.add("chord")
                geometry = Geometry(model, options)
                elliptic = only(draw(geometry,
                    LineCableModels.BaseGeometry.EllipticArc(
                        Placement(0.1, 0.2), 3.0e-3, 3.0e-3, φ0, span, offset)))
                radius = 3.0e-3 + offset
                circular = only(draw(geometry,
                    LineCableModels.BaseGeometry.Arc(Placement(0.1, 0.2, φ0), radius, span)))
                for sense in (1, -1), ds in (1e-5, 1e-4)

                    @test chord(geometry, sense * elliptic, ds) ≈
                          chord(geometry, sense * circular, ds) atol = 1e-12
                end
                gmsh.model.remove()
            end
        finally
            gmsh.finalize()
        end
    end

    @testset "the boundary walk turns at the elliptic arcs of the D and the duct" begin
        reference = last(first(REFERENCES))
        model = FieldModel(reference, FEM())
        D = Contour(
            LineCableModels.BaseGeometry.Curve{Float64}[
            LineCableModels.BaseGeometry.EllipticArc(Placement(0, 0), 20e-3, 5e-3, 0, π),
            LineCableModels.BaseGeometry.Segment((-20e-3, 0.0), (-20e-3, -1e-3)),
            LineCableModels.BaseGeometry.Segment((-20e-3, -1e-3), (20e-3, -1e-3)),
            LineCableModels.BaseGeometry.Segment((20e-3, -1e-3), (20e-3, 0.0))],
            Vector{LineCableModels.BaseGeometry.Curve{Float64}}[])
        gmsh.initialize()
        try
            gmsh.model.add("D")
            geometry = Geometry(model, options)
            outline = draw(geometry, D.outer)
            surfaces, outside = extension.walk!(geometry, (outline.curves,))
            @test length(surfaces) == 1
            @test isempty(outside)
            gmsh.model.remove()
            gmsh.model.add("duct")
            geometry = Geometry(model, options)
            duct = draw(geometry, Contour(LineCableModels.Ellipse(30e-3, 20e-3)).outer)
            core = draw(geometry, Contour(Circle(10e-3)).outer)
            surfaces, outside = extension.walk!(geometry, (duct.curves, -core.curves))
            @test length(surfaces) == 1
            @test isempty(outside)
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
        arc = LineCableModels.BaseGeometry.EllipticArc(
            Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, 0.5e-3)
        @test LineCableModels.BaseGeometry.turn(arc).fraction ≈ 1 / 4 rtol = 1e-12
        @test LineCableModels.BaseGeometry.turn(arc).length ≈
              0.003 * 1.3221199657742158 + 0.5e-3 * π / 2 rtol = 1e-12
        @test LineCableModels.BaseGeometry.turn(LineCableModels.BaseGeometry.EllipticArc(
            Placement(0, 0), 3.0e-3, 2.0e-3, 0, π / 2, 0)).length ≈
              0.003 * 1.3221199657742158 rtol = 1e-12
        problem = LineParametersProblem(system, homogeneous(rho = 100.0, eps_r = 10.0);
            frequencies = [0.1])
        wall = extension.section(FieldModel(problem, FEM()), length(system.regions),
            options.overrides)
        @test wall.target ≈ 1.5e-3 rtol = 1e-12
        @test wall.width ≈ 3e-3 rtol = 1e-12
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
                @test geometry.curves[abs(tag)] isa LineCableModels.BaseGeometry.EllipticArc
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
        model = FieldModel(problem, FEM())
        strand = extension.section(model, 8, options.overrides)
        @test (strand.target, strand.width, strand.depth) == (0.001, 0.002, 0.45 * 0.002)
        container = extension.section(model, length(system.regions), options.overrides)
        @test container.target ≈ 0.002 rtol = 1e-12
        @test container.width ≈ 0.004 rtol = 1e-12
        gmsh.initialize()
        try
            geometry = build(Geometry, FieldModel(problem, FEM()), options)
            @test length.(geometry.surfaces) == [fill(1, 7); fill(5, 3); 9]
            @test length(gmsh.model.getEntities(2)) == 31
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
