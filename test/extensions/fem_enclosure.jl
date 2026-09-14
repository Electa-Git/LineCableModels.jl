@testitem "Gmsh FEM / complete strand envelopes preserve enclosing material" tags=[:extension] begin
    using Gmsh
    const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    const DM = LineCableModels.DataModel
    copper = Material(kind=:conductor, rho=1.72e-8)
    dielectric = Material(kind=:insulator, rho=1e10, eps_r=3.0)
    wire_radius = 0.5e-3
    envelope = Disk(sqrt(7) * wire_radius)
    exterior = Ellipse(4e-3, 3e-3)
    compacted = stranded(copper; shape=Disk(wire_radius), compact=true, boundary=envelope)
    designs = map((compacted, solid(copper, envelope))) do body
        build(CableDesign, "bounded-enclosure",
            pipe(terminal(:core, body); shape=exterior, fill=dielectric))
    end
    source_fill = last(first(designs).geometry.regions).primitive
    @test source_fill isa DM.DifferenceShape
    @test length(source_fill.holes) == 7
    for dielectric_formula in (:default, :Ametani2004)
        formulation = Formulation(:LineCableModelsFEM;
            insulation_admittance=dielectric_formula)
        models = map(designs) do design
            system = build(LineCableSystem, design, Pose2(0.0, -0.1);
                connections=Dict(:core=>1))
            problem = LineParametersProblem(system; frequencies=[50.0, 1000.0],
                earth_props=homogeneous(rho=100.0))
            FEM._resolved_fem_model(problem, formulation)
        end
        for model in models
            @test length(model.region_plans) == 2
            metal, matrix = model.region_plans
            @test metal.terminal_index == 1
            @test matrix.terminal_index == 0
            @test metal.shape isa Disk
            @test matrix.shape isa DM.DifferenceShape
            @test only(matrix.shape.holes) == metal.shape
            @test area(metal.shape) ≈ 7pi * wire_radius^2
            @test sum(plan -> area(plan.shape), model.region_plans) ≈ area(exterior)
            coefficients = model.material_plans[matrix.material_index].admittivity
            @test real.(coefficients) == fill(dielectric_formula === :default ? 0.0 : 1e-10, 2)
            @test imag.(coefficients) ≈ 2pi .* [50.0, 1000.0] .* 8.8541878128e-12 .* 3
        end
        @test [plan.admittivity for plan in first(models).material_plans] ==
            [plan.admittivity for plan in last(models).material_plans]
        @test area.(getproperty.(first(models).region_plans, :shape)) ≈
            area.(getproperty.(last(models).region_plans, :shape))
    end
    # Coalescing is a detached FEM adaptation, not a mutation of strand identity.
    @test length(first(designs).geometry.regions) == 8
    @test length(source_fill.holes) == 7
end

@testitem "Gmsh FEM / disjoint assembly boundaries retain terminal and mesh ownership" tags=[:extension] begin
    using Gmsh
    const FEM = Base.get_extension(LineCableModels, :LineCableModelsGmshExt)
    const DM = LineCableModels.DataModel
    copper = Material(kind=:conductor, rho=1.72e-8)
    design = build(CableDesign, "disjoint-assembly", assembly(
        at(terminal(:a, core(copper; r=0.5e-3)), -2e-3, 0),
        at(terminal(:b, core(copper; r=0.8e-3)), 2e-3, 0),
    ))
    for height in (-1.0, 1.0)
        system = build(LineCableSystem, design, Pose2(0.0, height);
            connections=Dict(:a=>1, :b=>2))
        problem = LineParametersProblem(system; frequencies=[50.0],
            earth_props=homogeneous(rho=100.0))
        model = FEM._resolved_fem_model(problem, LineCableModelsFEM())
        @test only(model.cable_boundaries) isa DM.AssemblyShape
        @test length(only(model.cable_boundaries).members) == 2
        @test model.cable_hosts == [height > 0 ? :air : :earth]
        @test getproperty.(model.region_plans, :terminal_index) == [1, 2]
        @test only(model.cable_outer_mesh_sizes) == minimum(plan -> plan.mesh_size, model.region_plans)
        @test sum(plan -> area(plan.shape), model.region_plans) ≈ pi * (0.5e-3^2 + 0.8e-3^2)
        boundary_shape = only(model.cable_boundaries)
        @test boundary(boundary_shape) === boundary_shape
        @test area(boundary_shape) ≈ pi * (0.5e-3^2 + 0.8e-3^2)
        @test perimeter(boundary_shape) ≈ 2pi * (0.5e-3 + 0.8e-3)
        centre = centroid(boundary_shape)
        @test centre[1] ≈ 2e-3 * (0.8e-3^2 - 0.5e-3^2) / (0.8e-3^2 + 0.5e-3^2)
        @test centre[2] ≈ height
        for angle in (0.0, pi / 4, pi / 2, pi)
            @test support(boundary_shape, angle) ≈ max(
                -2e-3 * cos(angle) + height * sin(angle) + 0.5e-3,
                2e-3 * cos(angle) + height * sin(angle) + 0.8e-3)
        end
        @test support(boundary_shape) ≈ hypot(2e-3, height) + 0.8e-3
    end
end
