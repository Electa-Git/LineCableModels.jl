@testitem "ParametricBuilder / v1 direct target materialization" tags=[:unit] setup=[
    EngineTestSupport,
    UseEngineSupport
] begin
    import LineCableModels.ParametricBuilder as PB

    conductor=PB.Material(
        kind = Grid((:conductor,)),
        rho = Grid((1.7e-8, 2.8e-8)),
        eps_r = 1.0,
        mu_r = 1.0
    )
    @test conductor isa Gridspace{Material}
    @test length(conductor) == 2
    @test all(material -> material.kind === :conductor, conductor)

    wire=PB.Conductor.Solid(:wire, conductor; r = 0.01)
    terminal=PB.Group(:core, wire)
    design=build(CableDesign, "parametric", terminal)
    @test design isa Gridspace{CableDesign}
    @test length(design) == 2
    @test all(item -> item isa CableDesign, design)
    @test all(item -> item.terminal_order == [:core], design)

    system=build(
        LineCableSystem,
        design,
        PB.at(x = 0.0, y = -1.0);
        connections = Dict(:core=>1),
        system_id = "parametric-system",
        line_length = 100.0
    )
    @test system isa Gridspace{LineCableSystem}
    @test length(system) == length(design)
    @test all(item -> item isa LineCableSystem, system)

    problems=LineParametersProblem(
        system;
        earth_props = PB.Earth(rho = 100.0, eps_r = 10.0),
        frequencies = [50.0]
    )
    @test problems isa Gridspace{LineParametersProblem}
    @test length(problems) == 2
    @test all(problem -> problem.system isa LineCableSystem, problems)
    @test (@inferred CableDesign first(design)) isa CableDesign
    @test (@inferred LineCableSystem first(system)) isa LineCableSystem

    identifiers=Grid(("first", "second"))
    direct=build(CableDesign, "first", first(terminal))
    design_space=build(CableDesign, identifiers, first(terminal))
    @test eltype(design_space) === CableDesign
    @test first(design_space).root == direct.root
    @test first(design_space).terminal_order == direct.terminal_order
    @test rand(design_space) isa CableDesign

    atomic=build(
        LineCableSystem,
        [first(design), first(design)],
        [Pose2(-0.1, -1.0), Pose2(0.1, -1.0)];
        connections = [(core = 1,), (core = 2,)]
    )
    @test atomic isa LineCableSystem
    @test ncables(atomic) == 2

    nested_identifier=Gridspace{String}(
        value->string("nested-", value),
        (Grid((1, 2)),)
    )
    nested_designs=build(CableDesign, nested_identifier, first(terminal))
    @test nested_designs isa Gridspace{CableDesign}
    @test first(nested_designs).cable_id == "nested-1"

    import LineCableModels.DataModel as DM
    const resolution_count=Ref(0)
    struct CountedDefinition{T <: Real}<:DM.AbstractPrimitiveDefinition{T}
        radius::T
    end
    struct CountedPrimitive{T <: Real}<:DM.AbstractPrimitive{T}
        radius::T
        at::DM.Pose2{T}
    end
    DM.resolve(::DM.EmptyBoundary, definition::CountedDefinition) = begin
        resolution_count[]+=1
        CountedPrimitive(definition.radius, DM.Pose2(0, 0, 0))
    end
    DM.resolve(at::DM.Pose2, primitive::CountedPrimitive) = CountedPrimitive(primitive.radius, at*primitive.at)
    DM.boundary(primitive::CountedPrimitive) = primitive
    DM.area(primitive::CountedPrimitive) = pi*primitive.radius^2
    DM.centroid(primitive::CountedPrimitive) = (primitive.at.x, primitive.at.y)
    DM.support(primitive::CountedPrimitive, φ::Real) = primitive.at.x*cos(φ)+primitive.at.y*sin(φ)+primitive.radius
    DM.support(primitive::CountedPrimitive) = hypot(primitive.at.x, primitive.at.y)+primitive.radius
    counted_root=Group(
        :core,
        Region(:counted, CountedDefinition(0.01), Material(
            kind = :conductor,
            rho = 1.7e-8
        ))
    )
    counted_space=build(CableDesign, Grid(("counted-a", "counted-b")), counted_root)
    counted_designs=collect(counted_space)
    @test resolution_count[] == 2
    LineParametersProblem(
        first(counted_designs),
        Pose2(0.0, -1.0);
        connections = (core = 1,),
        earth_props = PB.Earth(rho = 100.0),
        frequencies = [1.0, 10.0, 100.0]
    )
    @test resolution_count[] == 2
end

@testitem "ParametricBuilder / v1 physical conveniences preserve one grammar" tags=[:unit] begin
    using LineCableModels
    import LineCableModels.ParametricBuilder as PB

    copper=Material(kind = :conductor, rho = 1.7e-8)
    xlpe=Material(kind = :insulator, rho = 1.0e14, eps_r = 2.3)
    core=PB.Conductor.Solid(:core, copper; r = 0.01)
    insulation=PB.Insulator.Shell(:insulation, xlpe; t = Grid((0.003, 0.004)))
    roots=PB.Stack(PB.Group(:core, core), insulation)
    @test roots isa Gridspace{Stack}

    designs=build(CableDesign, "two-thicknesses", roots)
    @test outer_radius.(collect(designs)) ≈ [0.013, 0.014]
    @test all(design -> design.geometry isa CableGeometry, designs)

    @test PB.at(x = 1, y = -2) == Pose2(1, -2, 0)
    design = first(designs)
    connections = (core = (1, 2, 3),)
    @test length(PB.trefoil(
        design;
        center = PB.at(0, -1),
        spacing = 0.1,
        connections
    )) == 3
    @test length(PB.hflat(design; spacing = 0.1, connections)) == 3
    @test length(PB.vflat(design; spacing = 0.1, connections)) == 3
end

@testitem "Architecture / retired cable builders are absent" tags=[:quality] begin
    using LineCableModels
    import LineCableModels.DataModel as DM
    import LineCableModels.ParametricBuilder as PB

    for name in (
        :CableBuilder, :SystemBuilder, :PartBuilder, :PositionDefinition,
        :DesignMaterializer, :SystemMaterializer, :MemberPlacementSpec,
        :HeterogeneousAssembly, :AssemblySpec, :StrandDefinition,
        :RopeDefinition, :DuctBank, :ResolvedPart, :ResolvedRegion,
        :CableComponent, :CablePosition, :ConductorGroup, :InsulatorGroup,
        :CircStrands, :RectStrands, :Tubular, :Strip, :Semicon
    )
        @test !isdefined(LineCableModels, name)
        @test !isdefined(DM, name)
        @test !isdefined(PB, name)
    end

    retired=(
        "CableBuilder", "SystemBuilder", "PartBuilder", "GroupBuilder",
        "DesignMaterializer", "SystemMaterializer", "MemberPlacementSpec",
        "HeterogeneousAssembly", "AssemblySpec", "StrandDefinition",
        "RopeDefinition", "DuctBank", "ResolvedPart", "ResolvedRegion",
        "BuildContext", "BuildManager", "BuildPlan", "BuildSpec",
        "BuildResult", "CableConstantsProblem", "compute_cable_constants",
        "_current_radial_equivalence", "reference_frequency",
        "CableComponent", "CablePosition", "ConductorGroup", "InsulatorGroup",
        "CircStrands", "RectStrands", "ConductorPart", "InsulatorPart",
        "SolidCore", "TubularLayer", "EnclosureSpec", "StrandedSpec",
        "LayRatioSpec", "AbstractRole", "AbstractPattern", "AbstractPath",
        "AbstractCompaction", "AbstractNames", "PartGroup"
    )
    pattern=Regex("\\b(?:"*join(retired, "|")*")\\b")
    root=pkgdir(LineCableModels)
    maintained=String[
        joinpath(root, "README.md"),
        joinpath(root, "CHANGELOG.md")
    ]
    for directory in ("src", "examples", "dev", joinpath("docs", "src"))
        base=joinpath(root, directory)
        append!(maintained,
            [joinpath(path, file)
             for (path, _, files) in walkdir(base)
             for file in files
             if endswith(file, ".jl")||endswith(file, ".md")])
    end
    offenders=filter(maintained) do path
        occursin(pattern, read(path, String))
    end
    @test isempty(offenders)

    for path in (
        joinpath("src", "parametricbuilder", "cablebuilder.jl"),
        joinpath("src", "parametricbuilder", "systembuilder.jl"),
        joinpath("src", "parametricbuilder", "parts.jl"),
        joinpath("src", "datamodel", "cablecomponent"),
        joinpath("src", "datamodel", "conductorgroup.jl"),
        joinpath("src", "datamodel", "insulatorgroup.jl")
    )
        @test !ispath(joinpath(root, path))
    end
end
