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
    design=PB.CableDesign(terminal; cable_id = "parametric")
    @test design isa Gridspace{CableDesign}
    @test length(design) == 2
    @test all(item -> item isa CableDesign, design)
    @test all(item -> item.terminal_order == [:core], design)

    system=PB.LineCableSystem(
        design;
        position = PB.at(x = 0.0, y = -1.0),
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

    designs=PB.CableDesign(roots; cable_id = "two-thicknesses")
    @test outer_radius.(collect(designs)) ≈ [0.013, 0.014]
    @test all(design -> design.effective !== nothing, designs)

    @test PB.at(x = 1, y = -2) == Pose2(1, -2, 0)
    @test length(PB.trifoil(y = -1, spacing = 0.1)) == 3
    @test length(PB.hflat(spacing = 0.1, n = 4)) == 4
    @test length(PB.vflat(spacing = 0.1, n = 2)) == 2
end

@testitem "Architecture / retired cable builders are absent" tags=[:quality] begin
    using LineCableModels
    import LineCableModels.DataModel as DM
    import LineCableModels.ParametricBuilder as PB

    for name in (
        :CableBuilder, :SystemBuilder, :PartBuilder, :PositionDefinition,
        :CableComponent, :CablePosition, :ConductorGroup, :InsulatorGroup,
        :CircStrands, :RectStrands, :Tubular, :Strip, :Semicon
    )
        @test !isdefined(LineCableModels, name)
        @test !isdefined(DM, name)
        @test !isdefined(PB, name)
    end

    retired=(
        "CableBuilder", "SystemBuilder", "PartBuilder", "GroupBuilder",
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
