@testitem "ParametricBuilder / armor grids retain course area and helical resistance" tags=[:unit, :parametric] begin
    const DM = LineCableModels.DataModel
    copper = Material(kind=:conductor, rho=1.72e-8)
    steel = Material(kind=:conductor, rho=1.5e-7, mu_r=80.0)
    dielectric = Material(kind=:insulator, rho=Inf, eps_r=2.3)
    path = Helix(LayRatio(10); dir=-1, φ0=0.2)
    course = armor(steel; shape=Disk(0.1e-3), n=Grid((12, 18)), lay=path)
    @test course isa Gridspace{Group}
    designs = build(CableDesign, "armor-grid", Stack(
        terminal(:core, core(copper; r=1e-3), insulation(dielectric; t=0.5e-3)),
        course,
        shell(dielectric; t=0.2e-3),
    ))
    @test designs isa Gridspace{CableDesign}
    @test length(designs) == 2
    for (count, design) in zip((12, 18), designs)
        wires = filter(placed -> placed.terminal === :armor, design.geometry.regions)
        @test length(wires) == count
        @test all(placed -> placed.primitive isa Disk, wires)
        @test all(placed -> area(placed.primitive) ≈ pi * (0.1e-3)^2, wires)
        @test all(placed -> only(placed.paths).path === path, wires)
        @test all(placed -> only(placed.paths).radius ≈ 1.6e-3, wires)
        @test design.terminal_order == [:core, :armor]
        @test outer_radius(design) ≈ 1.9e-3
        blueprint = LineCableModels.DataModel.flatten(design)
        conductor = only(filter(row -> row.terminal === :armor, blueprint.conductors))
        @test sum(placed -> area(placed.primitive), wires) ≈ count * pi * (0.1e-3)^2
        @test DM.tubular_resistance(conductor.r_in, conductor.r_ex, conductor.material.rho) ≈
            steel.rho * sqrt(1 + (pi / 10)^2) / (count * pi * (0.1e-3)^2) rtol = 1e-12
        @test all(placed -> placed.source.material === steel, wires)
    end
    @test_throws ArgumentError armor(dielectric; shape=Disk(0.1e-3), n=12)
end
