@testitem "ParametricBuilder / armor grids retain course area and helical resistance" tags=[:unit] begin
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
        wires = filter(region -> region.terminal === :armor, design.geometry.regions)
        @test length(wires) == count
        @test all(region -> region.primitive isa Disk, wires)
        @test all(region -> area(region.primitive) ≈ pi * (0.1e-3)^2, wires)
        @test all(region -> only(region.paths).path === path, wires)
        @test all(region -> only(region.paths).radius ≈ 1.6e-3, wires)
        @test design.terminal_order == [:core, :armor]
        @test outer_radius(design) ≈ 1.9e-3
        layers = DM.flatten(design, 50.0)
        conductor = only(filter(component -> component.name === :armor, layers)).conductor
        @test conductor.cross_section ≈ count * pi * (0.1e-3)^2
        @test conductor.resistance ≈
            steel.rho * sqrt(1 + (pi / 10)^2) / (count * pi * (0.1e-3)^2)
        @test all(region -> region.source.material === steel, wires)
    end
    @test_throws ArgumentError armor(dielectric; shape=Disk(0.1e-3), n=12)
end
