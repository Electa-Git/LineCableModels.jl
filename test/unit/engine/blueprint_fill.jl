@testitem "Blueprint / finite-resistivity interstitial fill" tags=[:unit] begin
    import LineCableModels.Engine as Engine
    import LineCableModels.DataModel as DM

    copper = Material(kind=:conductor, rho=1.7241e-8)
    pe = Material(kind=:insulator, rho=1.97e14, eps_r=2.3)
    design = @cable "finite-pe-fill" begin
        @terminal :core begin
            stranded(copper; shape=Disk(0.5e-3), boundary=Disk(2.5e-3), fill=pe)
            insulation(pe; t=1e-3)
        end
        @terminal :sheath begin
            sheath(copper; t=0.2e-3)
        end
    end
    fills = filter(r -> r.primitive isa DM.DifferenceShape, design.geometry.regions)
    @test length(fills) == 1
    @test only(fills).source.material == pe
    @test isfinite(only(fills).source.material.rho)

    blueprint = Engine.flatten(LineCableModelsCoaxial(), design, Float64)
    @test getproperty.(blueprint.conductors, :terminal) == [:core, :sheath]
    @test blueprint.assembly_ranges == [1:2]
    @test only(blueprint.dielectrics).material == pe
    @test only(blueprint.dielectrics).r_in ≈ 2.5e-3
    @test only(blueprint.dielectrics).r_ex ≈ 3.5e-3

    constants = CableConstants(design; frequency=50.0)
    row = only(constants)
    @test all(isfinite, (row.R, row.L, row.C, row.G))
    @test row.G > 0
    @test row.G ≈ 2pi / (pe.rho * log(3.5 / 2.5))

    system = @system "finite-pe-fill" begin
        @at design (0.0, -1.0) core=1 sheath=2
    end
    problem = LineParametersProblem(system; earth_props=homogeneous(rho=100.0), frequencies=[50.0])
    parameters = compute(problem)
    @test size(parameters.Z) == (2, 2, 1)
    @test size(parameters.Y) == (2, 2, 1)
    @test all(isfinite, parameters.Z)
    @test all(isfinite, parameters.Y)
end
