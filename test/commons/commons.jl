@testset "Commons" begin
    @testset "validate returns what it checks" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        @test LineCableModels.validate(copper) === copper
        earth = homogeneous(rho = 100, eps_r = 10)
        @test LineCableModels.validate(earth) === earth
    end

    @testset "resolve gives a layer its place around a part" begin
        @test LineCableModels.resolve(Circle(1.0e-3), Shell(1.0e-3)) ==
              Annulus(1.0e-3, 2.0e-3)
    end

    @testset "build gives a cable its regions and terminals" begin
        copper = Material(:conductor, 1.7241e-8, 1, 1, 20, 0.00393)
        wire = build(CableDesign, "wire",
            Terminal(:core, Region(:core_metal, Circle(1.0e-3), copper)))
        @test [(row.name, row.regions) for row in wire.terminals] == [(:core, [1])]
    end
end
