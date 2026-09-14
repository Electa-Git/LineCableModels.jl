@testitem "Gauntlet / current touching disks / UQ clearance and retained geometry" tags=[:gauntlet_toolkit] setup=[GauntletSupport] begin
    using LineCableModels,Measurements,LinearAlgebra
    built=LineCableSystem[]
    function problem(radius)
        radius>0 || throw(DomainError(radius,"radius must be positive"))
        wire=build(CableDesign,"touching-current-disk",terminal(:core,
            Region(:metal,Disk(radius),Material(kind=:conductor,rho=2e-8))))
        system=build(LineCableSystem,[wire,wire],[(0.,-1.),(.02,-1.)];
            connections=[Dict(:core=>1),Dict(:core=>2)])
        push!(built,system)
        LineParametersProblem(system;frequencies=[50.,1000.],earth_props=homogeneous(rho=100.))
    end
    space=Gridspace{LineParametersProblem}(problem,(Grid(.01,AbsoluteError(.0002)),))
    inner=Formulation(options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    sampled=compute(ParametricProblem(space),MonteCarlo(inner;trials=8,seed=103,
        distribution=:uniform,return_samples=true,retain_details=true))
    @test sampled.trial_counts==[8]
    @test isempty(only(sampled.details.failures))
    @test length(built)==9
    @test only(sampled.details.clearance).adjustments>0
    for system in built
        centres=[centroid(resolve(pose,boundary(design.geometry))) for (pose,design) in zip(system.positions,system.designs)]
        radii=[outer_radius(design) for design in system.designs]
        gap=norm(collect(centres[1]).-collect(centres[2]))-sum(radii)
        @test nominal(gap)>=nominal(system.clearances[1,2])-64eps(Float64)
        @test nominal(gap)>0
    end
    linear=compute(ParametricProblem(space),LinearError(inner))
    @test all(isfinite,nominal.(Z(only(linear.values))))
    @test all(isfinite,uncertainty.(real.(Y(only(linear.values)))))
end
