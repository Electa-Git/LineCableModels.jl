@testitem "UQ / ordinary joint rings stacks and bounded cores" tags=[:integration, :extension] begin
    using Measurements, Random, Statistics
    const DM = LineCableModels.DataModel
    b0, d0, count = 0.0679,0.006,68
    gap(b,d) = 2(b+d/2)*sinpi(1/count)-d
    @test gap(b0,d0) > 0
    error = try
        DM.placements(Ring(count;r=b0+1.1d0/2),Disk(1.1d0/2),nothing)
        nothing
    catch exception
        exception
    end
    @test error isa DomainError
    @test error.val == count
    @test occursin("deficit=",error.msg)
    @test occursin("available chord=",error.msg)
    copper = Material(kind=:conductor,rho=1.72e-8)
    steel = Material(kind=:conductor,rho=2e-7,mu_r=10.)
    dielectric = Material(kind=:insulator,rho=Inf,eps_r=2.3)
    function problem(p,bounded)
        b,d = b0*p.scale,d0*p.scale
        core = bounded ? stranded(copper;shape=Disk(0.01p.scale),
            boundary=Disk(0.04p.scale),compact=true,fill=dielectric) :
            Region(:metal,Disk(0.04p.scale),copper)
        ring = Group(:armour,Region(:wires,Disk(d/2),steel);
            pattern=Ring(count;r=b+d/2))
        cable = build(CableDesign,"joint",Stack(
            terminal(:core,core),Region(:insulation,Shell(b-0.04p.scale),dielectric),
            Enclosure(:interstices,ring;primitive=Annulus(b,b+d),fill=dielectric),
            Region(:jacket,Shell(0.003p.scale),dielectric)))
        system = build(LineCableSystem,[cable],[Pose2(p.x,-1.)];
            connections=[Dict(:core=>1,:armour=>2)],system_id="joint")
        return LineParametersProblem(system;frequencies=[0.1,50.,1e7],
            earth_props=homogeneous(rho=100.,eps_r=10.))
    end
    inner = Formulation(options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    # One scale simultaneously controls the ring, the stack and its local
    # dimensions. A separate bounded coordinate remains an independent input.
    inputs = Gridspace{NamedTuple{(:scale,:x)}}((base,s,x)->(scale=base*s,x=x),
        (Grid((1.,1.05)),Grid(1.,10.),Grid(0.,AbsoluteError(0.002))))
    for bounded in (false,true)
        space = Gridspace{LineParametersProblem}(p->problem(p,bounded),(inputs,))
        parametric = ParametricProblem(space)
        lep = compute(parametric,LinearError(inner))
        mc = compute(parametric,MonteCarlo(inner;trials=16,seed=0x1234,
            distribution=:uniform,return_samples=true,retain_details=true))
        @test mc.trial_counts == [16,16]
        @test all(isempty,mc.details.failures)
        @test length(lep.values) == 2
        @test all(v->all(isfinite,observe(v,R)),lep.values)
        @test all(v->all(isfinite,observe(v,R)),mc.values)
        z = measurement(1.,0.01)
        measured = compute(problem((scale=z,x=0.),bounded),inner)
        @test any(v->uncertainty(v)>0,observe(measured,R))
        for h in (1e-5,3e-6)
            plus = compute(problem((scale=1+h,x=0.),bounded),inner)
            minus = compute(problem((scale=1-h,x=0.),bounded),inner)
            for quantity in (R,L,C,G)
                expected = (observe(plus,quantity).-observe(minus,quantity))./(2h)
                actual = map(observe(measured,quantity)) do v
                    v isa Measurement ? Measurements.derivative(v,z) : 0.
                end
                @test actual ≈ expected rtol=2e-4 atol=1e-12
            end
        end
        sampled = rand(Xoshiro(3),space;distribution=:uniform)
        @test length(first(sampled.system.designs).geometry.regions) ==
            length(first(first(space).system.designs).geometry.regions)
        @test uncertainty(first(first(space).system.positions).x) ≈ 0.002
    end
end
