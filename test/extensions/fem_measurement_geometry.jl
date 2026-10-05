@testitem "FEM / measurement anchors stop inside their own metal" tags=[:extension] begin
    using Gmsh
    copper=Material(kind=:conductor,rho=2e-8,eps_r=1.,mu_r=1.)
    dielectric=Material(kind=:insulator,rho=1e8,eps_r=3.,mu_r=1.)
    design=build(CableDesign,"measurement-coax",Stack(
        terminal(:core,Region(:metal,Disk(.005),copper)),
        Region(:insulation,Shell(.005),dielectric),
        terminal(:sheath,Region(:metal,Shell(.001),copper))))
    system=build(LineCableSystem,[design],[Pose2(0.,1.)];
        connections=[Dict(:core=>1,:sheath=>2)],system_id="measurement-coax",line_length=1.)
    problem=LineParametersProblem(system;temperature=20.,earth_props=homogeneous(rho=1000.),frequencies=[1e3])
    mktempdir() do dir
        export_data(:onelab,problem,Formulation(:LineCableModelsFEM);
            file_name=joinpath(dir,"model.pro"))
        source=read(joinpath(dir,"geometry","physical.geo"),String)
        values(name)=parse.(Float64,split(only(match(Regex(name*"\\(\\) = \\{([^}]+)\\};"),source).captures),','))
        x,y,d=values("FEMReceiverX"),values("FEMReceiverY"),values("FEMReceiverMetalDimension")
        @test length(x)==length(y)==length(d)==2
        shifted = x .+ min.(1e-5,0.01 .* d)
        @test hypot(shifted[1],y[1]-1.)<.005
        @test .010<hypot(shifted[2],y[2]-1.)<.011
        @test y[2]<.990 # Stops before entering the dielectric or core.
        @test d≈[.010,.001]
    end
end
