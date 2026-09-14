@testitem "Gmsh FEM / quasi-TEM / independent annular equation and mesh refinement" tags=[:extension,:integration,:fem_numerical] setup=[TestFixtures] begin
    using Gmsh, LinearAlgebra, TOML
    include(joinpath(pkgdir(LineCableModels),"test/support/radial_control.jl"))
    const FEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    BLAS.set_num_threads(1)
    started=time(); jobs=Ref(0)
    directory=mktempdir(;prefix="lcm-fem-provisional-",cleanup=false)
    println("FEM provisional evidence: ",directory);flush(stdout)
    controls=(gmsh_verbosity=0,getdp_verbosity=0,frequency_workers=1,solver_threads=1,
        keep_run_directory=true,trace=true)
    design=TestFixtures.coaxial_design()
    for f in (50.0,10000.0), grounded in (false,true)
        connections=Dict(:core=>1,:sheath=>(grounded ? 0 : 2))
        system=build(LineCableSystem,design,Pose2(0.0,-1.0);connections)
        problem=LineParametersProblem(system;frequencies=[f],earth_props=homogeneous(rho=100.0,eps_r=10.0))
        formulation=Formulation(:LineCableModelsFEM;insulation_admittance=:Ametani2004,
            options=(physics=:quasi_tem,reduce_bundle=false,kron_reduction=grounded,ideal_transposition=false))
        references=[setprecision(BigFloat,bits) do
            s=2big(pi)*im*f; kappa=inv(big"1e8")+s*3big"8.8541878128e-12"
            RadialControl.annular_admittance(big".005",big".01",kappa,4big(pi)*big"1e-7",s)
        end for bits in (128,256,512)]
        reference=last(references)
        uncertainty=reference.bound+abs(reference.value-references[2].value)
        samples=ComplexF64[]; meshes=String[]; triangles=Int[]; resolved=false
        for level in 0:3
            time()-started<45*60 && jobs[]<128 || error("inconclusive FEM: batch resource limit")
            options=controls
            if level>0
                refined=joinpath(directory,"f$(Int(f))-grounded$grounded-level$level.msh")
                session=FEM._start_gmsh(0)
                try
                    gmsh.open(last(meshes))
                    kinds,tags,_=gmsh.model.mesh.get_elements(2)
                    all(==(2),kinds) || error("control requires first-order triangles")
                    count=sum(length,tags)
                    4count<=500000 || error("inconclusive FEM: refinement exceeds 500000 triangles")
                    gmsh.model.mesh.refine()
                    gmsh.write(refined)
                finally
                    FEM._finish_gmsh(session)
                end
                options=merge(controls,(mesh_path=refined,mesh_policy=:reuse))
            end
            timed=@timed compute(problem,formulation;options)
            result=timed.value;jobs[]+=1
            timed.time<=300 || error("inconclusive FEM: solve exceeded five minutes")
            @test size(Y(result))==(grounded ? (1,1,1) : (2,2,1))
            @test length(unique(string.(axes(result.Y,1))))==size(result.Y,1)
            run=result.details.fem.run.run_directory
            mesh=joinpath(run,"mesh","model.msh");push!(meshes,mesh)
            session=FEM._start_gmsh(0)
            try
                gmsh.open(mesh)
                kinds,tags,_=gmsh.model.mesh.get_elements(2)
                @test all(==(2),kinds)
                push!(triangles,sum(length,tags))
                @test last(triangles)<=500000
            finally
                FEM._finish_gmsh(session)
            end
            value=Y(result)[1,1,1];push!(samples,value)
            # The closed grounded sheath separates the interior Dirichlet
            # annulus from the exterior. For this entry the domain error is zero
            # analytically; the independent BVP directly bounds mesh error.
            resolved=all(component->begin
                budget=.01abs(component(reference.value))
                uncertainty<=budget/4 && abs(component(value-reference.value))+uncertainty<=budget/4
            end,(real,imag))
            open(joinpath(directory,"annular-$(Int(f))-$grounded.toml"),"w") do io
                TOML.print(io,Dict("status"=>(resolved ? "resolved control" : "inconclusive refinement"),
                    "frequency"=>f,"grounded"=>grounded,"triangles"=>triangles,"meshes"=>meshes,
                    "real"=>real.(samples),"imag"=>imag.(samples),
                    "reference_real"=>string(real(reference.value)),"reference_imag"=>string(imag(reference.value)),
                    "reference_uncertainty"=>string(uncertainty)))
            end
            resolved && break
        end
        @test resolved
        # A resolved refined solve cannot certify an inaccurate default mesh.
        for component in (real,imag)
            @test abs(component(first(samples)-reference.value))+uncertainty<=.01abs(component(reference.value))
        end
    end
end

@testitem "Gmsh FEM / quasi-TEM / differential low-frequency magnetic limit" tags=[:extension,:integration,:fem_numerical] setup=[TestFixtures] begin
    using Gmsh,LinearAlgebra,QuadGK,TOML
    include(joinpath(pkgdir(LineCableModels),"test/support/radial_control.jl"))
    FEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    BLAS.set_num_threads(1)
    directory=mktempdir(;prefix="lcm-fem-magnetic-provisional-",cleanup=false)
    println("FEM magnetic evidence: ",directory);flush(stdout)
    # For unit differential current, integrate the independently specified
    # enclosed-current law. This limit includes both finite metal regions.
    limit=setprecision(BigFloat,256) do
        a,b,c,rho=big".005",big".01",big".011",big"2e-8"
        mu=4big(pi)*big"1e-7"
        energy,error=quadgk(r->((c^2-r^2)/(c^2-b^2))^2/r,b,c;
            rtol=big"1e-12",maxevals=10^6)
        (R=rho/(big(pi)*a^2)+rho/(big(pi)*(c^2-b^2)),
            L=mu/(2big(pi))*(big".25"+log(b/a)+energy),
            uncertainty=mu/(2big(pi))*error)
    end
    design=TestFixtures.coaxial_design()
    system=build(LineCableSystem,design,Pose2(0.,-1.);
        connections=Dict(:core=>1,:sheath=>2))
    formulation=Formulation(:LineCableModelsFEM;insulation_admittance=:Ametani2004,
        options=(physics=:quasi_tem,reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    controls=(gmsh_verbosity=0,getdp_verbosity=0,frequency_workers=1,solver_threads=1,
        keep_run_directory=true,trace=true)
    rows=Dict{String,Any}[]
    previous=Ref(Inf)
    for frequency in (1.0,.1,.01)
        problem=LineParametersProblem(system;frequencies=[frequency],
            earth_props=homogeneous(rho=100.,eps_r=10.))
        references=[setprecision(BigFloat,bits) do
            s=2big(pi)*im*frequency;mu=4big(pi)*big"1e-7"
            core=RadialControl.surfaces(big"0",big".005",big"2e-8",mu,s)
            sheath=RadialControl.surfaces(big".01",big".011",big"2e-8",mu,s)
            (value=core.outer+sheath.inner+s*mu/(2big(pi))*log(big"2"),
                bound=core.bound+sheath.bound)
        end for bits in (128,256,512)]
        reference=last(references);omega=2pi*frequency
        reference_error=reference.bound+abs(reference.value-references[2].value)
        finite_R=abs(real(reference.value)-limit.R)+reference_error
        finite_L=abs(imag(reference.value)/omega-limit.L)+reference_error/omega+limit.uncertainty
        @test finite_R/limit.R+finite_L/limit.L < previous[]
        previous[]=Float64(finite_R/limit.R+finite_L/limit.L)
        result=compute(problem,formulation;options=controls)
        drive=[1.,-1.];response=dot(drive,Z(result)[:,:,1]*drive)
        @test size(Z(result))==(2,2,1)
        @test isfinite(response)
        push!(rows,Dict("frequency_Hz"=>frequency,"run"=>result.details.fem.run.run_directory,
            "R_Ohm_per_m"=>real(response),"L_H_per_m"=>imag(response)/omega,
            "R_limit"=>string(limit.R),"L_limit"=>string(limit.L),
            "metal_diffusion_R_bound"=>string(finite_R),"metal_diffusion_L_bound"=>string(finite_L),
            "status"=>"inconclusive: exterior conduction/displacement and mesh/domain errors not yet bounded"))
        open(joinpath(directory,"low-frequency.toml"),"w") do io
            TOML.print(io,Dict("cases"=>rows,
                "reference_scope"=>"uniform-current limit and finite metal diffusion; no exterior field bound"))
        end
    end
    include(joinpath(pkgdir(LineCableModels),"test/support/fem_mesh_control.jl"))
    include(joinpath(pkgdir(LineCableModels),"test/support/edge_field_control.jl"))
    low=LineParametersProblem(system;frequencies=[.01],earth_props=homogeneous(rho=100.,eps_r=10.))
    drive=[1.,-1.]
    study=FEMMeshControl.study(FEM,low,formulation,controls;directory,
        observable=(z,p)->[dot(drive,z*drive)])
    for component in (real,imag)
        scale=component===real ? Float64(limit.R) : 2pi*.01Float64(limit.L)
        mesh=EdgeFieldControl.richardson([component(only(v)) for v in study.mesh_values],.01scale)
        @test mesh.resolved
        differences=abs.(diff([component(only(v)) for v in study.domain_values]))
        @test length(differences)==2 && last(differences)<=first(differences)/2
        @test 2last(differences)<=.01scale/4
    end
    # A shielded radial metal solution cannot certify the full native model:
    # its conducting exterior and displacement terms require their own bound.
    error("inconclusive S7 magnetic evidence: exterior conduction/displacement and mesh/domain allocations remain unresolved; see $directory")
end

@testitem "Gmsh FEM / quasi-TEM / current voltage and mixed-excitation extraction" tags=[:extension,:integration,:fem_numerical] begin
    using Gmsh,LinearAlgebra,Printf
    FEM=Base.get_extension(LineCableModels,:LineCableModelsGmshExt)
    metal=Material(kind=:conductor,rho=2e-8)
    designs=[build(CableDesign,"electrode-$i",terminal(:core,Region(:metal,Disk(r),metal)))
        for (i,r) in enumerate((.005,.007))]
    system=build(LineCableSystem,designs,[(0.0,-1.0),(.4,-1.5)];
        connections=[Dict(:core=>1),Dict(:core=>2)])
    problem=LineParametersProblem(system;frequencies=[10000.0],earth_props=homogeneous(rho=100.0,eps_r=10.0))
    formulation=Formulation(:LineCableModelsFEM;options=(physics=:quasi_tem,
        reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    controls=(gmsh_verbosity=0,getdp_verbosity=0,solver_threads=1,frequency_workers=1,
        keep_run_directory=true,trace=true)
    actual=compute(problem,formulation;options=controls)
    run_directory=actual.details.fem.run.run_directory
    println("FEM extraction evidence: ",run_directory);flush(stdout)
    directory=mktempdir(;prefix="lcm-fem-excitation-provisional-",cleanup=false)
    println("FEM native excitation evidence: ",directory);flush(stdout)
    begin
        for name in readdir(joinpath(run_directory,"input/getdp"))
            cp(joinpath(run_directory,"input/getdp",name),joinpath(directory,name))
        end
        operator_path=joinpath(directory,"quasi-tem.pro");source=Ref(read(operator_path,String))
        edits=("Value \$FEM_Q~{t};"=>"Value -\$FEM_Q~{t};",
            "NameOfCoef Q; EntityType Auto; NameOfConstraint FEMTransverseCurrent;"=>
                "NameOfCoef V; EntityType Auto; NameOfConstraint FEMTransverseCurrent;",
            "Re[{V} / UnitTransverseSource]"=>"-Re[{Q} / UnitTransverseSource]",
            "Im[{V} / UnitTransverseSource]"=>"-Im[{Q} / UnitTransverseSource]")
        for (before,after) in edits
            @test count(before,source[])==1
            source[]=replace(source[],before=>after)
        end
        write(operator_path,source[])
        command=FEM._getdp_command(actual.details.fem.inputs.getdp_provenance.path,
            joinpath(directory,"model.pro"),joinpath(run_directory,"mesh/model.msh"),
            (path=run_directory,),formulation,computation_options(LineCableModelsFEM,controls),
            only(actual.details.fem.inputs.mesh_plans),[1,2],directory;reuse_factorization=false)
        open(joinpath(directory,"getdp.log"),"w") do io
            run(pipeline(command;stdout=io,stderr=io))
        end
        currents=zeros(ComplexF64,2,2)
        for basis in 1:2
            path=joinpath(directory,"raw/jobs",@sprintf("getdp-f0001-b%04d-P.tsv",basis))
            @test FEM._valid_job_raw(path,2,1,10000.0,basis)
            for line in eachline(path)
                cells=split(line)
                currents[parse(Int,cells[3]),basis]=complex(parse(Float64,cells[5]),parse(Float64,cells[6]))
            end
        end
        for index in eachindex(currents)
            @test currents[index] ≈ Y(actual)[index] rtol=1e-10 atol=0
        end
        # Solve the native operator again with mixed voltage columns [1 1;0 2].
        # This changes the drive only; no operator coefficients are regenerated.
        voltage=[1.0 1.0;0.0 2.0]
        drive="-UnitTransverseSource * (\$FEMBasisTerminal == t), 0.]"
        mixed="-UnitTransverseSource * ((t == 1) + 2*(t == 2)*(\$FEMBasisTerminal == 2)), 0.]"
        @test count(drive,source[])==1
        write(operator_path,replace(source[],drive=>mixed))
        mixed_directory=joinpath(directory,"mixed");mkpath(mixed_directory)
        mixed_command=FEM._getdp_command(actual.details.fem.inputs.getdp_provenance.path,
            joinpath(directory,"model.pro"),joinpath(run_directory,"mesh/model.msh"),
            (path=run_directory,),formulation,computation_options(LineCableModelsFEM,controls),
            only(actual.details.fem.inputs.mesh_plans),[1,2],mixed_directory;reuse_factorization=false)
        open(joinpath(mixed_directory,"getdp.log"),"w") do io
            run(pipeline(mixed_command;stdout=io,stderr=io))
        end
        mixed_currents=zeros(ComplexF64,2,2)
        for basis in 1:2
            raw=joinpath(mixed_directory,"raw/jobs",@sprintf("getdp-f0001-b%04d-P.tsv",basis))
            @test FEM._valid_job_raw(raw,2,1,10000.0,basis)
            for line in eachline(raw)
                cells=split(line)
                mixed_currents[parse(Int,cells[3]),basis]=complex(parse(Float64,cells[5]),parse(Float64,cells[6]))
            end
        end
        @test mixed_currents/voltage ≈ Y(actual)[:,:,1] rtol=1e-10

    end
end
