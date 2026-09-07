@testitem "Engine / literature assimilation / common pipe matrix assembly" begin
    using LineCableModels, LinearAlgebra
    E=LineCableModels.Engine; DM=LineCableModels.DataModel
    copper=Material(:conductor,1.7241e-8,1.0,1.0,20.0,0.0)
    steel=Material(:conductor,1e-6,1.0,10.0,20.0,0.0)
    oil=Material(:insulator,Inf,2.3,1.0,20.0,0.0)
    function core(name,r)
        Stack(Group(name,Region(name,Disk(r),copper)),
            Region(Symbol(name,:_coat),Shell(0.002),oil))
    end
    function pipe_design(;shell=false)
        contents=DM.Assembly(Pose2(0.0,0.0,0.0),(
            DM.AssemblyMember(core(:a,0.01),Pose2(0.04,0.0,0.0)),
            DM.AssemblyMember(core(:b,0.012),Pose2(-0.02,0.03,0.0))),
            nothing,nothing,nothing,nothing)
        wall=Stack(Group(:pipe,Region(:pipe_wall,
            shell ? Shell(0.006) : Annulus(0.127,0.133),steel)),
            Region(:jacket,Shell(0.002),oil))
        return build(CableDesign,"common-pipe",Enclosure(:cavity,contents;
            primitive=Disk(0.127),fill=oil,wall))
    end
    design=pipe_design()
    blueprint=E.flatten(LineCableModelsCoaxial(),design)
    localdata=E.LocalCableData(blueprint)
    @test design.terminal_order == [:a,:b,:pipe]
    @test blueprint.assembly_ranges == [1:1,2:2,3:3]
    @test length(blueprint.pipes)==1
    @test only(blueprint.pipes).conductor==3
    @test only(blueprint.pipes).children==[1,2]
    function with_pipes(pipes)
        E.CableBlueprint{Float64}(blueprint.cable_id,blueprint.conductors,
            blueprint.dielectrics,blueprint.dielectric_ranges,blueprint.assembly_ranges,pipes)
    end
    @test_throws ArgumentError with_pipes([E.PipeAssembly(4,[1,2],oil)])
    @test_throws ArgumentError with_pipes([E.PipeAssembly(3,[1,4],oil)])
    @test_throws ArgumentError with_pipes([E.PipeAssembly(3,[1,1],oil)])
    @test_throws ArgumentError with_pipes([E.PipeAssembly(3,[1,3],oil)])
    @test_throws ArgumentError with_pipes([E.PipeAssembly(3,[1],oil)])
    @test_throws ArgumentError with_pipes([E.PipeAssembly(1,[2],oil)])
    @test_throws ArgumentError with_pipes([E.PipeAssembly(3,[1,2],steel)])
    @test_throws ArgumentError with_pipes(fill(E.PipeAssembly(3,[1,2],oil),2))
    @test with_pipes(copy(blueprint.pipes)).pipes==blueprint.pipes
    @test E._external_groups(localdata)==(indices=[[1,2,3]],representatives=[3])
    contextual=E.flatten(LineCableModelsCoaxial(),pipe_design(shell=true))
    @test [c.r_in for c in contextual.conductors] ≈ localdata.r_in
    @test [c.r_ex for c in contextual.conductors] ≈ localdata.r_ext
    @test_throws ArgumentError compute(CableConstantsProblem(design),CableConstantsFormulation())

    methods=Formulation(insulation_admittance=:Ametani1980).methods
    layers=zeros(ComplexF64,length(blueprint.dielectrics))
    terms=zeros(ComplexF64,3); tails=similar(terms)
    K=[1.0 0.0;0.0 1.0;-1.0 -1.0]
    ε=8.8541878128e-12*oil.eps_r; μ=4π*1e-7
    r=[0.01,0.012]; b=r .+ 0.002
    z=[complex(0.04,0.0),complex(-0.02,0.03)]
    geometric=[i==j ? log(0.127/b[i])+log1p(-abs2(z[i])/0.127^2) :
        log(abs(0.127^2-z[i]*conj(z[j]))/(0.127*abs(z[i]-z[j]))) for i in 1:2,j in 1:2]
    for (selector,frequency) in ((:DaSilva2006,1.0),(:DaSilva2006,50.0),
            (:DaSilva2006,10000.0),(:Hoidalen2013,1e-4),(:Hoidalen2013,1e-3),
            (:Yang2001,50.0))
        selected_methods=Formulation(pipe_impedance=selector,insulation_admittance=:Ametani1980).methods
        s=complex(0.0,2π*frequency)
        wall=E.PipeImpedance.Formula(selector)(0.127,0.133,steel.rho,steel.mu_r,s)
        metal=E.InternalImpedance.Formula(:Schelkunoff1934)
        core_z=[metal(0.0,r[i],copper.rho,1.0,s)(Val(:outer))+
            s*μ/(2π)*log(b[i]/r[i]) for i in 1:2]
        H=[wall(i==j ? Val(:self) : Val(:mutual),
            E.PipeImpedance.Pair(i,j,((real(z[i]),imag(z[i])),(real(z[j]),imag(z[j]))),
                (b[i],b[j]))) for i in 1:2,j in 1:2]
        a=wall(Val(:outer))+
            s*μ/(2π)*log(0.135/0.133)
        transfer=wall(Val(:mutual))
        reference=fill(a-transfer,3,3)
        reference[1:2,1:2]=Diagonal(core_z)+H .+ a .- 2transfer
        reference[3,3]=a
        Z=zeros(ComplexF64,3,3)
        E.cable_impedance!(Z,localdata,localdata.rho0_cond,selected_methods,s)
        @test Z ≈ reference rtol=1e-11
        @test transpose(K)*Z*K ≈ Diagonal(core_z)+H rtol=1e-11
        @test Z ≈ transpose(Z)
        @test minimum(eigvals(Symmetric(real.(Z))))>=-1e-12

        exterior=log(0.135/0.133)/(2π*ε)
        expected=fill(complex(exterior),3,3)
        inside=Diagonal(log.(b./r)./(2π*ε))+geometric/(2π*ε)
        expected[1:2,1:2]=inside .+ exterior
        P=zeros(ComplexF64,3,3)
        E.cable_potential!(P,localdata,selected_methods,frequency,20.0,s,layers,terms,tails)
        @test P ≈ expected rtol=1e-12
        @test transpose(K)*P*K ≈ inside rtol=1e-12
        @test minimum(eigvals(Symmetric(real.(P))))>0
        @test_throws ArgumentError E.cable_admittance!(similar(P),localdata,
            selected_methods,frequency,20.0,s,layers)
    end

    for (selector,frequencies) in ((:DaSilva2006,[50.0,10000.0]),
            (:Hoidalen2013,[1e-4,1e-3]),(:Yang2001,[50.0,10000.0])), count in (1,2)
        system=build(LineCableSystem,fill(design,count),
            [Pose2(2.0*(i-1),-1.0,0.4*(i-1)) for i in 1:count];
            connections=[Dict("a"=>3i-2,"b"=>3i-1,"pipe"=>3i) for i in 1:count])
        problem=LineParametersProblem(system;earth_props=EarthModel(100.0,10.0,1.0),
            frequencies,temperature=20.0)
        formulation=Formulation(pipe_impedance=formula(selector;rtol=1e-11),
            insulation_admittance=:Ametani1980,earth_impedance=:Pollaczek1926,
            earth_admittance=:IdealGround,
            options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        result=compute(problem,formulation;options=(trace=true,))
        trace=details(result).trace
        @test size(result.Z)==(3count,3count,2)
        @test size(trace.Zg)==(count,count,2)
        @test trace.cable_map==repeat(1:count;inner=3)
        @test all(isfinite,result.Z)
        @test all(isfinite,result.Y)
        for k in 1:2
            s=complex(0.0,2π*problem.frequencies[k])
            @test result.Y[:,:,k] ≈ s*inv(trace.P[:,:,k]) rtol=1e-11
            @test result.Z[:,:,k] ≈ transpose(result.Z[:,:,k])
            for left in 1:count,right in 1:count
                rows=3left-2:3left; columns=3right-2:3right
                @test trace.Z[rows,columns,k]-trace.Zin[rows,columns,k] ≈
                    fill(trace.Zg[left,right,k],3,3) rtol=1e-11
                @test trace.P[rows,columns,k]-trace.Pin[rows,columns,k] ≈
                    fill(trace.Pg[left,right,k],3,3) rtol=1e-11
            end
            for unit in 1:count
                rows=3unit-2:3unit
                @test transpose(K)*trace.Z[rows,rows,k]*K ≈
                    transpose(K)*trace.Zin[rows,rows,k]*K rtol=1e-11
            end
        end
    end
end
