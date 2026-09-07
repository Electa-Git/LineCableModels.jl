@testitem "Engine / literature assimilation / bonded conductor and semiconductor surfaces" begin
    using LineCableModels,LinearAlgebra
    E=LineCableModels.Engine; II=E.InternalImpedance
    # Independent unscaled Bessel equations and complete circuit elimination.
    function source_surfaces(a,b,rho,mur,s)
        T=typeof(a); piT=one(T)*π
        m=sqrt(s*(4piT*T(10)^(-7))*mur/rho)
        I(n,z)=E.special_besselix(n,z)*exp(abs(real(z)))
        K(n,z)=E.special_besselkx(n,z)*exp(-z)
        if iszero(a)
            return (inner=zero(s),outer=rho*m*I(0,m*b)/(2piT*b*I(1,m*b)),mutual=zero(s))
        end
        D=I(1,m*b)*K(1,m*a)-I(1,m*a)*K(1,m*b)
        (inner=rho*m/(2piT*a*D)*(I(0,m*a)*K(1,m*b)+K(0,m*a)*I(1,m*b)),
         outer=rho*m/(2piT*b*D)*(I(0,m*b)*K(1,m*a)+K(0,m*b)*I(1,m*a)),
         mutual=rho/(2piT*a*b*D))
    end
    function circuit_reference(layers,owners,gaps,s)
        setprecision(BigFloat,512) do
            surfaces=[source_surfaces(BigFloat(l.r_in),BigFloat(l.r_ex),
                BigFloat(l.rho),BigFloat(l.mu_r),Complex{BigFloat}(s)) for l in layers]
            n=length(layers); loop=zeros(Complex{BigFloat},n,n)
            for k in 1:n
                loop[k,k]=surfaces[k].outer+BigFloat(gaps[k])*Complex{BigFloat}(s)
                if k<n
                    loop[k,k]+=surfaces[k+1].inner
                    loop[k,k+1]=loop[k+1,k]=-surfaces[k+1].mutual
                end
            end
            transform=BigFloat[i>=j for i in 1:n,j in 1:n]
            primitive=transpose(transform)*loop*transform
            bonds=BigFloat[owners[i]==j for i in 1:n,j in 1:maximum(owners)]
            inv(transpose(bonds)*(primitive\bonds))
        end
    end

    for T in (Float32,Float64,BigFloat),solid in (false,true),
            f in (50,100000),rho2 in (1e-7,1e3)
        a=solid ? zero(T) : T(.005); b=T(.01); c=T(.011)
        rho1=T(1.7e-8); mu1=T(2); mu2=one(T)
        s=complex(zero(T),T(2)*T(π)*T(f))
        layers=[(r_in=a,r_ex=b,rho=rho1,mu_r=mu1),
                (r_in=b,r_ex=c,rho=T(rho2),mu_r=mu2)]
        value=II.Formula(:Ametani2004)(layers,s)
        expected=circuit_reference(layers,[1,1],[0.,0.],s)[1,1]
        tolerance=T===Float32 ? 3e-6 : 3e-13
        @test value(Val(:outer)) ≈ expected rtol=tolerance
        @test value(Val(:outer)) isa Complex{T}
        @test value(Val(:outer))==II.Formula(:Ametani2004)(
            a,b,b,c,rho1,T(rho2),mu1,mu2,s)(Val(:outer))
        @test real(value(Val(:outer)))>0
        negative=II.Formula(:Ametani2004)(layers,-s)
        for mode in (:inner,:outer,:mutual)
            @test negative(Val(mode)) ≈ conj(value(Val(mode))) rtol=tolerance
        end
        if !solid
            # Direct two-material Maxwell continuity, expressed as a boundary
            # solve in surface variables, tests both current excitations.
            setprecision(BigFloat,512) do
                l1,l2=layers
                z1=source_surfaces(BigFloat(a),BigFloat(b),BigFloat(rho1),BigFloat(mu1),Complex{BigFloat}(s))
                z2=source_surfaces(BigFloat(b),BigFloat(c),BigFloat(rho2),BigFloat(mu2),Complex{BigFloat}(s))
                for (inside,outside) in ((1,0),(0,1))
                    interface=(z1.mutual*inside+z2.mutual*outside)/(z1.outer+z2.inner)
                    ein=z1.inner*inside-z1.mutual*interface
                    eout=z2.outer*outside-z2.mutual*interface
                    @test value(Val(:inner))*inside-value(Val(:mutual))*outside ≈ ein rtol=tolerance
                    @test value(Val(:outer))*outside-value(Val(:mutual))*inside ≈ eout rtol=tolerance
                end
            end
        end
        dc=II.Formula(:Ametani2004)(layers,zero(s))
        resistance=inv(sum(l->(T(π)*(l.r_ex^2-l.r_in^2))/l.rho,layers))
        @test dc(Val(:outer)) ≈ resistance rtol=tolerance
        @test iszero(imag(dc(Val(:outer))))
    end

    # Splitting a homogeneous annulus does not change the physical formula.
    for a in (0.,.003),f in (50.,1e6)
        radii=[a,.006,.008,.011]
        layers=[(r_in=radii[k],r_ex=radii[k+1],rho=1.7e-8,mu_r=1.) for k in 1:3]
        s=2π*f*im; combined=II.Formula(:Ametani2004)(layers,s)
        homogeneous=II.Formula(:Schelkunoff1934)(a,last(radii),1.7e-8,1.,s)
        for mode in (:inner,:outer,:mutual)
            @test combined(Val(mode)) ≈ homogeneous(Val(mode)) rtol=3e-13
        end
    end
    @test_throws ArgumentError II.Formula(:Ametani2004)(
        .001,.002,.0021,.003,1e-8,1.,1.,1.,100π*im)
    @test_throws DomainError II.Formula(:Ametani2004)(
        .001,.002,.002,.003,1e-8,-1.,1.,1.,100π*im)

    copper=Material(:conductor,1.7e-8,1.,1.,20.,.004)
    metal2=Material(:conductor,2.8e-8,1.,2.,20.,.003)
    semi=Material(:semicon,100.,100.,1.,20.,.001)
    dielectric=Material(:insulator,1e12,2.3,1.,20.,0.)
    # N-screen source geometry: screen outside the core, inside each later metal.
    for n in (1,2,4),screened in (false,true),temperature in (20.,80.)
        parts=LineCableModels.DataModel.AbstractCablePart[]
        push!(parts,Group(:core,Stack(
            Region(:inner_core,Disk(.006),copper),
            Region(:outer_core,Shell(.004),metal2))))
        screened && push!(parts,Region(:screen1,Shell(.001),semi))
        push!(parts,Region(:insulation1,Shell(.004),dielectric))
        for k in 2:n
            screened && push!(parts,Region(Symbol(:screen,k),Shell(.001),semi))
            push!(parts,Group(Symbol(:conductor,k),Region(Symbol(:metal,k),Shell(.001),copper)))
            push!(parts,Region(Symbol(:insulation,k),Shell(.004),dielectric))
        end
        design=build(CableDesign,"bonded-$n-$screened",Stack(parts))
        blueprint=E.flatten(LineCableModelsCoaxial(),design)
        @test length(blueprint.conductors[1].layers)==2
        @test blueprint.conductors[1].layers[1].material.rho==copper.rho
        input=E.LocalCableData(blueprint)
        formula=Formulation(internal_impedance=:Ametani2004,
            options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
        s=100π*im; computed=zeros(ComplexF64,n,n)
        E.cable_impedance!(computed,input,input.rho0_cond,formula.methods,s;temperature)
        layers=NamedTuple[];owners=Int[];gaps=Float64[]
        physical=design.geometry.regions
        nextmetal=1
        for (r,region) in pairs(physical)
            material=region.source.material
            a=E.DataModel.r_in(region.primitive);b=E.DataModel.r_ex(region.primitive)
            if material.kind===:insulator
                gaps[end]+=4π*1e-7*material.mu_r/(2π)*log(b/a)
            else
                owner=if region.terminal!==nothing
                    findfirst(==(region.terminal),input.terminals)
                elseif r>1 && physical[r-1].terminal===:core
                    1
                else
                    target=findnext(p->p.terminal!==nothing,physical,r+1)
                    findfirst(==(physical[target].terminal),input.terminals)
                end
                rho=material.rho*(1+material.alpha*(temperature-material.T0))
                push!(layers,(r_in=a,r_ex=b,rho=rho,mu_r=material.mu_r))
                push!(owners,owner);push!(gaps,0.)
            end
        end
        reference=circuit_reference(layers,owners,gaps,s)
        @test computed ≈ reference rtol=4e-12
        @test computed ≈ transpose(computed)
        @test minimum(eigvals(Symmetric(real.(computed))))>0
        # Shunt layers are unchanged by the choice of series-surface formula.
        result=compute(CableConstantsProblem(design;temperature),
            CableConstantsFormulation(internal_impedance=:Ametani2004))
        plain=compute(CableConstantsProblem(design;temperature),CableConstantsFormulation())
        @test result.C ≈ plain.C rtol=1e-12
        @test result.G ≈ plain.G rtol=1e-12
    end
end
