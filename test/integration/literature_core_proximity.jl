@testitem "Engine / literature assimilation / pipe core proximity and infinite wall" begin
    using LineCableModels, SpecialFunctions, LinearAlgebra
    E=LineCableModels.Engine; PI=E.PipeImpedance; DM=LineCableModels.DataModel
    mu0=4π*1e-7
    function literal(id,r,d,rho,mur,s)
        z=r*sqrt(s*mu0*mur/rho); g=(r/d)^2
        if id===:Kane1995
            # Original Kane (10),(17) and Da Silva (4), without derivative recurrence.
            value=sum(1:70) do n
                I=besselix(n,z); derivative=(besselix(n-1,z)+besselix(n+1,z))/2
                g^n*I/(n*I+z/mur*derivative)
            end
            return s*mu0/π*value
        end
        k1=g*besselix(2,z)/besselix(0,z)
        value=sum(g^n/z*besselix(n,z)/besselix(n-1,z)*
            abs2(1+n*k1/(1-k1)) for n in 1:70)
        return s*mu0/(2π)*(value+log1p(-g)/2)
    end
    for f in (.01,1.,50.,1e4,1e7), d in (.021,.04,.2),mur in (1.,5.),
            id in (:Kane1995,:Hoidalen2013)
        id===:Hoidalen2013 && mur!=1 && continue
        s=2π*f*im
        value=PI.core_proximity(Val(id),.01,d,1.7e-8,mur,s;rtol=1e-12)
        @test value ≈ literal(id,.01,d,1.7e-8,mur,s) rtol=1e-7
        @test real(value)>=0
        @test PI.core_proximity(Val(id),.01,d,1.7e-8,mur,-s;rtol=1e-12) ≈ conj(value)
    end
    for T in (Float32,Float64,BigFloat),id in (:Kane1995,:Hoidalen2013)
        r=T(.01); d=T(.04); rho=T(1.7e-8); mur=one(T)
        for f in (T(1e-9),T(50))
            s=complex(zero(T),2*T(π)*f)
            v=PI.core_proximity(Val(id),r,d,rho,mur,s)
            @test v isa Complex{T}
            @test isfinite(v)
            @test real(v)>0
            if f<T(.01)
                g=(r/d)^2
                # The removed static term and leading dissipative O(omega^2) term.
                pref=s*(T(4)*T(π)/T(10)^7)/(2*T(π))
                z2=s*(T(4)*T(π)/T(10)^7)*r^2/rho
                leading=-pref*z2*sum(g^n/(8n^2*(n+1)) for n in 1:50)
                if id===:Kane1995
                    @test imag(v/pref) ≈ imag(2leading/pref) rtol=1e-4
                    @test real(v/pref) ≈ -log1p(-g) rtol=1e-5
                else
                    @test real(v) ≈ real(leading) rtol=1e-5
                    @test abs(imag(v))<=abs(real(v))*T(1e-6)
                end
            end
        end
    end
    @test_throws DomainError PI.core_proximity(Val(:Hoidalen2013),.01,.04,1e-8,2.,im*100.)
    @test_throws DomainError PI.core_proximity(Val(:Kane1995),.04,.01,1e-8,1.,im*100.)
    @test_throws ErrorException PI.core_proximity(Val(:Kane1995),.01,.021,1e-8,1.,im*100.;max_terms=1)

    a=.127; b=.133; rho=1e-6
    positions=[(.03,0.),(-.015,.03sqrt(3)/2),(-.015,-.03sqrt(3)/2)]
    cores=PI.Cores([1,2,3],positions,fill(.01,3),fill(1.7e-8,3),ones(3))
    for f in (1.,50.,1e5),id in (:Kane1995,:Hoidalen2013)
        s=2π*f*im
        parent=PI.Formula(:DaSilva2006)(a,b,rho,10.,s)
        candidate=id===:Kane1995 ? PI.Formula(:Kane1995) :
            PI.Formula(:DaSilva2006;proximity=:Hoidalen2013)
        wall=PI.with_cores(candidate(a,b,rho,10.,s),cores)
        d=hypot((positions[1].-positions[2])...)
        delta=PI.core_proximity(Val(id),.01,d,1.7e-8,1.,s)
        difference=zeros(ComplexF64,3,3)
        for i in 1:3,j in 1:3
            pair=PI.Pair(i,j,(positions[i],positions[j]),(.012,.012))
            kind=i==j ? Val(:self) : Val(:mutual)
            difference[i,j]=wall(kind,pair)-parent(kind,pair)
        end
        @test difference ≈ delta*(ones(3,3)+I) rtol=1e-8
        @test difference*[1.,-1.,0.] ≈ delta*[1.,-1.,0.] atol=1e-12
        @test difference*ones(3) ≈ 4delta*ones(3) rtol=1e-8
        @test wall(Val(:outer))==parent(Val(:outer))
        @test wall(Val(:mutual))==parent(Val(:mutual))
        unequal=PI.Cores([1,2,3],positions,[.01,.009,.01],fill(1.7e-8,3),ones(3))
        @test_throws ArgumentError PI.with_cores(candidate(a,b,rho,10.,s),unequal)
    end

    for f in (1e-7,1e-5,1e-3),mur in (1.,10.),i in 1:3,j in i:3
        s=2π*f*im; x=a*sqrt(s*mu0*mur/rho); pref=s*mu0/(2π)
        z=complex.(first.(positions),last.(positions))
        w=z[i]*conj(z[j])/a^2
        pair=PI.Pair(i,j,(positions[i],positions[j]),(.012,.012))
        kind=i==j ? Val(:self) : Val(:mutual)
        geometric=i==j ? log((a^2-abs2(z[i]))/(a*.012)) :
            log(abs(a^2-z[i]*conj(z[j]))/(a*abs(z[i]-z[j])))
        # |w|<0.06 here: 30 literal terms leave a negligible geometric tail
        # without overflowing the independent unscaled-order K evaluation.
        harmonic=sum(1:30) do n
            K=besselkx(n,x); derivative=-(besselkx(n-1,x)+besselkx(n+1,x))/2
            2mur*real(w^n)*K/(n*mur*K-x*derivative)
        end
        full=PI.Formula(:DaSilva2006;wall=:infinite,proximity=:none)(a,Inf,rho,mur,s)
        expected=pref*(mur*besselkx(0,x)/(x*besselkx(1,x))+geometric+harmonic)
        @test full(kind,pair) ≈ expected rtol=1e-10
        @test_throws ArgumentError full(Val(:outer))
        @test_throws ArgumentError full(Val(:mutual))
        low=PI.Formula(:Hoidalen2013;wall=:infinite)(a,Inf,rho,mur,s)
        L=log(2/x)-Base.MathConstants.eulergamma
        scale=2mur/(1+mur); Delta=scale*real(w)*L*x^2/(1+mur+L*x^2)
        # Parent (3) gives minus Delta, although printed (15) has plus Delta.
        href=-scale*log(abs(1-w))-Delta
        @test low(kind,pair) ≈ pref*(mur*L+geometric+href) rtol=1e-10
        # Source (14) retains only the leading K0 logarithm. At the largest
        # tested magnetic skin argument its omitted terms cause about 0.3%.
        @test low(kind,pair) ≈ full(kind,pair) rtol=5e-3
        @test abs(href-harmonic) <= abs(-scale*log(abs(1-w))+Delta-harmonic)+1e-14
        @test_throws ArgumentError PI.Formula(:DaSilva2006;wall=:infinite)(a,b,rho,mur,s)
    end
    s=100π*im
    raw=PI.Formula(:DaSilva2006;wall=:infinite)(a,Inf,rho,10.,s)
    coupled=PI.with_cores(raw,cores)
    pair=PI.Pair(1,1,(positions[1],positions[1]),(.012,.012))
    bare=PI.Formula(:DaSilva2006;wall=:infinite,proximity=:none)(a,Inf,rho,10.,s)
    @test coupled(Val(:self),pair)-bare(Val(:self),pair) ≈
        sum(PI.core_proximity(Val(:Kane1995),.01,hypot((positions[1].-positions[k])...),
            1.7e-8,1.,s) for k in 2:3) rtol=1e-8
    @test_throws ArgumentError raw(Val(:self),pair)

    copper=Material(:conductor,1.7e-8,1.,1.,20.,0.)
    steel=Material(:conductor,rho,1.,10.,20.,0.)
    oil=Material(:insulator,Inf,2.3,1.,20.,0.)
    function core(name)
        Stack(Group(name,Region(name,Disk(.01),copper)),Region(Symbol(name,:_coat),Shell(.002),oil))
    end
    contents=DM.Assembly(Pose2(0.,0.,0.),Tuple(
        DM.AssemblyMember(core(Symbol(:c,i)),Pose2(positions[i]...,0.)) for i in 1:3),
        nothing,nothing,nothing,nothing)
    wall=Stack(Group(:pipe,Region(:wall,Annulus(a,b),steel)),Region(:jacket,Shell(.002),oil))
    design=build(CableDesign,"proximity-pipe",Enclosure(:cavity,contents;primitive=Disk(a),fill=oil,wall))
    system=build(LineCableSystem,[design],[(0.,-1.)];
        connections=[Dict("c1"=>1,"c2"=>2,"c3"=>3,"pipe"=>4)])
    problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),
        frequencies=[1.,50.,1e5],temperature=20.)
    function calculate(selection)
        compute(problem,Formulation(pipe_impedance=selection,insulation_admittance=:Ametani1980,
            earth_impedance=:Pollaczek1926,earth_admittance=:IdealGround,
            options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false));options=(trace=true,))
    end
    baseline=calculate(:DaSilva2006)
    for id in (:Kane1995,:Hoidalen2013)
        selection=id===:Kane1995 ? formula(:Kane1995) :
            formula(:DaSilva2006;proximity=:Hoidalen2013)
        result=calculate(selection)
        @test result.Y ≈ baseline.Y
        @test size(result.Z)==(4,4,3)
        for k in 1:3
            local s=2π*problem.frequencies[k]*im
            delta=PI.core_proximity(Val(id),.01,hypot((positions[1].-positions[2])...),1.7e-8,1.,s)
            expected=zeros(ComplexF64,4,4)
            expected[1:3,1:3]=delta*(ones(3,3)+I)
            @test result.Z[:,:,k]-baseline.Z[:,:,k] ≈ expected rtol=1e-8 atol=1e-14
            @test minimum(eigvals(Symmetric(real.(result.Z[:,:,k]))))>=-1e-12
        end
    end
end
