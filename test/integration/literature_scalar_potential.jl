@testitem "Engine / literature assimilation / Kikuchi and Mariscotti scalar potentials" begin
    using LineCableModels, SpecialFunctions, QuadGK, LinearAlgebra
    E=LineCableModels.Engine; EA=E.EarthAdmittance
    mu0=4π*1e-7; eps0=8.8541878128e-12
    # Source integrals evaluated independently on a logarithmic spectral grid.
    function spectral_reference(s,rho,epsilon,mu,kx,h,x,y;source=1,hankel=false)
        kappa=inv.(rho).+s.*epsilon
        bulk=s.*mu.*kappa; chi2=bulk.+kx^2
        root(z)=iszero(imag(z)) && real(z)<0 ?
            complex(0.,sign(imag(s))*sqrt(-real(z))) : sqrt(z)
        chi=root(chi2[source])
        d=hypot(x,y-h); D=hypot(x,y+h)
        points=[-100.,0.,80.]
        for a in chi2
            q=abs(root(a)); q>0 && -100<log(q)<80 && push!(points,log(q))
        end
        sort!(unique!(points))
        integral=quadgk(points...;rtol=1e-10) do t
            lambda=exp(t); q1=root(lambda^2+chi2[1]); q2=root(lambda^2+chi2[2])
            qsource=source==1 ? q1 : q2
            decay=y>=0 ? exp(-(h+y)*qsource) : exp(y*q2-h*q1)
            bulk[source]*decay*cos(x*lambda)*lambda/(bulk[1]*q2+bulk[2]*q1)
        end
        if y<0
            return s/(π*kappa[1])*integral[1]
        end
        direct=if abs(chi2[source])<1e-25
            log(D/d)
        elseif hankel
            im*π/2*(hankelh1(0,im*chi*d)-hankelh1(0,im*chi*D))
        else
            besselk(0,chi*d)-besselk(0,chi*D)
        end
        return s/(2π*kappa[source])*(direct+2integral[1])
    end

    for f in (50.,1e4,1e6),rhoearth in (10.,1000.),h in (1.,10.)
        s=2π*f*im; rho=[Inf,rhoearth]; epsilon=[eps0,10eps0]; mu=[mu0,mu0]
        default=EA.Formula(:Kikuchi1957)(rho,epsilon,mu,s,nothing)
        k0=sqrt(complex(-s*s*mu0*eps0))
        pair=E.EarthPair(1,1,(h,h),.01,(1,1))
        for kx in (k0,complex(.03,.005))
            functor=EA.Formula(:Kikuchi1957)(rho,epsilon,mu,s,kx)
            value=functor(Val(:self),pair)
            @test value ≈ spectral_reference(s,rho,epsilon,mu,kx,h,0.,h-.01;hankel=true) rtol=2e-7
            @test value ≈ EA.earth_potential_coefficient(Val(:Kikuchi1957),
                Val(:observation),functor,h,0.,h-.01) rtol=1e-12
            negative=EA.Formula(:Kikuchi1957)(rho,epsilon,mu,-s,-conj(kx))
            @test negative(Val(:self),pair) ≈ conj(value) rtol=2e-7
            for y in (-.3,0.,2.)
                observation=EA.earth_potential_coefficient(Val(:Kikuchi1957),
                    Val(:observation),functor,h,.4,y)
                @test observation ≈ spectral_reference(s,rho,epsilon,mu,kx,h,.4,y;hankel=true) rtol=2e-7
            end
        end
        # Same air-reference kernel, with the source's lower-surface geometry.
        observation_pair=E.EarthPair(1,2,(h,h-.01),0.,(1,1))
        wise=EA.Formula(:Wise1948)(rho,epsilon,mu,s,nothing)(Val(:mutual),observation_pair)
        @test default(Val(:self),pair) ≈ wise rtol=2e-7
        @test real(s/default(Val(:self),pair))>=0
        # The air-side normal field scales as 1/epsilon_air: a fixed relative
        # potential tolerance at a fixed distance is not a continuity test.
        # Bound the first-order spatial change by its electrostatic field scale
        # and approach the interface from both sides at two distances.
        for y in (-1e-8,1e-8,-1e-10,1e-10)
            boundary=EA.earth_potential_coefficient(Val(:Kikuchi1957),
                Val(:observation),default,h,.4,0.)
            nearby=EA.earth_potential_coefficient(Val(:Kikuchi1957),
                Val(:observation),default,h,.4,y)
            field_bound=4abs(y)/(π*eps0*h)
            @test abs(nearby-boundary)<=field_bound+1e-7abs(boundary)
        end
    end

    for f in (50.,1e4,1e6),source in (1,2),mur in (1.,2.),self in (false,true)
        s=2π*f*im; rho=[Inf,100.]; epsilon=[eps0,10eps0]; mu=[mu0,mur*mu0]
        h=source==1 ? 10. : 1.; y=self ? h : h+.2; x=self ? .01 : .7
        heights=source==1 ? (h,y) : (-h,-y)
        pair=E.EarthPair(1,self ? 1 : 2,heights,x,(source,source))
        k0=sqrt(complex(-s*s*mu0*eps0))
        for kx in (k0,complex(.03,.005))
            functor=EA.Formula(:Mariscotti2019)(rho,epsilon,mu,s,kx)
            value=functor(self ? Val(:self) : Val(:mutual),pair)
            @test value ≈ spectral_reference(s,rho,epsilon,mu,kx,h,x,y;source) rtol=2e-7
            if !self
                reverse_pair=E.EarthPair(2,1,reverse(heights),x,(source,source))
                @test functor(Val(:mutual),reverse_pair) ≈ value rtol=2e-7
            end
        end
    end
    # Identical media away from the light line reduce to the unbounded K0 field.
    for source in (1,2)
        s=100π*im; rho=[100.,100.]; epsilon=[10eps0,10eps0]; mu=[mu0,mu0]
        kx=.03+.005im; h=2.; y=3.; x=.7
        pair=E.EarthPair(1,2,source==1 ? (h,y) : (-h,-y),x,(source,source))
        value=EA.Formula(:Mariscotti2019)(rho,epsilon,mu,s,kx)(Val(:mutual),pair)
        kappa=.01+s*10eps0; chi=sqrt(s*mu0*kappa+kx^2)
        @test value ≈ s/(2π*kappa)*besselk(0,chi*hypot(x,y-h)) rtol=2e-7
    end

    for T in (Float32,Float64,BigFloat),id in (:Kikuchi1957,:Mariscotti2019)
        mu=T(4)*T(π)/T(10)^7; epsilon=T(88541878128)/T(10)^22
        s=complex(zero(T),T(100)*T(π))
        args=(T[Inf,100],T[epsilon,10epsilon],T[mu,mu],s,nothing)
        pair=E.EarthPair(1,1,(T(10),T(10)),T(.01),(1,1))
        value=EA.Formula(id)(args...)(Val(:self),pair)
        @test value isa Complex{T}
        @test isfinite(value)
        expected=EA.Formula(id)([Inf,100.],[eps0,10eps0],[mu0,mu0],100π*im,nothing)(
            Val(:self),E.EarthPair(1,1,(10.,10.),.01,(1,1)))
        @test value ≈ expected rtol=(T===Float32 ? 2e-5 : 2e-7)
    end
    args=([Inf,100.],[eps0,10eps0],[mu0,mu0],100π*im,nothing)
    mixed=E.EarthPair(1,2,(1.,-1.),.1,(1,2))
    @test_throws ArgumentError EA.Formula(:Mariscotti2019)(args...)(Val(:mutual),mixed)
    @test_throws ArgumentError EA.Formula(:Kikuchi1957)(args...)(Val(:mutual),mixed)
    @test_throws DomainError EA.Formula(:Kikuchi1957)([100.,100.],args[2:end]...)
    @test_throws DomainError EA.Formula(:Kikuchi1957)(args[1:2]...,[mu0,2mu0],args[4:end]...)

    copper=Material(:conductor,1.7e-8,1.,1.,20.,0.)
    design=build(CableDesign,"scalar-potential-wire",Group(:core,Region(:core,Disk(.01),copper)))
    for id in (:Kikuchi1957,:Mariscotti2019),source in (1,2)
        id===:Kikuchi1957 && source==2 && continue
        count=id===:Kikuchi1957 ? 1 : 2
        positions=[(i-1.,source==1 ? 10. : -1.) for i in 1:count]
        system=build(LineCableSystem,fill(design,count),positions;
            connections=[Dict("core"=>i) for i in 1:count])
        problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),
            frequencies=[50.,1e5],temperature=20.)
        result=compute(problem,Formulation(earth_impedance=:Pollaczek1926,earth_admittance=id,
            options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false));options=(trace=true,))
        trace=details(result).trace
        for k in 1:2
            s=2π*problem.frequencies[k]*im
            @test result.Y[:,:,k] ≈ s*inv(trace.P[:,:,k]) rtol=1e-10
            @test result.Y[:,:,k] ≈ transpose(result.Y[:,:,k])
            @test minimum(eigvals(Symmetric(real.(result.Y[:,:,k]))))>=-1e-12
            @test trace.P[:,:,k] ≈ trace.Pg[:,:,k]
        end
    end
end
