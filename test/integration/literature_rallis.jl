@testitem "Engine / literature assimilation / Rallis complex images" begin
    using LineCableModels,QuadGK,LinearAlgebra
    E=LineCableModels.Engine;EI=E.EarthImpedance
    # Recover a known exponential sum using the source pencil, then test it
    # at points that were not used to form either Hankel matrix.
    poles=ComplexF64[-.2-.3im,-1.1+.4im,-4-.5im]
    residues=ComplexF64[2+im,-.2+.4im,.3-.1im]
    original=t->sum(residues.*exp.(poles*t))
    fit=EI._rallis_gpof(original,3.0,0.0,100,14)
    @test length(fit.poles)==3
    for t in (.01,.137,.77,1.291,3.5,5.)
        @test EI._rallis_spectrum(fit,t) ≈ original(t) rtol=2e-11
    end
    # The final rational and K1 sums must equal integrals of the fitted
    # spectral functions, independently of their error against the parent.
    kap=(1+im)/sqrt(2)
    for levels in (1,2),kind in (:overhead,:mixed,:underground)
        depth=kind===:mixed ? .2 : 0.
        fitted=EI._rallis_fit(kind,kap,depth,levels,100,14)
        @test all(p->real(p)<0,fitted.poles)
        for H in (.3,2.),x in (0.,.5,2.)
            integrand=kind===:underground ?
                t->exp(-H*sqrt(t^2+kap^2))*EI._rallis_spectrum(fitted,sqrt(t^2+kap^2))*cos(x*t) :
                t->exp(-H*t)*EI._rallis_spectrum(fitted,t)*cos(x*t)
            integral=2quadgk(integrand,0.,Inf;rtol=2e-10)[1]
            @test EI._rallis_correction(fitted,kap,H,x,Val(kind)) ≈ integral rtol=2e-8
        end
    end
    μ0=4π*1e-7;ε0=8.8541878128e-12
    # Parent comparisons in a stated finite validation domain. DCIM has no
    # universal error bound, particularly for widely separated buried wires.
    for frequency in (50.,5000.),evaluation in (:dcim,:dcim_two_level)
        s=2π*frequency*im
        args=([Inf,100.],[ε0,10ε0],[μ0,μ0],s,nothing)
        parent=EI.Formula(:Pollaczek1926)(args...)
        approximate=EI.Formula(:Pollaczek1926;evaluation)(args...)
        negative=EI.Formula(:Pollaczek1926;evaluation)(
            args[1],args[2],args[3],-s,nothing)
        @test EI.formula_id(approximate.state.formula)===:Pollaczek1926
        for kind in (:overhead,:mixed,:underground),x in (1.,10.,100.)
            heights,layers=kind===:overhead ? ((10.,2.),(1,1)) :
                kind===:mixed ? ((10.,-2.),(1,2)) : ((-10.,-2.),(2,2))
            pair=E.EarthPair(1,2,heights,x,layers)
            reversed_pair=E.EarthPair(2,1,reverse(heights),x,reverse(layers))
            value=approximate(Val(:mutual),pair)
            reference=parent(Val(:mutual),pair)
            tolerance=kind===:underground ? .04 :
                evaluation===:dcim ? .002 : 2e-6
            @test value ≈ reference rtol=tolerance
            @test approximate(Val(:mutual),reversed_pair) ≈ value rtol=2e-10
            @test negative(Val(:mutual),pair) ≈ conj(value) rtol=2e-8
        end
    end
    for T in (Float32,Float64,BigFloat)
        args=(T[Inf,100],T[ε0,10ε0],T[μ0,μ0],
            complex(zero(T),T(100)*T(π)),nothing)
        fitted=EI.Formula(:Pollaczek1926;evaluation=:dcim_two_level)(args...)
        for (height,layer) in ((T(10),1),(T(-10),2))
            pair=E.EarthPair(1,1,(height,height),T(.01),(layer,layer))
            value=fitted(Val(:self),pair)
            @test value isa Complex{T}
            @test isfinite(value) && real(value)>0
        end
    end
    @test_throws ArgumentError EI.Formula(:Pollaczek1926;evaluation=:dcim,dcim_terms=51)
    @test_throws ArgumentError EI.Formula(:Pollaczek1926;evaluation=:dcim,dcim_samples=3)
    for s in (0im,1.0+im)
        @test_throws DomainError EI.Formula(:Pollaczek1926;evaluation=:dcim)(
            [Inf,100.],[ε0,10ε0],[μ0,μ0],ComplexF64(s),nothing)
    end

    copper=Material(:conductor,1.7241e-8,1.,1.,20.,0.)
    design=build(CableDesign,"Rallis-wire",Group(:core,Region(:metal,Disk(.01),copper)))
    system=build(LineCableSystem,[design,design],[(0.,10.),(3.,-10.)];
        connections=[Dict("core"=>1),Dict("core"=>2)])
    problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),
        frequencies=[50.,5000.],temperature=20.)
    selected=EI.Formula(:Pollaczek1926;evaluation=:dcim_two_level)
    result=compute(problem,Formulation(earth_impedance=selected,
        earth_admittance=:MartinsBritto2024,
        options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false));
        options=(trace=true,))
    trace=details(result).trace
    for k in 1:2
        s=2π*problem.frequencies[k]*im
        leaf=selected([Inf,100.],[ε0,10ε0],[μ0,μ0],s,nothing)
        pair=E.EarthPair(1,2,(10.,-10.),3.,(1,2))
        @test trace.Zg[1,2,k] ≈ leaf(Val(:mutual),pair) rtol=1e-10
        @test result.Z[:,:,k] ≈ transpose(result.Z[:,:,k])
        @test result.Y[:,:,k] ≈ transpose(result.Y[:,:,k])
        @test minimum(eigvals(Symmetric(real.(result.Z[:,:,k]))))>=-1e-12
    end
end
