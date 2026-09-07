@testitem "Engine / literature assimilation / Wait circuit reductions" begin
    using LineCableModels, QuadGK, SpecialFunctions
    E=LineCableModels.Engine; EI=E.EarthImpedance; EA=E.EarthAdmittance
    μ0=4π*1e-7; ε0=8.8541878128e-12
    for frequency in (1.0,50.0,10000.0,1e6), rho in (10.0,1000.0),
        (radius,height) in ((0.005,1.0),(0.05,20.0))
        s=complex(0.0,2π*frequency)
        args=([Inf,rho],[ε0,10ε0],[μ0,μ0],s,nothing)
        buried=E.EarthPair(1,1,(-height,-height),radius,(2,2))
        # Direct source equation (12), at higher precision so the divergent
        # K1/alpha and exponential terms can cancel without accuracy loss.
        literal,smallliteral=setprecision(256) do
            sb=Complex{BigFloat}(s); mub=BigFloat(μ0); eb=BigFloat(ε0)
            k=sqrt(-sb*mub*(inv(BigFloat(rho))+sb*10eb))
            α=2im*k*BigFloat(height)
            image=E.special_besselk(0,α)+2E.special_besselk(1,α)/α-
                2*(1+α)*exp(-α)/α^2
            direct=E.special_besselk(0,im*k*BigFloat(radius))
            z=sb*mub/(2BigFloat(π))*(direct+image)
            zsmall=-sb*mub/(2BigFloat(π))*(log(BigFloat(0.89)*k*BigFloat(radius))+
                im*BigFloat(π)/2-image)
            ComplexF64(z),ComplexF64(zsmall)
        end
        f=EI.Formula(:Wait1978)(args...)
        fs=EI.Formula(:Wait1978;approximation=:small_argument)(args...)
        @test f(Val(:self),buried) ≈ literal rtol=1e-10
        @test fs(Val(:self),buried) ≈ smallliteral rtol=1e-10
        @test real(f(Val(:self),buried))>=0
        @test_throws ArgumentError f(Val(:mutual),buried)
        overhead=E.EarthPair(1,1,(height,height),radius,(1,1))
        potential=EA.Formula(:Wait1972a)(args...)
        ysource=s*2π*ε0/log(2height/radius)
        @test s/potential(Val(:self),overhead) ≈ ysource rtol=1e-12
        @test_throws ArgumentError potential(Val(:mutual),overhead)
        # Wait (31),(33), conduction-only Carson qTEM reduction.
        k2=-s*μ0/rho
        jc=2/k2*quadgk(0.0,Inf;rtol=1e-11) do λ
            # Rationalization of u-lambda, not a change of the source root.
            (-k2)/(sqrt(λ^2-k2)+λ)*exp(-2height*λ)
        end[1]
        expected=s*μ0/(2π)*(log(2height/radius)-jc)
        @test EI.Formula(:Carson1926)(args...)(Val(:self),overhead) ≈ expected rtol=3e-8
    end
    copper=Material(:conductor,1.7241e-8,1.0,1.0,20.0,0.0)
    design=build(CableDesign,"Wait-wire",Group(:core,Region(:core,Disk(0.01),copper)))
    system=build(LineCableSystem,[design],[(0.0,10.0)];connections=[Dict("core"=>1)])
    EP=LineCableModels.EarthProps
    earth=build(EP.EarthModel,(EP.EarthLayer(100.0,10.0,1.0),))
    problem=LineParametersProblem(system;earth_props=earth,frequencies=[50.0,1e5])
    formulation=Formulation(earth_impedance=:Carson1926,earth_admittance=:Wait1972a,
        options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false))
    result=compute(problem,formulation)
    @test all(isfinite,result.Z)
    @test all(isfinite,result.Y)
    @test size(result.Y)==(1,1,2)
end
