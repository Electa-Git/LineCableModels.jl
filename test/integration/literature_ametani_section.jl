@testitem "Engine / literature assimilation / Ametani physical cross-section" begin
    using LineCableModels,LinearAlgebra
    E=LineCableModels.Engine;DM=LineCableModels.DataModel;II=E.InternalImpedance
    μ0=4π*1e-7
    for T in (Float32,Float64,BigFloat),angle in (.5,2.),f in (0.,50.,1e8)
        r=T(.01);θ=T(angle);S=θ*r^2/2;ell=r*(2+θ)
        rho=T(1.7241e-8);mur=T(3)
        s=complex(zero(T),T(2)*T(π)*T(f))
        leaf=II.Formula(:Ametani1992)(Val(:section),S,ell,rho,mur,s)
        reference=rho/S*sqrt(1+s*T(μ0)*mur*S/((rho/S)*ell^2))
        tolerance=T===Float32 ? 2e-6 : 1e-12
        @test leaf(Val(:outer)) ≈ reference rtol=tolerance
        @test leaf(Val(:outer)) isa Complex{T}
        @test iszero(leaf(Val(:inner))) && iszero(leaf(Val(:mutual)))
        @test T(π)*(leaf.state.r_ex^2-leaf.state.r_in^2) ≈ S rtol=tolerance
        @test 2T(π)*leaf.state.r_ex ≈ ell rtol=tolerance
        reversed=II.Formula(:Ametani1992)(Val(:section),S,ell,rho,mur,-s)
        @test reversed(Val(:outer)) ≈ conj(leaf(Val(:outer))) rtol=tolerance
    end
    @test_throws DomainError II.Formula(:Ametani1992)(
        Val(:section),1.,.1,1e-8,1.,100π*im)
    copper=Material(:conductor,1.7241e-8,1.,2.,20.,.003)
    for primitive in (Disk(.01),Annulus(.005,.01),
            Sector(2π/3,0.,.02),Sector(2π/3,.003,.02,.001)),
            temperature in (20.,80.),correction in (false,true)
        design=build(CableDesign,"physical-section",
            Group(:core,Region(:metal,primitive,copper)))
        shape=only(design.geometry.regions).primitive
        section=only(E.flatten(LineCableModelsCoaxial(),design).conductors).section
        @test section.area ≈ DM.area(shape)
        expected_perimeter=shape isa DM.Annulus ? 2π*shape.ro : DM.perimeter(shape)
        @test section.perimeter ≈ expected_perimeter
        @test section.material==copper
        system=build(LineCableSystem,[design,design],[(-1.,10.),(1.,10.)];
            connections=[Dict("core"=>1),Dict("core"=>2)])
        problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),
            frequencies=[50.,1e5],temperature)
        result=compute(problem,Formulation(internal_impedance=:Ametani1992,
            earth_impedance=:Carson1926,earth_admittance=:Wise1948,
            options=(temperature_correction=correction,reduce_bundle=false,
                kron_reduction=false,ideal_transposition=false));options=(trace=true,))
        trace=details(result).trace
        rho=copper.rho*(correction ? 1+copper.alpha*(temperature-20) : 1)
        for k in 1:2
            s=2π*problem.frequencies[k]*im
            value=II.Formula(:Ametani1992)(Val(:section),section.area,
                section.perimeter,rho,copper.mu_r,s)(Val(:outer))
            @test trace.Zin[:,:,k] ≈ value*Matrix{Float64}(I,2,2) rtol=1e-11
            @test result.Z[:,:,k] ≈ transpose(result.Z[:,:,k])
            @test result.Y[:,:,k] ≈ transpose(result.Y[:,:,k])
        end
    end
    # Sector contours must not collapse back to equal-area circular ones.
    S=2π/3*.02^2/2;ell=.02*(2+2π/3)
    physical=II.Formula(:Ametani1992)(Val(:section),S,ell,copper.rho,2.,2e5π*im)
    circular=II.Formula(:Ametani1992)(0.,sqrt(S/π),copper.rho,2.,2e5π*im)
    @test abs(physical(Val(:outer))/circular(Val(:outer))-1)>.1
end
