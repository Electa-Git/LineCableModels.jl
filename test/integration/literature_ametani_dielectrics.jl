@testitem "Engine / literature assimilation / Ametani annular magnetic and dielectric terms" begin
    using LineCableModels,LinearAlgebra
    E=LineCableModels.Engine;IZ=E.InsulationImpedance
    IA=E.InsulationAdmittance;SA=E.SemiconAdmittance
    for T in (Float32,Float64,BigFloat),frequency in (T(0),T(50),T(1e6)),
            a in (T(1e-8),T(.01)),ratio in (T(1.001),T(2)),mur in (T(1),T(5))
        mu0=T(4)*T(π)/T(10)^7;eps0=T(88541878128)/T(10)^22
        s=complex(zero(T),2T(π)*frequency);b=a*ratio
        z=IZ.Formula(:Ametani1980)(a,b,mur,s)
        expected=s*mu0*mur/(2T(π))*log(b/a)
        @test z isa Complex{T}
        @test z ≈ expected rtol=500eps(T)
        @test IZ.Formula(:Ametani1980)(a,a,mur,s)==0
        for rho in (T(100),T(Inf))
            material=Material(:semicon,rho,T(4),one(T),T(20),zero(T))
            # The closed constitutive limit exists at DC; the public
            # frequency-domain matrix API still requires positive f.
            first=iszero(frequency) ? IA.insulation_material(
                Val(:Ametani2004),material,frequency,T(20),(;)) :
                IA.Formula(:Ametani2004)(material,frequency,T(20))
            second=iszero(frequency) ? SA.semicon_material(
                Val(:Ametani2004),material,frequency,T(20),(;)) :
                SA.Formula(:Ametani2004)(material,frequency,T(20))
            @test first==second
            @test first ≈ inv(rho)+s*eps0*material.eps_r
            @test E.layer_admittance(a,b,second) ≈
                2T(π)*(inv(rho)+s*eps0*material.eps_r)/log(b/a)
        end
    end
    copper=Material(:conductor,1.724e-8,1.,1.,20.,0.)
    function solve(mur)
        ins=Material(:insulator,Inf,2.3,mur,20.,0.)
        design=build(CableDesign,"magnetic-annuli",Stack(
            Group(:core,Region(:core,Disk(.01),copper)),
            Region(:inner,Shell(.002),ins),
            Group(:sheath,Region(:sheath,Shell(.001),copper)),
            Region(:outer,Shell(.003),ins)))
        system=build(LineCableSystem,[design],[(0.,-1.)];
            connections=[Dict("core"=>1,"sheath"=>2)])
        problem=LineParametersProblem(system;earth_props=EarthModel(100.,10.,1.),
            frequencies=[50.,1e5],temperature=20.)
        result=compute(problem,Formulation(
            options=(reduce_bundle=false,kron_reduction=false,ideal_transposition=false)))
        return result
    end
    reference=solve(1.);permeable=solve(5.)
    for (k,f) in enumerate((50.,1e5))
        factor=2π*f*im*4π*1e-7/(2π)*(5-1)
        z1=factor*log(.012/.01);z2=factor*log(.016/.013)
        @test permeable.Z[:,:,k]-reference.Z[:,:,k] ≈ [z1+z2 z2;z2 z2] rtol=2e-11
        @test permeable.Y[:,:,k] ≈ reference.Y[:,:,k] rtol=2e-11
    end
end
