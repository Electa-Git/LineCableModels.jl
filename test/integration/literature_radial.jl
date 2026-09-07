@testitem "Engine / literature assimilation / radial admittance sources" begin
    using LineCableModels, LinearAlgebra
    const E=LineCableModels.Engine
    const IA=E.InsulationAdmittance
    const SA=E.SemiconAdmittance
    @test formula_id(IA.Formula(:Gustavsen2013)) === :Ametani1980
    for T in (Float32,Float64,BigFloat)
        ε0=T(88541878128)*T(10)^(-22)
        tolerance=T===Float32 ? T(2e-5) : T(1e-12)
        for frequency in T.((1,50,100000))
            s=complex(zero(T),T(2)*T(π)*frequency)
            material=Material(:insulator,T(200),T(3),one(T),T(20),zero(T))
            a,b=T(1)/100,T(2)/100
            for tangent in T.((0,0.01,0.2))
                law=IA.Formula(:Weeks1984;loss_tangent=tangent)
                y=E.layer_admittance(a,b,law(material,frequency,T(20)))
                C=T(2)*T(π)*ε0*material.eps_r/log(b/a)
                expected=s*C*(one(T)-complex(zero(T),tangent))
                @test y ≈ expected rtol=tolerance
                @test real(y)>=zero(T)
                @test y isa Complex{T}
            end
            lossless=IA.Formula(:Ametani1980)(material,frequency,T(20))
            @test iszero(real(lossless))
            @test E.potential_coefficient(a,b,lossless,s) ≈
                log(b/a)/(T(2)*T(π)*ε0*material.eps_r) rtol=tolerance

            # Pawlik (13): full conductivity, not a constant loss tangent.
            coated=IA.Formula(:Ametani2004)(material,frequency,T(20))
            yi=T(2)*T(π)*(inv(material.rho)+s*ε0*material.eps_r)/log(b/a)
            pi=E.potential_coefficient(a,b,coated,s)
            @test s/pi ≈ yi rtol=tolerance
            ye=complex(T(1e-7),T(3e-6))*frequency
            @test s/(pi+s/ye) ≈ inv(inv(yi)+inv(ye)) rtol=tolerance

            # Weeks' three-annulus sum and Ghosh's N-screen network use the
            # same algebra once each physical layer admittance is specified.
            # This does not endorse the missing 2π / relative-permittivity
            # normalization in the printed Ghosh layer coefficients.
            radii=T.((10,11,15,16,18,19))./T(1000)
            materials=[Material(k,T(r),T(e),one(T),T(20),zero(T))
                for (k,r,e) in ((:semicon,1000,1000),(:insulator,1e12,2.3),
                    (:semicon,500,500),(:insulator,1e11,3.0),(:semicon,800,800))]
            p=Complex{T}[]
            y=Complex{T}[]
            for n in eachindex(materials)
                m=materials[n]
                κ=(m.kind===:semicon ? SA.Formula(:Ametani2004) :
                    IA.Formula(:Ametani2004))(m,frequency,T(20))
                push!(p,E.potential_coefficient(radii[n],radii[n+1],κ,s))
                push!(y,T(2)*T(π)*κ/log(radii[n+1]/radii[n]))
            end
            @test s/E.radial_coefficient(p,1:3) ≈
                inv(inv(y[1])+inv(y[2])+inv(y[3])) rtol=tolerance
            @test s/E.radial_coefficient(p,4:5) ≈ inv(inv(y[4])+inv(y[5])) rtol=tolerance
            @test s/E.radial_coefficient(p,5:5) ≈ y[5] rtol=tolerance
            @test s/E.radial_coefficient(p,1:5) ≈ inv(sum(inv,y)) rtol=tolerance
        end
    end
    material=Material(:insulator,Inf,2.3,1.0,20.0,0.0)
    for invalid in (-0.1,Inf,NaN,"bad")
        law=IA.Formula(:Weeks1984;loss_tangent=invalid)
        @test_throws DomainError law(material,50.0,20.0)
    end

    copper=Material(:conductor,1.724e-8,1.0,1.0,20.0,0.0)
    dielectric=Material(:insulator,Inf,2.3,1.0,20.0,0.0)
    design=build(CableDesign,"radial-source-matrix",Stack(
        Group(:core,Region(:core,Disk(0.01),copper)),
        Region(:inner,Shell(0.002),dielectric),
        Group(:sheath,Region(:sheath,Shell(0.001),copper)),
        Region(:middle,Shell(0.003),dielectric),
        Group(:armor,Region(:armor,Shell(0.001),copper)),
        Region(:outer,Shell(0.002),dielectric)))
    blueprint=E.flatten(LineCableModelsCoaxial(),design)
    input=E.LocalCableData(blueprint)
    s=complex(0.0,100π)
    methods=Formulation(insulation_admittance=:Ametani1980).methods
    layers=zeros(ComplexF64,3); terms=similar(layers); tails=similar(layers)
    P=zeros(ComplexF64,3,3)
    E.cable_potential!(P,input,methods,50.0,20.0,s,layers,terms,tails)
    p=[log(layer.r_ex/layer.r_in)/(2π*8.8541878128e-12*layer.material.eps_r)
        for layer in blueprint.dielectrics]
    source=[sum(p) p[2]+p[3] p[3]; p[2]+p[3] p[2]+p[3] p[3]; p[3] p[3] p[3]]
    @test P ≈ source rtol=1e-12
    Y=zeros(ComplexF64,3,3)
    E.cable_admittance!(Y,input,methods,50.0,20.0,s,layers)
    @test Y ≈ s*inv(source) rtol=1e-12
end
