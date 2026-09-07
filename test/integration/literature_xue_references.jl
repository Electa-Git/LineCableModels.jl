@testitem "Engine / literature assimilation / Xue and Magalhaes voltage references" begin
    using LineCableModels,SpecialFunctions,QuadGK
    E=LineCableModels.Engine;EI=E.EarthImpedance;EA=E.EarthAdmittance
    mu0=4π*1e-7;eps0=8.8541878128e-12
    root(z,s)=iszero(imag(z))&&real(z)<0 ?
        complex(0.,sign(imag(s))*sqrt(-real(z))) : sqrt(z)
    for frequency in (50.,1e5,1e7),airrho in (Inf,1e7)
        s=2π*frequency*im;rho=[airrho,100.];epsilon=eps0*[1.,10.];mu=[mu0,mu0]
        g=s.*mu.*(inv.(rho).+s.*epsilon);gamma=sqrt(g[2])
        cut=iszero(imag(g[1]))&&real(g[1])<0 ? sqrt(-real(g[1])) : 0.
        bounds=cut>0 ? (0.,cut,Inf) : (0.,Inf)
        for pair in (E.EarthPair(1,1,(-.5,-.5),.02,(2,2)),
                E.EarthPair(1,2,(-.5,-1.5),.7,(2,2)))
            H=sum(abs,pair.heights);x=pair.separation
            d=hypot(x,pair.heights[1]-pair.heights[2]);D=hypot(x,H)
            function integrals(height)
                quadgk(bounds...;rtol=2e-10) do u
                    a0=root(u^2+g[1],s);a1=root(u^2+g[2],s)
                    weight=exp(-height*a1)*cos(x*u)/a1^2
                    # Separate S11, S12, S13 as printed; no identity
                    # contraction from either engine evaluator is used.
                    weight*[u^2/(a0+a1),u^2/(a0+g[1]/g[2]*a1),1/(a0+a1)]
                end[1]
            end
            q=integrals(H)
            B=besselk(0,gamma*d)-besselk(0,gamma*D)
            expectedZ=s*mu0/(2π)*(B+2q[1]+2g[2]*q[3])
            factor=s/(2π*(1/rho[2]+s*epsilon[2]))
            infinite=factor*(B+2q[2]+2g[2]*q[3])
            surface_terms=integrals(H/2)
            surface=infinite-2factor*(surface_terms[2]+g[2]*surface_terms[3])
            # Use the literal printed radical in an ordinary lossy range.
            omega=imag(s)
            delta=inv(sqrt(omega^2*epsilon[2]*mu0/2*
                (sqrt(1+(1/rho[2]/(omega*epsilon[2]))^2)-1)))
            penetration_terms=integrals(H/2+delta)
            Bd=besselk(0,gamma*hypot(x,-delta-H/2))-
                besselk(0,gamma*hypot(x,-delta+H/2))
            penetration=infinite-factor*(Bd+2penetration_terms[2]+2g[2]*penetration_terms[3])
            fi=EI.Formula(:Xue2018b)(rho,epsilon,mu,s,nothing)
            fa=EA.Formula(:Xue2018b)(rho,epsilon,mu,s,nothing)
            fm=EA.Formula(:Magalhaes2018)(rho,epsilon,mu,s,nothing)
            kind=pair.row==pair.column ? Val(:self) : Val(:mutual)
            @test fi(kind,pair) ≈ expectedZ rtol=2e-7
            @test fa(kind,pair) ≈ infinite rtol=2e-7
            @test fa.routes.surface(fa,pair) ≈ surface rtol=2e-7
            @test fa.routes.penetration(fa,pair) ≈ penetration rtol=2e-7
            @test fm(kind,pair) ≈ surface rtol=2e-7
            @test formula_id(EI.Formula(:Magalhaes2018))==:Xue2018b
            @test EI.Formula(:Magalhaes2018)(rho,epsilon,mu,s,nothing)(kind,pair)==fi(kind,pair)
            for (family,id) in ((EI,:Xue2018b),(EA,:Xue2018b),(EA,:Magalhaes2018))
                f=family.Formula(id)(rho,epsilon,mu,s,nothing)
                neg=family.Formula(id)(rho,epsilon,mu,-s,nothing)
                @test neg(kind,pair) ≈ conj(f(kind,pair)) rtol=2e-8
            end
        end
    end
    for T in (Float32,Float64,BigFloat),family in (EI,EA),id in
            (family==EI ? (:Xue2018b,) : (:Xue2018b,:Magalhaes2018))
        rho=T[Inf,100];epsilon=T(eps0)*T[1,10];mu=T(mu0)*ones(T,2)
        s=complex(zero(T),100T(π))
        pair=E.EarthPair(1,2,(T(-.5),T(-1.5)),T(.7),(2,2))
        leaf=family.Formula(id)(rho,epsilon,mu,s,nothing)
        value=leaf(Val(:mutual),pair)
        @test value isa Complex{T}
        @test isfinite(value)
    end
end
